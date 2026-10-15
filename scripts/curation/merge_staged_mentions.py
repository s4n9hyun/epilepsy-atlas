#!/usr/bin/env python3
"""재추출 결과(스테이징 JSON)를 기존 DB 행과 대조해 빠진 환자만 추가한다.

기존 mention은 지우거나 덮지 않는다 (CP 고정·병합 규칙·정정 규칙이 mention_id에 걸려 있음).
1. `batch_extract_patient_mentions.py --stage-dir DIR`로 만든 {GENE}_{PMID}.json을 읽는다.
2. 같은 PMID의 기존 original/manual_exclude mention과 새 original 환자 행을 맞춘다.
   우선순위: 환자 ID+변이 > 변이가 같은 자리표시 행(환자 ID 없음, "patient with p.X" 등) > 환자 ID가 같고 한쪽 변이 없음
   > 변이 같고 한쪽 환자 ID 없음. 같은 순위면 현재 canonical에 있는 행(CP 번호 있음)을 우선한다.
   - matched / new / ambiguous(변이만 같고 환자 ID가 다름, 또는 후보 여럿)
   인용·요약 행(cited, aggregate_only 등)은 다른 환자이므로 대조하지 않는다.
3. 보고서 CSV: matched 행과의 값 차이.
4. 채우기 규칙(fill_rules_*.csv, config/manual_value_corrections.csv 형식):
   - matched 기존 행의 빈 칸(변이·성별·나이 등)과 자리표시 patient_id를 새 값으로
   - 새 행이 자리표시 행과 짝지어져 남게 된, 같은 환자 ID의 변이 없는 기존 행은 중복으로 manual_exclude
   기존에 값이 있는 칸은 바꾸지 않는다 (차이는 보고서에서 검토).
5. --apply일 때만 new 행 추가 + 채우기 규칙을 config에 추가하고 DB에 적용한다.

사용:
  python3 scripts/curation/merge_staged_mentions.py --stage-dir runs/X/staging --gene SCN1A            # 대조 보고만
  python3 scripts/curation/merge_staged_mentions.py --stage-dir runs/X/staging --gene SCN1A --apply
"""
import argparse
import csv
import json
import re
import sqlite3
import sys
from datetime import datetime
from pathlib import Path

WORKSPACE_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(WORKSPACE_ROOT / "scripts" / "core"))
from batch_extract_patient_mentions import PatientGeneDatabase, flatten_patients  # noqa: E402

MISSING = {"", "not specified", "n/a", "na", "none", "unknown", "nan"}
MATCHABLE_TYPES = {"original_patient", "manual_exclude", "unclear"}
NON_ORIGINAL = {"previously_reported_patient", "reanalyzed_previous_patient", "cited_patient"}
TABLE_OVERWRITE_FIELDS = ["sex", "age_of_onset", "last_followup_age", "seizure_types", "phenotype_review_category",
                          "response_to_asm", "epilepsy_type", "other_phenotypes", "variant_c_original", "variant_p"]
PLACEHOLDER_PID = re.compile(r"(patients?|individuals?|cases?|probands?|carriers?) with (p\.|c\.)|mutated patient|"
                             r"^(not specified|unknown)$|aggregate|cohort patient", re.I)
VARIANT_FIELDS = ["variant_c_original", "variant_c_gnomad", "variant_p", "nm_transcript", "nm_number", "gnomad_transcript",
                  "gnomad_conversion_status", "clinvar", "clinvar_pathogenicity", "clinvar_id"]
FILL_FIELDS = ["sex", "age_of_onset", "last_followup_age", "seizure_presence", "seizure_types", "phenotype_review_category",
               "response_to_asm", "asm_used", "responsive_asm", "targeted_treatment_response", "epilepsy_type",
               "other_phenotypes"]
COMPARE_FIELDS = ["sex", "age_of_onset", "last_followup_age", "seizure_types", "phenotype_review_category",
                  "response_to_asm", "targeted_treatment_response", "epilepsy_type", "variant_p", "mention_source_type"]
AA3 = {"A": "ala", "R": "arg", "N": "asn", "D": "asp", "C": "cys", "Q": "gln", "E": "glu", "G": "gly", "H": "his",
       "I": "ile", "L": "leu", "K": "lys", "M": "met", "F": "phe", "P": "pro", "S": "ser", "T": "thr", "W": "trp",
       "Y": "tyr", "V": "val", "X": "ter", "*": "ter"}


def clean(value):
    value = "" if value is None else str(value).strip()
    return "" if value.lower() in MISSING else value


def norm_pid(value):
    value = clean(value).lower()
    value = re.sub(r"\b(patient|pt|case|subject|individual|proband|no\.?|#)\b", " ", value)
    return re.sub(r"[\s_\-:.()]+", "", value)


def c_key(*values):
    for value in values:
        value = clean(value)
        if "c." in value:
            core = "c." + value.split("c.", 1)[1]
            core = re.sub(r"\s", "", core).lower()
            return re.sub(r"(del|dup)[acgt]+$", r"\1", core)
    return ""


def p_key(value):
    value = clean(value)
    m = re.match(r"p\.\(?([A-Za-z]{3}|[A-Z*])(\d+)", value)
    if not m:
        return ""
    ref = m.group(1).lower() if len(m.group(1)) == 3 else AA3.get(m.group(1).upper(), m.group(1).lower())
    rest = value[m.end():].lower()
    alt = "fs" if "fs" in rest else re.sub(r"[()\s]", "", rest)[:3]
    alt = AA3.get(alt.upper(), alt) if len(alt) == 1 else alt
    return f"p.{ref}{m.group(2)}{'ter' if alt in ('*', 'x') else alt}"


def variant_keys(c_values, p_value):
    return {k for k in (c_key(*c_values), p_key(p_value)) if k}


def is_placeholder(pid):
    return not clean(pid) or bool(PLACEHOLDER_PID.search(str(pid)))


def match_rows(new_rows, existing, canonical_mids, resolutions=None):
    """새 행마다 (status, 기존 mention_id, 이유)를 정한다. 기존 행 하나는 한 번만 짝지어진다.

    resolutions: {patient_id: ["same:123" | "different", ...]} 원문 대조로 정한 판정 (순서대로 소비)
    """
    used, out = set(), []
    resolutions = {k: list(v) for k, v in (resolutions or {}).items()}
    by_id = {str(e["mention_id"]): e for e in existing}
    # 판정으로 짝이 정해진 기존 행은 먼저 예약해 다른 새 행의 후보에서 뺀다
    reserved = {v.split(":", 1)[1].strip() for vals in resolutions.values() for v in vals if v.startswith("same:")}
    reserved_ids = {by_id[m]["mention_id"] for m in reserved if m in by_id}
    used |= reserved_ids
    claimed = set()
    # 환자 ID와 변이가 정확히 같은 짝도 먼저 예약한다 (처리 순서 때문에 모호해지지 않게)
    exact = {}
    for i, row in enumerate(new_rows):
        if resolutions.get(str(row.get("patient_id", "")).strip()):
            continue
        nkeys = variant_keys([row.get("variant_c_gnomad"), row.get("variant_c_original")], row.get("variant_p"))
        npid = "" if is_placeholder(row.get("patient_id")) else norm_pid(row.get("patient_id"))
        for ex in existing:
            if (npid and ex["mention_id"] not in used and ex["mention_id"] not in exact.values()
                    and norm_pid(ex["patient_id"]) == npid
                    and nkeys & variant_keys([ex["variant_c_gnomad"], ex["variant_c_original"]], ex["variant_p"])):
                exact[i] = ex["mention_id"]
                break
    used |= set(exact.values())
    for i, row in enumerate(new_rows):
        if i in exact:
            out.append(("matched", str(exact[i]), "patient_id+variant"))
            continue
        pending = resolutions.get(str(row.get("patient_id", "")).strip(), [])
        forced = pending.pop(0) if pending else ""
        if forced.startswith("same:"):
            mid = forced.split(":", 1)[1].strip()
            ex = by_id.get(mid)
            if ex is not None and ex["mention_id"] not in claimed:
                claimed.add(ex["mention_id"])
                out.append(("matched", mid, "verified same patient (resolution)"))
                continue
            if ex is not None:
                out.append(("new", "", f"additional variant row of patient already matched to {mid}"))
                continue
        nkeys = variant_keys([row.get("variant_c_gnomad"), row.get("variant_c_original")], row.get("variant_p"))
        npid = "" if is_placeholder(row.get("patient_id")) else norm_pid(row.get("patient_id"))
        cands = []
        for ex in existing:
            if ex["mention_id"] in used or (forced == "different" and str(ex["mention_id"]) in forced_block(row, existing)):
                continue
            ekeys = variant_keys([ex["variant_c_gnomad"], ex["variant_c_original"]], ex["variant_p"])
            epid = "" if is_placeholder(ex["patient_id"]) else norm_pid(ex["patient_id"])
            same_var = bool(nkeys & ekeys)
            same_pid = bool(npid) and npid == epid
            sex_conflict = clean(row.get("sex")) and clean(ex["sex"]) and clean(row.get("sex")) != clean(ex["sex"])
            if same_pid and same_var:
                rank, why = 0, "patient_id+variant"
            elif same_var and not epid and not sex_conflict:
                rank, why = 1, "placeholder row with same variant"
            elif same_pid and (not nkeys or not ekeys):
                rank, why = 2, "patient_id, one side without variant"
            elif same_var and not npid and not sex_conflict:
                rank, why = 3, "variant, new row without patient_id"
            elif same_var and not sex_conflict:
                rank, why = 4, "variant only (patient_id differs)"
            else:
                continue
            # 같은 순위면 제외되지 않은 행, 그다음 canonical에 있는 행을 우선
            tie = (1 if ex["mention_source_type"] == "manual_exclude" else 0) * 2 + (0 if str(ex["mention_id"]) in canonical_mids else 1)
            cands.append((rank, tie, ex, why))
        cands.sort(key=lambda c: (c[0], c[1]))
        best = [c for c in cands if c[:2] == cands[0][:2]] if cands else []
        if not best:
            out.append(("new", "", "no existing row with same patient_id or variant"))
        elif best[0][0] == 4:
            # 변이만 같고 환자 ID가 다름: 재발 변이를 가진 다른 환자일 수 있어 사람이 본다
            out.append(("ambiguous", ";".join(str(c[2]["mention_id"]) for c in best), best[0][3]))
        elif len(best) > 1:
            out.append(("ambiguous", ";".join(str(c[2]["mention_id"]) for c in best), best[0][3] + ", multiple candidates"))
        else:
            used.add(best[0][2]["mention_id"])
            out.append(("matched", str(best[0][2]["mention_id"]), best[0][3]))
    return out, used


def forced_block(row, existing):
    """'different' 판정 행은 같은 변이를 가진 기존 행과 짝짓지 않는다."""
    nkeys = variant_keys([row.get("variant_c_gnomad"), row.get("variant_c_original")], row.get("variant_p"))
    return {str(e["mention_id"]) for e in existing
            if nkeys & variant_keys([e["variant_c_gnomad"], e["variant_c_original"]], e["variant_p"])}


def canonical_mention_ids(gene):
    import glob
    import pandas as pd
    files = glob.glob(str(WORKSPACE_ROOT / "outputs" / "current" / "canonical_patients" / f"{gene}_canonical_patients" / "*.csv"))
    mids = set()
    for f in files:
        for value in pd.read_csv(f, dtype=str)["supporting_mention_ids"].fillna(""):
            mids |= {m.strip() for m in re.split(r"[;,]", value) if m.strip()}
    return mids


def register_article(db_path, gene, pmid, stage_file, n_added):
    """논문 테이블에 없는 PMID면 등록한다 (다른 유전자 DB의 같은 논문 정보를 복사)."""
    con = sqlite3.connect(db_path)
    if con.execute("SELECT 1 FROM pubmed_articles WHERE pmid=?", (pmid,)).fetchone():
        con.close()
        return
    now = datetime.now().isoformat(timespec="seconds")
    base = {}
    for other in sorted(Path(db_path).parent.glob("*.sqlite")):
        o = sqlite3.connect(other)
        o.row_factory = sqlite3.Row
        try:
            r = o.execute("SELECT * FROM pubmed_articles WHERE pmid=?", (pmid,)).fetchone()
        except sqlite3.Error:
            r = None
        o.close()
        if r:
            base = dict(r)
            break
    row = {**base, "pmid": pmid, "gene": gene, "source_csv": "merge_staged_mentions.py (article added by re-extraction)",
           "imported_at": now, "extracted": "Yes", "extraction_status": "success", "extracted_at": now}
    cols = [c[1] for c in con.execute("PRAGMA table_info(pubmed_articles)")]
    con.execute(f"INSERT INTO pubmed_articles ({','.join(cols)}) VALUES ({','.join('?' * len(cols))})",
                [row.get(c, "") for c in cols])
    pcols = [c[1] for c in con.execute("PRAGMA table_info(processed_articles)")]
    prow = {"article_id": pmid, "pmcid": base.get("pmcid", ""), "pmid": pmid, "gene": gene, "status": "success",
            "patient_count": n_added, "used_supplementary": "True", "files_processed": json.dumps([stage_file]),
            "error": "", "processed_at": now}
    con.execute(f"INSERT OR IGNORE INTO processed_articles ({','.join(pcols)}) VALUES ({','.join('?' * len(pcols))})",
                [prow.get(c, "") for c in pcols])
    con.commit()
    con.close()


def fill_rule(gene, ex, field, value, evidence):
    return {"gene": gene, "mention_id": str(ex["mention_id"]), "pmid": str(ex["pmid"]), "patient_id": ex["patient_id"],
            "field": field, "value": str(value), "previous_value": ex.get(field) or "", "evidence": evidence,
            "source_run": "", "status": "active"}


def main():
    parser = argparse.ArgumentParser(description="재추출 스테이징 결과를 기존 DB와 대조·병합")
    parser.add_argument("--stage-dir", required=True)
    parser.add_argument("--gene", nargs="+", required=True)
    parser.add_argument("--pmid", nargs="*", help="특정 PMID만")
    parser.add_argument("--apply", action="store_true", help="new 행을 DB에 추가 (없으면 보고만)")
    parser.add_argument("--db-dir", default=str(WORKSPACE_ROOT / "databases"))
    parser.add_argument("--resolutions", help="원문 대조 판정 CSV (gene,pmid,patient_id,resolution=same:<mention_id>|different)")
    args = parser.parse_args()
    resolutions = {}
    if args.resolutions:
        for r in csv.DictReader(open(args.resolutions, encoding="utf-8")):
            if r["patient_id"].startswith("("):
                continue
            resolutions.setdefault((r["gene"].upper(), r["pmid"]), {}).setdefault(r["patient_id"].strip(), []).append(r["resolution"])
    corrected = {(r["gene"].upper(), r["mention_id"], r["field"]) for r in
                 csv.DictReader(open(WORKSPACE_ROOT / "config" / "manual_value_corrections.csv", encoding="utf-8-sig"))
                 if r.get("status", "active") == "active"}

    stage_dir = Path(args.stage_dir)
    if not stage_dir.is_absolute():
        stage_dir = WORKSPACE_ROOT / stage_dir
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    report_path = stage_dir / f"merge_report_{stamp}.csv"
    added_path = stage_dir / f"added_mentions_{stamp}.csv"
    report, added, fills, canonical_mids = [], [], [], {}
    for gene in [g.upper() for g in args.gene]:
        db_path = Path(args.db_dir) / f"{gene}.sqlite"
        for path in sorted(stage_dir.glob(f"{gene}_*.json")):
            result = json.loads(path.read_text(encoding="utf-8"))
            pmid = str(result.get("pmid", ""))
            if args.pmid and pmid not in args.pmid:
                continue
            if result.get("status") != "success":
                report.append({"gene": gene, "pmid": pmid, "status": f"stage_{result.get('status')}",
                               "reason": str(result.get("error", ""))[:200]})
                continue
            conn = sqlite3.connect(db_path)
            conn.row_factory = sqlite3.Row
            existing = [dict(r) for r in conn.execute(
                "SELECT * FROM patient_mentions WHERE pmid=?", (pmid,)) if r["mention_source_type"] in MATCHABLE_TYPES]
            db_columns = {r[1] for r in conn.execute("PRAGMA table_info(patient_mentions)")}
            conn.close()
            all_rows = flatten_patients(result)
            # 인용·요약 행은 다른 환자이거나 환자가 아니므로 대조·추가 대상이 아니다 (보고서에만 남김)
            new_rows = [r for r in all_rows if clean(r.get("mention_source_type")) == "original_patient"]
            for r in all_rows:
                if r not in new_rows:
                    report.append({"gene": gene, "pmid": pmid, "status": f"skipped_{clean(r.get('mention_source_type')) or 'unknown'}",
                                   "patient_id": r.get("patient_id", ""), "variant_c": r.get("variant_c_original", ""),
                                   "variant_p": r.get("variant_p", "")})
            decisions, used = match_rows(new_rows, existing, canonical_mids.setdefault(gene, canonical_mention_ids(gene)),
                                         resolutions.get((gene, pmid)))
            # 원문상 이미 발표된 환자인데 기존 DB에 original로 들어간 행 → 분류 변경
            for r in all_rows:
                src = clean(r.get("mention_source_type"))
                rpid = "" if is_placeholder(r.get("patient_id")) else norm_pid(r.get("patient_id"))
                if src not in NON_ORIGINAL or not rpid:
                    continue
                rkeys = variant_keys([r.get("variant_c_gnomad"), r.get("variant_c_original")], r.get("variant_p"))
                for ex in existing:
                    ekeys = variant_keys([ex["variant_c_gnomad"], ex["variant_c_original"]], ex["variant_p"])
                    if (ex["mention_source_type"] in ("original_patient", "unclear") and norm_pid(ex["patient_id"]) == rpid
                            and (not ekeys or rkeys & ekeys)):
                        fills.append(fill_rule(gene, ex, "mention_source_type", src,
                                               f"PMID {pmid} {ex['patient_id']}: {r.get('originality_evidence', '')} ({path.name})"))
            by_id = {str(e["mention_id"]): e for e in existing}
            to_add = []
            for i, (row, (status, mid, reason)) in enumerate(zip(new_rows, decisions)):
                ex = by_id.get(mid) if status == "matched" else None
                diffs = []
                if ex:
                    for field in COMPARE_FIELDS:
                        old, new = clean(ex.get(field)), clean(row.get(field))
                        if new and old != new:
                            diffs.append(f"{field}: {old or '(empty)'} -> {new}")
                if ex:
                    ev = f"PMID {pmid} {row.get('patient_id', '')}: re-extraction {path.name} ({row.get('location_in_paper', '')})"
                    if is_placeholder(ex["patient_id"]) and not is_placeholder(row.get("patient_id")):
                        fills.append(fill_rule(gene, ex, "patient_id", row.get("patient_id"), ev))
                    if ex["mention_source_type"] == "unclear":
                        fills.append(fill_rule(gene, ex, "mention_source_type", "original_patient", ev))
                    from_table = json.loads(row.get("raw_patient_json") or "{}").get("parsed_from_table")
                    if from_table:
                        # 표를 직접 읽은 값은 원문 그대로이므로 기존 값과 다르면 고친다 (사람이 정정한 칸은 제외)
                        for field in TABLE_OVERWRITE_FIELDS:
                            new = clean(row.get(field))
                            if (field in db_columns and new and clean(ex.get(field)) and clean(ex.get(field)) != new
                                    and (gene, str(ex["mention_id"]), field) not in corrected
                                    and not (field in ("variant_c_original", "variant_p")
                                             and variant_keys([row.get("variant_c_gnomad"), row.get("variant_c_original")], row.get("variant_p"))
                                             & variant_keys([ex["variant_c_gnomad"], ex["variant_c_original"]], ex["variant_p"]))):
                                fills.append(fill_rule(gene, ex, field, row.get(field), ev + " [table value]"))
                    # 짝지어진 행은 같은 변이이거나 기존 쪽 변이가 없으므로, 빈 변이 칸은 항목별로 채운다
                    for field in VARIANT_FIELDS + FILL_FIELDS:
                        if field in db_columns and not clean(ex.get(field)) and clean(row.get(field)):
                            fills.append(fill_rule(gene, ex, field, row.get(field), ev))
                report.append({
                    "gene": gene, "pmid": pmid, "staged_row": i + 1, "status": status, "matched_mention_id": mid,
                    "existing_source_type": ex["mention_source_type"] if ex else "",
                    "reason": reason, "patient_id": row.get("patient_id", ""),
                    "variant_c": row.get("variant_c_gnomad") or row.get("variant_c_original", ""),
                    "variant_p": row.get("variant_p", ""), "sex": row.get("sex", ""),
                    "age_of_onset": row.get("age_of_onset", ""), "location_in_paper": row.get("location_in_paper", ""),
                    "differences_vs_existing": " | ".join(diffs),
                })
                if status == "new":
                    to_add.append(row)
            # 자리표시 행과 짝지어진 새 환자의, 같은 환자 ID를 가진 변이 없는 기존 행 → 중복 제외
            for row, (status, mid, reason) in zip(new_rows, decisions):
                if status != "matched" or reason != "placeholder row with same variant":
                    continue
                npid = norm_pid(row.get("patient_id"))
                for ex in existing:
                    if (ex["mention_id"] not in used and ex["mention_source_type"] == "original_patient" and npid
                            and norm_pid(ex["patient_id"]) == npid
                            and not variant_keys([ex["variant_c_gnomad"], ex["variant_c_original"]], ex["variant_p"])):
                        fills.append(fill_rule(gene, ex, "mention_source_type", "manual_exclude",
                                               f"PMID {pmid} {ex['patient_id']}: duplicate of mention {mid} (placeholder row "
                                               f"now carrying this patient's re-extracted data; {path.name})"))
                        used.add(ex["mention_id"])
            if args.apply and to_add:
                register_article(db_path, gene, pmid, path.name, len(to_add))
                db = PatientGeneDatabase(db_path)
                with db.conn:
                    for row in to_add:
                        raw = json.loads(row.get("raw_patient_json") or "{}")
                        raw["reextract_merge"] = {"stage_file": path.name, "merged_at": stamp}
                        row["raw_patient_json"] = json.dumps(raw, ensure_ascii=False)
                        db._insert_mention(row)
                        mid = db.conn.execute("SELECT last_insert_rowid()").fetchone()[0]
                        added.append({"action": "add", "mention_id": str(mid), "field": "(row)", "old": "", "new": "added",
                                      "evidence": f"PMID {pmid} {row.get('patient_id', '')}: re-extraction ({path.name}), "
                                                  f"{row.get('location_in_paper', '')}; not in DB before"})
                db.close()
            print(f"{gene} {pmid}: staged {len(new_rows)}, matched {sum(d[0] == 'matched' for d in decisions)}, "
                  f"new {len(to_add)}, ambiguous {sum(d[0] == 'ambiguous' for d in decisions)}"
                  f"{' -> added' if args.apply and to_add else ''}")
    if report:
        fields = sorted({k for r in report for k in r}, key=lambda k: list(report[0]).index(k) if k in report[0] else 99)
        with report_path.open("w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=fields)
            w.writeheader()
            w.writerows(report)
        print(f"report -> {report_path}")
    if fills:
        run_rel = str(stage_dir.relative_to(WORKSPACE_ROOT)) if stage_dir.is_relative_to(WORKSPACE_ROOT) else str(stage_dir)
        for r in fills:
            r["source_run"] = run_rel
        fill_path = stage_dir / f"fill_rules_{stamp}.csv"
        with fill_path.open("w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=list(fills[0]))
            w.writeheader()
            w.writerows(fills)
        print(f"fill/exclude rules {len(fills)} -> {fill_path}")
        production_db = Path(args.db_dir).resolve() == (WORKSPACE_ROOT / "databases").resolve()
        if args.apply and not production_db:
            # 시험용 DB 사본: config는 건드리지 않고 사본 DB에만 직접 적용
            for r in fills:
                conn = sqlite3.connect(Path(args.db_dir) / f"{r['gene']}.sqlite")
                conn.execute(f"UPDATE patient_mentions SET {r['field']}=? WHERE mention_id=?", (r["value"], int(r["mention_id"])))
                conn.commit()
                conn.close()
            print(f"   --db-dir가 운영 DB가 아님: 채우기 규칙 {len(fills)}건을 사본 DB에만 적용 (config 미변경)")
        if args.apply and production_db:
            with (WORKSPACE_ROOT / "config" / "manual_value_corrections.csv").open("a", newline="", encoding="utf-8") as f:
                csv.DictWriter(f, fieldnames=list(fills[0])).writerows(fills)
            import subprocess
            for gene in sorted({r["gene"] for r in fills}):
                subprocess.run([sys.executable, str(WORKSPACE_ROOT / "scripts" / "curation" / "apply_manual_corrections.py"),
                                "--gene", gene], check=True)
            for r in fills:
                added.append({"action": "exclude" if r["value"] == "manual_exclude" else "value_fix",
                              "mention_id": r["mention_id"], "field": r["field"], "old": r["previous_value"],
                              "new": r["value"], "evidence": r["evidence"]})
    if added:
        with added_path.open("w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=list(added[0]))
            w.writeheader()
            w.writerows(added)
        print(f"added {len(added)} -> {added_path} (rebuild_canonical.py --evidence-log에 사용)")


if __name__ == "__main__":
    main()
