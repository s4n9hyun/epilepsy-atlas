#!/usr/bin/env python3
"""canonical 수정 전/후를 비교해 교수님 CRF 통합본 기준의 행 단위 변경 기록을 만든다.

한 batch(예: 20260923_KCNQ2_feedback)마다 실행하며, 결과는 누적 CSV
(기본: outputs/current/curation_change_log.csv)에 추가된다. 같은 batch_id를 다시
실행하면 그 batch·유전자의 기존 행을 지우고 새로 쓴다 (여러 유전자가 같은 batch를 공유 가능).

변경 유형:
- value_corrected: 같은 행의 CRF 필드 값이 바뀜
- merged_duplicate: 중복 환자로 병합되어 이 CP가 없어지고 다른 CP로 합쳐짐
- row_excluded: canonical에서 제외됨
- db_only_value_corrected: DB 값은 고쳤으나 병합 등으로 canonical 화면에는 다른 행 값이 표시됨
- row_added: 수정 후 canonical에 새로 들어온 행 (복원·재분류 등)
- same_variant_row_collapsed: 같은 환자의 같은 변이가 논문마다 다른 표기로 두 행이던 것을 한 행으로 정리

근거(evidence)는 config/manual_merge_rules.csv, DB의 manual_exclude 사유,
그리고 --evidence-log로 준 적용 로그(action, mention_id, field, evidence)에서 가져온다.
"""

import argparse
import csv
import re
import sqlite3
from pathlib import Path

import pandas as pd

WORKSPACE_ROOT = Path(__file__).resolve().parents[2]
MERGE_RULES = WORKSPACE_ROOT / "config" / "manual_merge_rules.csv"
VARIANT_EXCLUSIONS = WORKSPACE_ROOT / "config" / "manual_variant_exclusions.csv"
CORRECTIONS = WORKSPACE_ROOT / "config" / "manual_value_corrections.csv"
DEFAULT_CRF = WORKSPACE_ROOT / "feedback" / "20260923" / "유전자별_Phenotype_CRF_통합_0904.xlsx"
DEFAULT_LOG = WORKSPACE_ROOT / "outputs" / "current" / "curation_change_log.csv"

# canonical 컬럼 -> CRF 컬럼
FIELDS = {
    "variant_c": "Variant cDNA change",
    "variant_p": "Variant protein change",
    "clinvar_pathogenicity": "ClinVar pathogenicity",
    "ClinVar ID": "ClinVar ID",
    "sex": "Sex",
    "age_of_onset": "Age at seizure onset",
    "last_followup_age": "Age at last follow-up",
    "phenotype_review_category": "Phenotype",
    "targeted_treatment_response": "Response to targeted ASM",
}
# 적용 로그의 DB 필드명 -> canonical 필드명
LOG_FIELD_TO_CANONICAL = {"variant_c_gnomad": "variant_c", "variant_c_original": "variant_c", "variant_p": "variant_p",
                          "age_of_onset": "age_of_onset", "last_followup_age": "last_followup_age", "sex": "sex",
                          "targeted_treatment_response": "targeted_treatment_response",
                          "phenotype_review_category": "phenotype_review_category"}
COLUMNS = [
    "change_id", "change_date", "batch_id", "gene", "crf_patient_id", "current_patient_id",
    "pmid", "patient_in_paper", "variant_c", "variant_p", "change_type", "field",
    "crf_value", "before_value", "after_value", "reason", "evidence", "mention_id", "run_dir",
]


def clean(value):
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return ""
    return re.sub(r"\.0$", "", str(value).strip())


def load_crf(path, gene):
    # 0904 통합본에는 처음 6개 유전자 시트만 있다 (새 7개 유전자는 CRF 값 없이 기록)
    if gene not in pd.ExcelFile(path).sheet_names:
        return {}
    crf = pd.read_excel(path, sheet_name=gene, header=1, dtype=str)
    crf = crf[crf["Patient ID"].notna()]
    return {pid: rows for pid, rows in crf.groupby("Patient ID")}


def crf_value(crf_by_id, cp, variant_c, variant_p, crf_field):
    """같은 CP 안에서 변이가 맞는 CRF 행의 값을 찾는다 (CP가 하나뿐이면 그 행)."""
    rows = crf_by_id.get(cp)
    if rows is None or crf_field not in rows.columns:
        return ""
    if len(rows) > 1:
        key = lambda v: re.sub(r"\s+", "", clean(v)).lower()
        hit = rows[(rows["Variant cDNA change"].map(key) == key(variant_c)) |
                   (rows["Variant protein change"].map(key) == key(variant_p))]
        rows = hit if len(hit) else rows
    return clean(rows.iloc[0][crf_field])


def load_evidence_logs(paths):
    by_key = {}
    for path in paths:
        for row in csv.DictReader(open(path, encoding="utf-8-sig")):
            by_key.setdefault((row["mention_id"], row["action"], row["field"]), row)
    return by_key


def merge_rule_for(mention_id, gene):
    for row in csv.DictReader(open(MERGE_RULES, encoding="utf-8-sig")):
        if row["gene"] == gene and mention_id in row["mention_ids"].split(";"):
            return row
    return None


def variant_exclusion_for(mention_id, gene):
    for row in csv.DictReader(open(VARIANT_EXCLUSIONS, encoding="utf-8-sig")):
        if row["gene"] == gene and row.get("mention_id") == mention_id and row.get("status", "active") == "active":
            return row
    return None


def main():
    parser = argparse.ArgumentParser(description="CRF 통합본 기준 행 단위 변경 기록 생성")
    parser.add_argument("--gene", required=True)
    parser.add_argument("--before", required=True, type=Path, help="수정 전 canonical CSV")
    parser.add_argument("--after", required=True, type=Path, help="수정 후 canonical CSV")
    parser.add_argument("--db", type=Path, default=None, help="유전자 SQLite (기본: databases/{GENE}.sqlite)")
    parser.add_argument("--evidence-log", nargs="*", type=Path, default=[], help="적용 로그 CSV들")
    parser.add_argument("--lookup-canonical", nargs="*", type=Path, default=[],
                        help="--before에 없는 mention의 CP/환자 정보를 찾을 이전 canonical CSV들")
    parser.add_argument("--batch-id", required=True)
    parser.add_argument("--change-date", required=True)
    parser.add_argument("--reason", required=True, help="이 batch 변경의 요청 출처/사유 (한 줄)")
    parser.add_argument("--run-dir", default="")
    parser.add_argument("--crf", type=Path, default=DEFAULT_CRF)
    parser.add_argument("--output", type=Path, default=DEFAULT_LOG)
    args = parser.parse_args()

    gene = args.gene.upper()
    before = pd.read_csv(args.before, dtype=str).fillna("")
    after = pd.read_csv(args.after, dtype=str).fillna("")
    after_by_mid = after.set_index("mention_id")
    crf_by_id = load_crf(args.crf, gene)
    logs = load_evidence_logs(args.evidence_log)
    conn = sqlite3.connect(args.db or WORKSPACE_ROOT / "databases" / f"{gene}.sqlite")
    conn.row_factory = sqlite3.Row

    def base(row, change_type, field, before_value, after_value, reason, evidence, current_cp):
        cp = row["canonical_patient_id"]
        return {
            "change_date": args.change_date, "batch_id": args.batch_id, "gene": gene,
            "crf_patient_id": cp, "current_patient_id": current_cp,
            "pmid": clean(row["pmid"]), "patient_in_paper": row["patient_id"],
            "variant_c": row["variant_c"], "variant_p": row["variant_p"],
            "change_type": change_type, "field": field,
            "crf_value": crf_value(crf_by_id, cp, row["variant_c"], row["variant_p"], field) if field in FIELDS.values() else "",
            "before_value": before_value, "after_value": after_value,
            "reason": reason, "evidence": evidence, "mention_id": row["mention_id"], "run_dir": args.run_dir,
        }

    out = []
    for _, row in before.iterrows():
        mid = row["mention_id"]
        if mid in after_by_mid.index:
            new = after_by_mid.loc[mid]
            if row["canonical_patient_id"] != new["canonical_patient_id"]:
                rule = merge_rule_for(mid, gene)
                out.append(base(row, "merged_duplicate", "Patient ID", row["canonical_patient_id"], new["canonical_patient_id"],
                                args.reason, f"{rule['manual_cluster_id']}: {rule['evidence']}" if rule else "",
                                new["canonical_patient_id"]))
            for col, crf_field in FIELDS.items():
                if clean(row[col]) != clean(new[col]):
                    log = next((logs[k] for k in logs if k[0] == mid and LOG_FIELD_TO_CANONICAL.get(k[2]) == col), None)
                    if log:
                        evidence = log["evidence"]
                    elif col in ("clinvar_pathogenicity", "ClinVar ID"):
                        evidence = "변이 정정 후 ClinVar 재조회 결과"
                    else:
                        evidence = ""
                    out.append(base(row, "value_corrected", crf_field, row[col], new[col], args.reason, evidence,
                                    new["canonical_patient_id"]))
            continue
        collapsed = after[(after["canonical_patient_id"] == row["canonical_patient_id"]) &
                          after["supporting_mention_ids"].str.split("; ").map(lambda ids: mid in ids)]
        if len(collapsed):
            kept = collapsed.iloc[0]
            out.append(base(row, "same_variant_row_collapsed", "(row)", f"{row['variant_c']} / {row['variant_p']}",
                            f"{kept['variant_c']} / {kept['variant_p']}", args.reason,
                            "같은 환자의 같은 변이가 논문마다 다른 표기(cDNA 누락, frameshift 표기 차이 등)로 나뉘어 있던 행을 한 행으로 정리; "
                            f"PMID {row['pmid']}는 supporting_pmids에 유지", kept["canonical_patient_id"]))
            # 대표 행이 다른 mention으로 바뀌어 화면에 보이는 값이 달라진 경우도 남긴다.
            for col, crf_field in FIELDS.items():
                if clean(row[col]) != clean(kept[col]):
                    out.append(base(row, "value_corrected", crf_field, row[col], kept[col], args.reason,
                                    f"병합 후 대표 행이 PMID {clean(kept['pmid'])} ({kept['patient_id']})로 바뀜",
                                    kept["canonical_patient_id"]))
            continue
        rule = merge_rule_for(mid, gene)
        if rule:
            others = [m for m in rule["mention_ids"].split(";") if m != mid and m in after_by_mid.index]
            kept_cp = after_by_mid.loc[others[0]]["canonical_patient_id"] if others else ""
            out.append(base(row, "merged_duplicate", "Patient ID", row["canonical_patient_id"], kept_cp,
                            args.reason, f"{rule['manual_cluster_id']}: {rule['evidence']}", kept_cp))
            continue
        vex = variant_exclusion_for(mid, gene)
        if vex:
            out.append(base(row, "row_excluded", "(row)", "included", "excluded", args.reason,
                            f"{vex['reason']}: {vex['evidence']}", row["canonical_patient_id"]))
            continue
        db_row = conn.execute("SELECT mention_source_type, originality_evidence FROM patient_mentions WHERE mention_id=?",
                              (int(mid),)).fetchone()
        # 근거는 이번 작업 로그 → 정정 규칙(manual_exclude) → DB 추출 문구 순으로 쓴다.
        log = next((logs[k] for k in logs if k[0] == mid), None)
        rule_row = next((r for r in csv.DictReader(open(CORRECTIONS, encoding="utf-8-sig"))
                         if r["gene"] == gene and r["mention_id"] == mid and r["value"] == "manual_exclude"), None)
        if log:
            evidence = log["evidence"]
        elif rule_row:
            evidence = rule_row["evidence"]
        else:
            evidence = db_row["originality_evidence"].split(" | previous:")[0] if db_row else ""
        out.append(base(row, "row_excluded", "(row)", "included", "excluded", args.reason, evidence, ""))

    # 수정 후 새로 들어온 행
    before_mids = set(before["mention_id"])
    for _, new in after.iterrows():
        if new["mention_id"] in before_mids:
            continue
        log = next((logs[k] for k in logs if k[0] == new["mention_id"]), None)
        out.append(base(new, "row_added", "(row)", "", "included", args.reason,
                        log["evidence"] if log else "", new["canonical_patient_id"]))

    # canonical에 남지 않은 mention의 DB 값 정정 (병합되어 다른 행 값이 표시되는 경우)
    lookup = pd.concat([before] + [pd.read_csv(p, dtype=str).fillna("") for p in args.lookup_canonical])
    before_by_mid = lookup.drop_duplicates("mention_id").set_index("mention_id")
    for (mid, action, field), log in logs.items():
        col = LOG_FIELD_TO_CANONICAL.get(field)
        if action not in ("value_fix", "variant_fix") or col is None or mid in after_by_mid.index or mid not in before_by_mid.index:
            continue
        row = before_by_mid.loc[mid].copy()
        row["mention_id"] = mid
        rule = merge_rule_for(mid, gene)
        kept = ""
        if rule:
            others = [m for m in rule["mention_ids"].split(";") if m != mid and m in after_by_mid.index]
            kept = after_by_mid.loc[others[0]]["canonical_patient_id"] if others else ""
        out.append(base(row, "db_only_value_corrected", FIELDS[col], log["old"], log["new"], args.reason,
                        log["evidence"] + (f" (병합되어 canonical에는 {kept} 행 값이 표시됨)" if kept else ""), kept))

    existing = []
    if args.output.exists():
        existing = [r for r in csv.DictReader(open(args.output, encoding="utf-8-sig")) if not (r["batch_id"] == args.batch_id and r["gene"] == gene)]
    order = {"row_excluded": 0, "row_added": 1, "merged_duplicate": 2, "same_variant_row_collapsed": 3,
             "value_corrected": 4, "db_only_value_corrected": 5}
    out.sort(key=lambda r: (order[r["change_type"]], r["crf_patient_id"], r["field"]))
    # 기존 batch의 change_id는 유지하고, 새 행은 그 뒤 번호를 받는다.
    next_id = max((int(r["change_id"][3:]) for r in existing), default=0) + 1
    for i, r in enumerate(out, next_id):
        r["change_id"] = f"CHG{i:05d}"
    rows = existing + out
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with open(args.output, "w", newline="", encoding="utf-8-sig") as f:
        writer = csv.DictWriter(f, fieldnames=COLUMNS)
        writer.writeheader()
        writer.writerows(rows)
    counts = pd.Series([r["change_type"] for r in out]).value_counts().to_dict() if out else {}
    print(f"{args.batch_id}: {len(out)} rows {counts} -> {args.output} (total {len(rows)})")


if __name__ == "__main__":
    main()
