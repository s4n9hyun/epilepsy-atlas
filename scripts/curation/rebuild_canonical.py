#!/usr/bin/env python3
"""유전자 canonical을 한 번에 다시 만들고, 확인 후 운영 산출물로 교체한다.

단계:
  1. config/manual_value_corrections.csv 정정을 DB에 적용 (재추출 후에도 정정 유지)
  2. build_review_tables.py로 runs/rebuild/{시각}/{GENE}/에 빌드 (최근 3회분만 유지) (lock은 config/canonical_id_locks/{GENE}.csv)
  3. 현재 outputs/current와 대조: 행/환자 수, 없어진·새 CP, 남은 행의 CP 변경, 바뀐 컬럼
  4. --promote일 때만:
     - 현재 canonical/auxiliary/lock을 backups/output_history/{GENE}/{시각}/로 보관
     - 새 canonical → outputs/current/canonical_patients/{GENE}_canonical_patients/
     - 보조 산출물 → outputs/auxiliary/{GENE}/
     - lock 갱신: 기존 lock + 이번 빌드의 mention→CP (새 환자 번호 고정)
     - --batch-id가 있으면 outputs/current/curation_change_log.csv에 변경 기록 추가

사용 예:
  python3 scripts/curation/rebuild_canonical.py --gene KCNQ2                  # 대조만
  python3 scripts/curation/rebuild_canonical.py --gene KCNQ2 --promote \
      --batch-id 20261001_KCNQ2_xxx --reason "사유" --evidence-log runs/.../applied_changes.csv
"""
import argparse
import csv
import shutil
import subprocess
import sys
from datetime import datetime
from pathlib import Path

import pandas as pd

WORKSPACE_ROOT = Path(__file__).resolve().parents[2]
CORE = WORKSPACE_ROOT / "scripts" / "core"
CURATION = Path(__file__).resolve().parent
CURRENT = WORKSPACE_ROOT / "outputs" / "current" / "canonical_patients"
AUXILIARY = WORKSPACE_ROOT / "outputs" / "auxiliary"
LOCKS = WORKSPACE_ROOT / "config" / "canonical_id_locks"
HISTORY = WORKSPACE_ROOT / "backups" / "output_history"
AUX_KINDS = ["patient_mentions", "excluded_mentions", "patient_dedup_links"]


def current_csv(gene):
    files = sorted((CURRENT / f"{gene}_canonical_patients").glob(f"{gene}_canonical_patients_*.csv"))
    return files[-1] if files else None


def build(gene, out_root, stamp):
    subprocess.run([sys.executable, "build_review_tables.py", str(WORKSPACE_ROOT / "databases" / f"{gene}.sqlite"),
                    "--output-dir", str(out_root), "--gene", gene, "--experimental-strict-source-filter",
                    "--timestamp", stamp], cwd=CORE, check=True, capture_output=True, text=True)
    return out_root / f"{gene}_canonical_patients"


def compare(old_path, new_path):
    old = pd.read_csv(old_path, dtype=str).fillna("") if old_path else pd.DataFrame(columns=["mention_id", "canonical_patient_id"])
    new = pd.read_csv(new_path, dtype=str).fillna("")
    both = old.merge(new, on="mention_id", suffixes=("_o", "_n"))
    cols = [c for c in old.columns if c in new.columns and c not in ("mention_id", "extraction_timestamp")]
    changed = {c: int((both[c + "_o"] != both[c + "_n"]).sum()) for c in cols}
    return {
        "rows": (len(old), len(new)),
        "patients": (old["canonical_patient_id"].nunique(), new["canonical_patient_id"].nunique()),
        "retired_cp": sorted(set(old["canonical_patient_id"]) - set(new["canonical_patient_id"])),
        "new_cp": sorted(set(new["canonical_patient_id"]) - set(old["canonical_patient_id"])),
        "surviving_cp_changed": int((both["canonical_patient_id_o"] != both["canonical_patient_id_n"]).sum()),
        "changed_columns": {c: n for c, n in changed.items() if n},
        "removed_rows": len(set(old["mention_id"]) - set(new["mention_id"])),
        "added_rows": len(set(new["mention_id"]) - set(old["mention_id"])),
    }


def update_locks(gene, links_path, stamp):
    lock_path = LOCKS / f"{gene}.csv"
    rows = list(csv.DictReader(lock_path.open(encoding="utf-8-sig"))) if lock_path.exists() else []
    fields = list(rows[0].keys()) if rows else ["gene", "mention_id", "canonical_patient_id", "pmid", "pmcid",
                                                "patient_id", "lock_source", "locked_at"]
    by_mid = {r["mention_id"]: r for r in rows}
    links = pd.read_csv(links_path, dtype=str).fillna("")
    added = updated = 0
    for _, link in links.iterrows():
        mid, cp = link["mention_id"], link["canonical_patient_id"]
        if not mid or not cp:
            continue
        if mid not in by_mid:
            by_mid[mid] = {**{k: "" for k in fields}, "gene": gene, "mention_id": mid, "canonical_patient_id": cp,
                           "pmid": link.get("pmid", ""), "pmcid": link.get("pmcid", ""), "patient_id": link.get("patient_id", ""),
                           "lock_source": f"rebuild_canonical_{stamp}", "locked_at": datetime.now().isoformat(timespec="seconds")}
            added += 1
        elif by_mid[mid]["canonical_patient_id"] != cp:
            by_mid[mid]["canonical_patient_id"] = cp
            by_mid[mid]["lock_source"] = f"rebuild_canonical_{stamp}"
            updated += 1
    with lock_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(by_mid.values())
    return added, updated


def main():
    parser = argparse.ArgumentParser(description="canonical 재생성·대조·교체")
    parser.add_argument("--gene", nargs="+", required=True)
    parser.add_argument("--promote", action="store_true", help="대조 후 운영 산출물로 교체")
    parser.add_argument("--skip-corrections", action="store_true", help="정정 규칙 적용 생략")
    parser.add_argument("--batch-id", help="변경 기록 batch ID (promote 시)")
    parser.add_argument("--reason", default="", help="변경 기록 사유")
    parser.add_argument("--evidence-log", nargs="*", default=[], help="변경 근거 로그 CSV")
    args = parser.parse_args()

    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    out_root = WORKSPACE_ROOT / "runs" / "rebuild" / stamp
    # 빌드 중간 결과는 최근 3회분만 남긴다 (반영 전 canonical은 backups/output_history/에 보관)
    for old in sorted(out_root.parent.glob("*_*"))[:-2]:
        shutil.rmtree(old)
    for gene in [g.upper() for g in args.gene]:
        if not args.skip_corrections:
            subprocess.run([sys.executable, str(CURATION / "apply_manual_corrections.py"), "--gene", gene], check=True)
        new_dir = build(gene, out_root, stamp)
        new_csv = new_dir / f"{gene}_canonical_patients_{stamp}.csv"
        old_csv = current_csv(gene)
        report = compare(old_csv, new_csv)
        print(f"== {gene}: rows {report['rows'][0]} -> {report['rows'][1]}, patients {report['patients'][0]} -> {report['patients'][1]}, "
              f"removed {report['removed_rows']}, added {report['added_rows']}, surviving CP changed {report['surviving_cp_changed']}")
        print(f"   retired CP: {report['retired_cp'][:20]}{' ...' if len(report['retired_cp']) > 20 else ''}")
        print(f"   new CP: {report['new_cp'][:20]}{' ...' if len(report['new_cp']) > 20 else ''}")
        print(f"   changed columns: {report['changed_columns']}")
        print(f"   build: {new_dir}")
        if not args.promote:
            continue

        history = HISTORY / gene / stamp
        history.mkdir(parents=True, exist_ok=True)
        if old_csv:
            shutil.move(str(old_csv), history / old_csv.name)
        if (AUXILIARY / gene).exists():
            shutil.move(str(AUXILIARY / gene), history / "auxiliary")
        if (LOCKS / f"{gene}.csv").exists():
            shutil.copy2(LOCKS / f"{gene}.csv", history / f"{gene}_canonical_id_locks.csv")
        (CURRENT / f"{gene}_canonical_patients").mkdir(parents=True, exist_ok=True)
        shutil.copy2(new_csv, CURRENT / f"{gene}_canonical_patients" / new_csv.name)
        (AUXILIARY / gene).mkdir(parents=True, exist_ok=True)
        for kind in AUX_KINDS:
            shutil.copy2(new_dir / f"{gene}_{kind}_{stamp}.csv", AUXILIARY / gene)
        shutil.copy2(new_dir / "dedup_summary.json", AUXILIARY / gene)
        added, updated = update_locks(gene, new_dir / f"{gene}_patient_dedup_links_{stamp}.csv", stamp)
        print(f"   promoted; previous kept in {history.relative_to(WORKSPACE_ROOT)}; locks added {added}, updated {updated}")
        if args.batch_id:
            subprocess.run([sys.executable, str(CURATION / "build_curation_change_log.py"), "--gene", gene,
                            "--before", str(history / old_csv.name), "--after", str(CURRENT / f"{gene}_canonical_patients" / new_csv.name),
                            "--batch-id", args.batch_id, "--change-date", datetime.now().strftime("%Y-%m-%d"),
                            "--reason", args.reason, "--run-dir", str(out_root.relative_to(WORKSPACE_ROOT)),
                            *(["--evidence-log", *args.evidence_log] if args.evidence_log else [])], check=True)


if __name__ == "__main__":
    main()
