#!/usr/bin/env python3
"""config/manual_value_corrections.csv의 원문 확인 정정을 유전자 DB에 적용한다.

재추출로 patient_mentions 값이 덮여도 이 스크립트를 다시 돌리면 정정이 복원된다.
- 이미 같은 값이면 건너뛴다 (여러 번 실행해도 안전).
- mention_id가 없거나 PMID가 다르면 적용하지 않고 'unresolved'로 보고한다
  (재추출로 mention_id가 바뀐 경우 → 규칙의 mention_id를 새 번호로 고쳐야 함).

사용:
  python3 scripts/curation/apply_manual_corrections.py --gene CDKL5 [--dry-run]
"""
import argparse
import csv
import json
import sqlite3
from datetime import datetime
from pathlib import Path

WORKSPACE_ROOT = Path(__file__).resolve().parents[2]
CORRECTIONS = WORKSPACE_ROOT / "config" / "manual_value_corrections.csv"
FIELDS = ["gene", "mention_id", "pmid", "patient_id", "field", "value", "previous_value", "evidence", "source_run", "status"]


def load_corrections(gene):
    with CORRECTIONS.open(newline="", encoding="utf-8-sig") as f:
        return [r for r in csv.DictReader(f) if r["gene"].upper() == gene and r.get("status", "active") == "active"]


def apply(gene, dry_run=False):
    rows = load_corrections(gene)
    conn = sqlite3.connect(WORKSPACE_ROOT / "databases" / f"{gene}.sqlite")
    conn.row_factory = sqlite3.Row
    columns = {r[1] for r in conn.execute("PRAGMA table_info(patient_mentions)")}
    applied, same, unresolved = [], 0, []
    for rule in rows:
        row = conn.execute("SELECT * FROM patient_mentions WHERE mention_id=?", (int(rule["mention_id"]),)).fetchone()
        if row is None or str(row["pmid"]) != rule["pmid"] or rule["field"] not in columns:
            unresolved.append(rule)
            continue
        if (row[rule["field"]] or "") == rule["value"]:
            same += 1
            continue
        applied.append((rule, row[rule["field"]]))
        if not dry_run:
            raw = json.loads(row["raw_patient_json"] or "{}")
            raw.setdefault("manual_corrections", []).append({
                "reapplied_at": datetime.now().isoformat(timespec="seconds"),
                "source": "config/manual_value_corrections.csv",
                "evidence": rule["evidence"],
                "previous": {rule["field"]: row[rule["field"]]},
            })
            conn.execute(f"UPDATE patient_mentions SET {rule['field']}=?, raw_patient_json=? WHERE mention_id=?",
                         (rule["value"], json.dumps(raw, ensure_ascii=False), int(rule["mention_id"])))
    if not dry_run:
        conn.commit()
    conn.close()
    return rows, applied, same, unresolved


def main():
    parser = argparse.ArgumentParser(description="원문 확인 정정 규칙을 DB에 적용")
    parser.add_argument("--gene", nargs="+", required=True)
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()
    for gene in [g.upper() for g in args.gene]:
        rows, applied, same, unresolved = apply(gene, args.dry_run)
        print(f"{gene}: rules {len(rows)}, already applied {same}, "
              f"{'would apply' if args.dry_run else 'applied'} {len(applied)}, unresolved {len(unresolved)}")
        for rule, old in applied[:20]:
            print(f"  {rule['mention_id']} {rule['field']}: {old!r} -> {rule['value']!r}")
        for rule in unresolved[:20]:
            print(f"  UNRESOLVED mention {rule['mention_id']} (PMID {rule['pmid']}, {rule['patient_id']}) {rule['field']}")


if __name__ == "__main__":
    main()
