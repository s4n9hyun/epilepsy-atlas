"""Append a batch of correction rules (config/manual_value_corrections.csv format) and write per-gene evidence logs.

Usage: python apply_rule_batch.py <batch_dir> <superseded_tag> <rule.csv> [<rule.csv> ...] [--reason "..."]
Later files win on the same (gene, mention_id, field). Existing active rules with the same key are superseded with the tag.
Rules whose mention_id is not in the gene database are dropped and reported. Evidence logs: <batch_dir>/evidence_<GENE>.csv
(input for scripts/curation/rebuild_canonical.py --evidence-log). Prints the genes to rebuild.
"""
import csv
import sqlite3
import sys
from pathlib import Path

import pandas as pd

W = Path(__import__("os").environ.get("EPILEPSY_WORKSPACE", Path(__file__).resolve().parents[2]))
COLS = ["gene", "mention_id", "pmid", "patient_id", "field", "value", "previous_value", "evidence", "source_run", "status"]


def main():
    args = [a for a in sys.argv[1:]]
    reason = ""
    if "--reason" in args:
        i = args.index("--reason")
        reason = args[i + 1]
        del args[i:i + 2]
    batch, tag, files = Path(args[0]), args[1], args[2:]
    frames = []
    for k, f in enumerate(files):
        d = pd.read_csv(f, dtype=str, keep_default_na=False)
        if len(d):
            d["_order"] = k
            frames.append(d)
    r = pd.concat(frames, ignore_index=True)
    for c in COLS:
        if c not in r.columns:
            r[c] = ""
    r = r[r.mention_id.astype(str).str.strip() != ""]
    r["status"] = "active"
    r = r.sort_values("_order", kind="stable").drop_duplicates(["gene", "mention_id", "field"], keep="last")
    keep = []
    for g, d in r.groupby("gene"):
        ids = {str(x) for (x,) in sqlite3.connect(W / f"databases/{g}.sqlite").execute("select mention_id from patient_mentions")}
        bad = d[~d.mention_id.isin(ids)]
        if len(bad):
            print(f"dropped {len(bad)} {g} rules with unknown mention_id")
        keep.append(d[d.mention_id.isin(ids)])
    r = pd.concat(keep)
    c = pd.read_csv(W / "config/manual_value_corrections.csv", dtype=str, keep_default_na=False)
    key = set(zip(r.gene, r.mention_id, r.field))
    mm = [(g, m, f) in key and s == "active" for g, m, f, s in zip(c.gene, c.mention_id, c.field, c.status)]
    c.loc[mm, "status"] = tag
    pd.concat([c, r[COLS]], ignore_index=True).to_csv(W / "config/manual_value_corrections.csv", index=False, quoting=csv.QUOTE_MINIMAL)
    batch.mkdir(parents=True, exist_ok=True)
    for g, d in r.groupby("gene"):
        pd.DataFrame(dict(action=["exclude" if v == "manual_exclude" else "correct" for v in d.value], mention_id=d.mention_id,
                          field=d.field, old=d.previous_value, new=d.value,
                          evidence=(reason + " | " if reason else "") + d.evidence)).to_csv(batch / f"evidence_{g}.csv", index=False)
    print("rules", len(r), "superseded", sum(mm), r.groupby(["gene"]).size().to_dict())
    print("GENES", " ".join(sorted(r.gene.unique())))


if __name__ == "__main__":
    main()
