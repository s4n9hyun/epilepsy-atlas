"""Deterministic inclusion rule: exclude patient-variant rows whose only reported variant in the gene is classified benign/likely
benign in ClinVar with a review status of at least two stars (multiple submitters without conflicts, expert panel, practice guideline).

Reads runs/review_fixes/clinvar/allele_map.csv (after clinvar_hgvs_fallback.py) and the canonical rows; writes
runs/review_fixes/apply/rules_blb.csv in config/manual_value_corrections.csv format (mention_source_type -> manual_exclude) for
every supporting mention of such a row. Apply with apply_rule_batch.py and rebuild the genes.
"""
import os
import re
import sys
from pathlib import Path

import pandas as pd

W = Path(__import__("os").environ.get("EPILEPSY_WORKSPACE", Path(__file__).resolve().parents[2]))
sys.path.insert(0, str(W / "paper/analysis"))
from common import load_rows  # noqa: E402

CV = Path(os.environ.get("CLINVAR_MAP_DIR", W / "runs/review_fixes/clinvar"))
OUT = W / "runs/review_fixes/apply/rules_blb.csv"


def two_star(review):
    r = str(review).lower()
    return "practice guideline" in r or "expert panel" in r or ("multiple submitters" in r and "no conflict" in r)


def main():
    am = pd.read_csv(CV / "allele_map.csv", dtype=str, keep_default_na=False)
    am = am[(am.status == "Benign/likely benign") & am.review_status.map(two_star)]
    blb = {(g, h): (v, s, r) for g, h, v, s, r in zip(am.gene, am.hgvs, am.vcv, am.significance, am.review_status)}
    rows = load_rows()
    rows["hgvs"] = [f"{t}:{re.sub(r'[ ]', '', str(c))}" for t, c in zip(rows.gnomad_transcript, rows.variant_c)]
    nvar = rows.groupby(["gene", "canonical_patient_id"]).size()
    out = []
    for r in rows.itertuples():
        key = (r.gene, r.hgvs)
        if key not in blb or nvar[(r.gene, r.canonical_patient_id)] > 1:  # another variant in the gene: keep for manual review
            continue
        vcv, sig, rev = blb[key]
        for mid in re.split(r"[;,| ]+", str(r.supporting_mention_ids)):
            if mid:
                out.append(dict(gene=r.gene, mention_id=mid, pmid=r.pmid, patient_id=r.patient_id, field="mention_source_type",
                                value="manual_exclude", previous_value="original_patient",
                                evidence=f"Inclusion rule: ClinVar {sig} ({rev}; VCV {vcv}; release 2026-09-29) for {r.hgvs}; only variant in the gene.",
                                source_run="runs/review_fixes/clinvar", status="active"))
    pd.DataFrame(out).to_csv(OUT, index=False)
    print("B/LB >=2-star rows", len({(o['gene'], o['pmid'], o['patient_id']) for o in out}), "mentions", len(out))


if __name__ == "__main__":
    main()
