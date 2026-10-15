"""Fallback for atlas variants that NCBI Variation Services could not place on GRCh38 (mostly intronic CDKL5 variants).

Unmapped rows of runs/review_fixes/clinvar/allele_map.csv are matched by transcript HGVS against ClinVar hgvs4variation
(same release), with deleted/duplicated bases dropped on both sides. Matched rows get the ClinVar classification and
status 'matched by HGVS'; the remaining rows are 'Not in ClinVar' with no genomic allele (counted by their HGVS).
Rewrites allele_map.csv in place, adding column match_method.
"""
import re
from pathlib import Path

import pandas as pd

W = Path(__import__("os").environ.get("EPILEPSY_WORKSPACE", Path(__file__).resolve().parents[2]))
OUT = Path(__import__("os").environ.get("CLINVAR_MAP_DIR", W / "runs/review_fixes/clinvar"))
SNAP = W / "runs/paper/clinvar_snapshot"
sys_path = str(W / "scripts/curation")


def key(h):
    h = re.sub(r"\s", "", str(h))
    h = re.sub(r"(del|dup)[ACGT]+(?=ins|$)", r"\1", h)
    return h


def main():
    import sys
    sys.path.insert(0, sys_path)
    from clinvar_allele_match import status
    am = pd.read_csv(OUT / "allele_map.csv", dtype=str, keep_default_na=False)
    if "match_method" not in am.columns:
        am["match_method"] = ["" if s == "Unmapped" else "allele (GRCh38)" for s in am.status]
    hv = pd.read_csv(SNAP / "hgvs4variation_13genes.tsv", sep="\t", dtype=str, keep_default_na=False,
                     usecols=["VariationID", "NucleotideExpression"])
    vs = pd.read_csv(SNAP / "variant_summary_13genes.tsv", sep="\t", dtype=str, keep_default_na=False)
    sig_col = "ClinicalSignificance" if "ClinicalSignificance" in vs.columns else "GermlineClassification"
    vs = vs[vs.Assembly == "GRCh38"].drop_duplicates("VariationID").set_index("VariationID")
    idx = {}
    for vid, ne in zip(hv.VariationID, hv.NucleotideExpression):
        if ne.startswith("NM_"):
            idx.setdefault(key(ne), set()).add(vid)
    n = 0
    # also re-check alleles placed on GRCh38 but not found in ClinVar (diagnostic for normalisation differences)
    for i in am.index[(am.status == "Not in ClinVar") & (am.match_method == "allele (GRCh38)")]:
        vids = sorted(v for v in idx.get(key(am.at[i, "hgvs"]), set()) if v in vs.index)
        if vids:
            sigs = sorted({vs.at[v, sig_col] for v in vids})
            am.loc[i, ["vcv", "significance", "review_status", "status", "match_method"]] = [
                ";".join(vids), " | ".join(sigs), " | ".join(sorted({vs.at[v, "ReviewStatus"] for v in vids})),
                status(" | ".join(sigs)), "transcript HGVS (allele not matched)"]
    for i in am.index[am.status == "Unmapped"]:
        vids = sorted(v for v in idx.get(key(am.at[i, "hgvs"]), set()) if v in vs.index)
        if vids:
            sigs = sorted({vs.at[v, sig_col] for v in vids})
            am.loc[i, ["vcv", "significance", "review_status", "status", "match_method"]] = [
                ";".join(vids), " | ".join(sigs), " | ".join(sorted({vs.at[v, "ReviewStatus"] for v in vids})),
                status(" | ".join(sigs)), "transcript HGVS"]
            n += 1
        else:
            am.loc[i, ["status", "match_method"]] = ["Not in ClinVar", "not placed on GRCh38"]
    am.to_csv(OUT / "allele_map.csv", index=False)
    print("allele not matched but HGVS matched", int((am.match_method == "transcript HGVS (allele not matched)").sum()))
    print("fallback matched", n, "; remaining without allele", int((am.match_method == "not placed on GRCh38").sum()))
    print(am.status.value_counts().to_string())


if __name__ == "__main__":
    main()
