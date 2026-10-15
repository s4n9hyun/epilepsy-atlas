"""Curation quality records reported in Methods and the Supporting Information (reviewer comments R2-M5, R4-M3).

- Reconciliation of clinician case report form (CRF) values with the source: tables/supp_crf_reconciliation.csv
- Transcript/variant quality control against the reference sequence: tables/supp_variant_qc.csv
Merges results into analysis/results.json ("crf_check", "variant_qc").
"""
import json

from common import PROJECT, TAB, WORKSPACE, pd

CRF = WORKSPACE / "runs/crf_0904/adopted_value_source_check/checked_all.csv"
VQ = WORKSPACE / "runs/review_fixes/variant_qc"
VERDICT = {"match": "Confirmed at the source (CRF value adopted)", "partial": "Partly confirmed (source value kept or merged)",
           "conflict": "Contradicted by the source (source value kept)", "not_reported": "Not reported in the source (not adopted)",
           "cohort_level_only": "Cohort-level statement only (not adopted)", "patient_not_found": "Patient not found in the article",
           "external_not_checked": "External source, not checked"}


def crf():
    c = pd.read_csv(CRF, dtype=str, keep_default_na=False)
    t = c.verdict.map(VERDICT).value_counts().rename_axis("Result of source check").reset_index(name="CRF values, n")
    t["%"] = (100 * t["CRF values, n"] / len(c)).round(1)
    by = pd.crosstab(c.field, c.verdict == "match").rename(columns={True: "confirmed", False: "not confirmed"})
    t.to_csv(TAB / "supp_crf_reconciliation.csv", index=False)
    by.reset_index().to_csv(TAB / "supp_crf_reconciliation_by_field.csv", index=False)
    nc = int((c.verdict != "match").sum())
    return {"values": len(c), "confirmed": int((c.verdict == "match").sum()), "not_confirmed": nc,
            "not_confirmed_pct": round(100 * nc / len(c), 1), "by_verdict": c.verdict.value_counts().to_dict()}


def variant_qc():
    r = pd.read_csv(VQ / "report.csv", dtype=str, keep_default_na=False)
    flag = pd.read_csv(VQ / "source_check_needed.csv", dtype=str, keep_default_na=False)
    agent = VQ / "agent/decisions.csv"
    ad = pd.read_csv(agent, dtype=str, keep_default_na=False) if agent.exists() else pd.DataFrame(columns=["decision"])
    cat = []
    for d, w in zip(r.decision, r.why):
        if not w:
            cat.append("Already on MANE numbering and consistent with the reference")
        elif w.startswith("fits MANE") or w.startswith("anchor positions"):
            cat.append("Transcript relabelled to MANE (identical coding sequence at the variant)")
        elif w.startswith("legacy"):
            cat.append("Legacy numbering converted to MANE")
        elif w.startswith("transcript label aligned") or w.startswith("c. prefix"):
            cat.append("Notation completed (c. prefix or transcript label)")
        elif w.startswith("SCN8A exon 5N") or w.startswith("fits only NM_003159"):
            cat.append("Kept on alternative transcript (exon absent from MANE)")
        elif w.startswith("SCN8A exon 5"):
            cat.append("SCN8A exon 5 isoform: kept on MANE (isoform not stated or exon 5A)")
        else:
            cat.append("Other")
    t = pd.Series(cat).value_counts().rename_axis("Check").reset_index(name="Patient-variant rows, n")
    extra = pd.DataFrame([{"Check": "cDNA and stated protein inconsistent: re-checked at the source", "Patient-variant rows, n": len(flag)}])
    t = pd.concat([t, extra], ignore_index=True)
    LAB = {"article_inconsistent": "article itself inconsistent; printed notation kept, not placed on MANE",
           "article_inconsistent_resolved": "resolved by independent evidence (genomic coordinate, ClinVar, repeated mention)",
           "cdna_ref_mismatch": "printed cDNA does not match the reference; not placed on MANE",
           "legacy_numbering": "legacy numbering; converted", "protein_ref_label_error": "protein label error in the article; cDNA confirmed",
           "checked; already matches its column": "confirmed as recorded"}
    if len(ad):
        ad = ad[ad.part.astype(str).str.upper().str.startswith("A")] if "part" in ad.columns else ad
        cats = ad.decision.str.split(":").str[0].map(lambda k: LAB.get(k, "confirmed as recorded"))
        for k, v in cats.value_counts().items():
            t.loc[len(t)] = {"Check": f"  of which: {k}", "Patient-variant rows, n": int(v)}
    t.to_csv(TAB / "supp_variant_qc.csv", index=False)
    n = len(r)
    return {"rows_checked": n, "relabelled": int(sum(x.startswith("Transcript relabelled") for x in cat)),
            "converted": int(sum(x.startswith("Legacy") for x in cat)), "notation": int(sum(x.startswith("Notation") for x in cat)),
            "alt_kept": int(sum(x.startswith("Kept on alternative") for x in cat)), "inconsistent": len(flag),
            "inconsistent_pct": round(100 * len(flag) / n, 1), "agent_decisions": ad.decision.str.split(":").str[0].value_counts().to_dict() if len(ad) else {}}


def review_fixes():
    """Duplicate adjudication, inclusion re-check, ClinVar B/LB exclusion, abstract-only articles, families, software."""
    import sys
    import matplotlib
    import numpy
    import scipy
    import statsmodels
    from common import family_ids, load_rows, load_patients
    RF = WORKSPACE / "runs/review_fixes"
    dup = pd.concat([pd.read_csv(f, dtype=str, keep_default_na=False) for f in sorted((RF / "dup/results").glob("U*.csv"))])
    inc = pd.concat([pd.read_csv(f, dtype=str, keep_default_na=False) for f in sorted((RF / "inclusion/results").glob("I*.csv"))])
    exc = pd.read_csv(RF / "apply/rules_inclusion.csv", dtype=str, keep_default_na=False)
    blb_f = RF / "apply/rules_blb.csv"
    blb = pd.read_csv(blb_f, dtype=str, keep_default_na=False) if blb_f.exists() else pd.DataFrame()
    r = load_rows()
    loc = (r.location_in_paper.astype(str) + " " + r.originality_evidence.astype(str)).str.lower()
    ab = r[loc.str.contains(r"only (?:the )?abstract|abstract only|abstract[- ]only|abstract accessible")]
    fid, _ = family_ids(load_patients(r))
    return {"dup_pairs": len(dup), "dup_same": int((dup.verdict == "same").sum()), "dup_different": int((dup.verdict == "different").sum()),
            "dup_uncertain": int((dup.verdict == "uncertain").sum()), "inc_checked": len(inc), "inc_excluded": len(exc),
            "blb_excluded": int(blb.drop_duplicates(["gene", "pmid", "patient_id"]).shape[0]) if len(blb) else 0, "abstract_only_articles": int(ab.pmid.nunique()),
            "abstract_only_patients": int(ab.canonical_patient_id.astype(str).add(ab.gene).nunique()),
            "families": int(pd.Series(fid.values).nunique()),
            "software": {"python": sys.version.split()[0], "pandas": pd.__version__, "numpy": numpy.__version__, "scipy": scipy.__version__,
                         "statsmodels": statsmodels.__version__, "matplotlib": matplotlib.__version__}}


def main():
    R = {"crf_check": crf(), "variant_qc": variant_qc(), "review_fixes": review_fixes()}
    res = json.loads((PROJECT / "analysis" / "results.json").read_text())
    res.update(R)
    (PROJECT / "analysis" / "results.json").write_text(json.dumps(res, indent=1, default=str), encoding="utf-8")
    print(json.dumps(R, indent=1))


if __name__ == "__main__":
    main()
