"""Accuracy of LLM extraction across all 13 genes, measured against the final source-verified dataset.

Every article with extracted patients was re-read against the source, so the final dataset is the reference.
- Recall: final patient-variant rows whose record came from LLM extraction, versus rows added during source review or by
  deterministic table parsing (raw_patient_json marked manual_review / parsed_from_table).
- Field values: for LLM-extracted rows, the value the model returned (raw_patient_json, kept unchanged by later corrections)
  versus the final value. Changes that are normalisation rather than error are not counted: formatting, genomic or protein-only
  notation converted to cDNA, legacy-transcript renumbering with the known offsets, renamed phenotype labels, and protein changes
  translated from the cDNA when the model returned none. Every other change (value corrected, missing value filled, value removed)
  counts as an error, so the estimate is an upper bound (it also counts typographical errors in the source that were corrected).
- A subset of 1,359 articles was additionally adjudicated field by field (runs/seven_genes/source_check); reported as Table S7b.
Writes tables/supp_S7a_extraction_accuracy.csv, tables/supp_S7b_adjudicated_subset.csv and merges results into results.json.
"""
import glob
import json
import re
import sqlite3

from common import GENES, PROJECT, TAB, WORKSPACE, pd

TOP = ["variant_c_original", "variant_p", "sex"]
PHE = ["age_of_onset", "last_followup_age", "seizure_presence", "phenotype_review_category", "targeted_treatment_response"]
LABEL = {"variant_c_original": "cDNA variant", "variant_p": "Protein change", "sex": "Sex", "age_of_onset": "Age at onset",
         "last_followup_age": "Age at last follow-up", "seizure_presence": "Seizure presence",
         "phenotype_review_category": "Phenotype category", "targeted_treatment_response": "Treatment response"}
NS = {"", "not specified", "unknown", "not available", "not applicable", "nan", "none", "na", "n/a", "nr", "not reported"}
OFFSETS = {0, 33, 111, 327}      # legacy-transcript cDNA offsets seen in the data (e.g. SCN1A NM_006920 +33)
AA_OFFSETS = {0, 11, 37, 109}


def norm(v):
    s = re.sub(r"\s+", " ", str(v if v is not None else "")).strip().lower().rstrip(".")
    return "" if s in NS else s


def cparts(v):
    v = re.sub(r"\s", "", re.sub(r"^[A-Z]{2}_\d+(\.\d+)?[:\s]*", "", str(v).strip()))
    return re.search(r"(\d+)(?:[+-]\d+)?(?:_\d+(?:[+-]\d+)?)?([ACGTacgt]*)(>|del|dup|ins)([ACGTacgt]*)", v)


def change_class(field, llm, final):
    """'same', 'normalisation' or 'error' (incl. omission and removal)."""
    a, b = norm(llm), norm(final)
    if a == b:
        return "same"
    if field == "variant_c_original" and a:
        if re.search(r"chr|g\.|^[a-z]\d+[a-z]$|^p\.|^[a-z]{3}\d+", a):
            return "normalisation"
        ma, mb = cparts(llm), cparts(final)
        if ma and mb and (ma.group(2).upper(), ma.group(3), ma.group(4).upper()) == (mb.group(2).upper(), mb.group(3), mb.group(4).upper()) \
                and abs(int(mb.group(1)) - int(ma.group(1))) in OFFSETS:
            return "normalisation"
    if field == "variant_p":
        if not a:
            return "normalisation"  # translated from the cDNA
        A, B = re.search(r"([a-z]{3})(\d+)", a), re.search(r"([a-z]{3})(\d+)", b)
        if A and B and A.group(1) == B.group(1) and abs(int(B.group(2)) - int(A.group(2))) in AA_OFFSETS and a[A.end():] == b[B.end():]:
            return "normalisation"
    if field == "phenotype_review_category" and a and b and (a in b or b in a):
        return "normalisation"
    return "error"


def main():
    rows, recs = [], []
    for g in GENES:
        can = pd.read_csv(sorted(glob.glob(str(WORKSPACE / f"outputs/current/canonical_patients/{g}_canonical_patients/*.csv")))[-1],
                          dtype=str, keep_default_na=False)
        con = sqlite3.connect(WORKSPACE / "databases" / f"{g}.sqlite")
        con.row_factory = sqlite3.Row
        for mid in can.mention_id.unique():
            r = con.execute("select * from patient_mentions where mention_id=?", (int(mid),)).fetchone()
            raw = json.loads(r["raw_patient_json"] or "{}")
            src = "review" if ("manual_review" in raw or "parsed_from_table" in raw) else "llm"
            rows.append({"gene": g, "source": src})
            if src != "llm":
                continue
            ph = raw.get("phenotype") if isinstance(raw.get("phenotype"), dict) else {}
            for f in TOP + PHE:
                llm = raw.get(f) if f in TOP else ph.get(f, raw.get(f))
                recs.append({"gene": g, "field": f, "cls": change_class(f, llm, r[f]), "llm_empty": norm(llm) == "",
                             "final_empty": norm(r[f]) == ""})
    rows, recs = pd.DataFrame(rows), pd.DataFrame(recs)
    t = []
    for f in TOP + PHE:
        x = recs[recs.field == f]
        err = int((x.cls == "error").sum())
        t.append({"Field": LABEL[f], "LLM-extracted values": len(x), "Changed on source review": err,
                  "Unchanged, %": round(100 * (1 - err / len(x)), 1)})
    t = pd.DataFrame(t)
    # non-empty view (reviewer comment R2-M6): changed among values the model returned; filled where the model returned nothing
    ne = []
    for f in TOP + PHE:
        x = recs[recs.field == f]
        nn = x[~x.llm_empty]
        ee = x[x.llm_empty & x.final_empty]
        ne.append({"Field": LABEL[f], "Non-empty LLM values": len(nn), "Changed, % of non-empty": round(100 * (nn.cls == "error").mean(), 1) if len(nn) else None,
                   "Filled on review (LLM empty)": int((x.llm_empty & ~x.final_empty).sum()), "Empty in both": len(ee)})
    ne = pd.DataFrame(ne)
    ne.to_csv(TAB / "supp_S7c_extraction_nonempty.csv", index=False)
    tot_n, tot_e = int(t["LLM-extracted values"].sum()), int(t["Changed on source review"].sum())
    t = pd.concat([t, pd.DataFrame([{"Field": "All fields", "LLM-extracted values": tot_n, "Changed on source review": tot_e,
                                     "Unchanged, %": round(100 * (1 - tot_e / tot_n), 1)}])])
    t.to_csv(TAB / "supp_S7a_extraction_accuracy.csv", index=False)
    byg = rows.groupby("gene").source.apply(lambda s: round(100 * (s == "review").mean(), 1)).to_dict()

    # adjudicated subset (field-by-field verdicts)
    SC = WORKSPACE / "runs/seven_genes/source_check"
    a = pd.read_csv(SC / "all_A.csv", dtype=str, keep_default_na=False)
    pat = a[a.field == "(patient)"]
    found = pat[pat.verdict == "found"]
    AF = {"variant_c": "cDNA variant", "variant_p": "Protein change", "sex": "Sex", "age_of_onset": "Age at onset",
          "last_followup_age": "Age at last follow-up", "seizure_presence": "Seizure presence",
          "phenotype_review_category": "Phenotype category", "targeted_treatment_response": "Treatment response",
          "mention_source_type": "Originality"}
    f = a[a.field.isin(AF) & a.verdict.isin(["wrong", "partial", "not_in_source", "cohort_level"])].drop_duplicates(["gene", "mention_id", "field"])
    s = pd.DataFrame([{"Field": lab, "Checked": len(found), "Wrong": int(((f.field == k) & (f.verdict == "wrong")).sum()),
                       "Partly correct": int(((f.field == k) & (f.verdict == "partial")).sum()),
                       "Not in source / cohort-level": int(((f.field == k) & f.verdict.isin(["not_in_source", "cohort_level"])).sum())}
                      for k, lab in AF.items()])
    s["Correct, %"] = (100 * (1 - s[["Wrong", "Partly correct", "Not in source / cohort-level"]].sum(axis=1) / s.Checked)).round(1)
    s.to_csv(TAB / "supp_S7b_adjudicated_subset.csv", index=False)
    sub_err = int(s[["Wrong", "Partly correct", "Not in source / cohort-level"]].sum().sum())
    articles = set(a.pmid) | set(pd.read_csv(SC / "all_B.csv", dtype=str).pmid) | set(pd.read_csv(SC / "all_C.csv", dtype=str).pmid)

    R = {"extraction_accuracy": {
        "final_rows": len(rows), "rows_llm": int((rows.source == "llm").sum()), "rows_added_on_review": int((rows.source == "review").sum()),
        "rows_added_pct": round(100 * (rows.source == "review").mean(), 1), "rows_added_pct_by_gene": byg,
        "field_values": tot_n, "field_changed": tot_e, "field_unchanged_pct": float(t.iloc[-1]["Unchanged, %"]),
        "field_unchanged_pct_by_field": dict(zip(t.Field, t["Unchanged, %"])),
        "nonempty_changed_pct": round(100 * (recs[~recs.llm_empty].cls == "error").mean(), 1),
        "nonempty_values": int((~recs.llm_empty).sum()), "empty_both": int((recs.llm_empty & recs.final_empty).sum()),
        "nonempty_changed_by_field": dict(zip(ne.Field, ne["Changed, % of non-empty"])),
        "subset_articles": len(articles), "subset_patients_checked": len(pat), "subset_patients_found": len(found),
        "subset_found_pct": round(100 * len(found) / len(pat), 1), "subset_not_in_paper": int((pat.verdict != "found").sum()),
        "subset_field_checks": int(s.Checked.sum()), "subset_field_correct_pct": round(100 * (1 - sub_err / s.Checked.sum()), 1)}}
    res = json.loads((PROJECT / "analysis" / "results.json").read_text())
    res.update(R)
    (PROJECT / "analysis" / "results.json").write_text(json.dumps(res, indent=1, default=str), encoding="utf-8")
    print(t.to_string(index=False))
    print(json.dumps(R, indent=1))


if __name__ == "__main__":
    main()
