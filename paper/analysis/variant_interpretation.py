"""ClinVar status of the atlas variants and case-level evidence for variants not classified pathogenic (Table 2, Tables S12-S14).

Revised after review (R2-M7, R3-1..4/7, R4-M1/M2):
- Unit: distinct genomic allele (GRCh38), so equivalent descriptions of one variant are counted once. Small variants on the MANE
  Select transcript and variants in exons absent from MANE (SCN8A/SCN2A exon 5N, CDKL5 NM_003159 3' exons) are included.
- Matching: allele-level against the frozen ClinVar release 2026-09-29 (scripts/curation/clinvar_allele_match.py).
- Recurrence: counted in unrelated families (one per family) across articles, with de novo families, gnomAD v4 allele frequency and
  ClinVar review status; variants above a maximum credible allele frequency are not offered as case-level evidence.
Writes tables/table2_variant_interpretation.csv, supp_variants_reclassification_candidates.csv, supp_variants_clinvar_all.csv,
supp_variants_author_vs_clinvar.csv and results.json entries.
"""
import json
import re
import sys
import time
import urllib.request
from pathlib import Path

from common import GENES, MANE, PROJECT, TAB, WORKSPACE, family_ids, load_rows, pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts/curation"))
CV = Path(__import__("os").environ.get("CLINVAR_MAP_DIR", WORKSPACE / "runs/review_fixes/clinvar"))
SNAP = WORKSPACE / "runs/paper/clinvar_snapshot"
ORDER = ["Pathogenic/likely pathogenic", "Conflicting", "Uncertain significance", "Benign/likely benign", "Other", "Not in ClinVar"]
ALT_OK = {"SCN8A": "NM_014191.4", "CDKL5": "NM_003159.3", "SCN2A": "NM_021007.3"}
HGVS_C = re.compile(r"^c\.[-*]?\d+(?:[+-]\d+)?(?:_[-*]?\d+(?:[+-]\d+)?)?(?:[ACGT]>[ACGT]|del[ACGT]*(?:ins[ACGT]+)?|dup[ACGT]*|ins[ACGT]+|delins[ACGT]+)$")
MAX_AF = 1e-4  # maximum credible population allele frequency for a dominant early-onset DEE allele (conservative)
CHR = {"NC_000002.12": "2", "NC_000012.12": "12", "NC_000020.11": "20", "NC_000023.11": "X", "NC_000006.12": "6", "NC_000009.12": "9",
       "NC_000004.12": "4", "NC_000014.9": "14", "NC_000022.11": "22", "NC_000001.11": "1", "NC_000005.10": "5"}


def stars(review):
    r = str(review).lower()
    if "practice guideline" in r:
        return 4
    if "expert panel" in r:
        return 3
    if "multiple submitters" in r and "no conflict" in r:
        return 2
    if "criteria provided" in r:
        return 1
    return 0


def ref_base_lookup():
    """Reference windows cached by clinvar_allele_match.py: (acc) -> list of (start, seq)."""
    win = {}
    for f in (CV / "ref").glob("*.txt"):
        acc, start, _ = f.stem.rsplit("_", 2)
        win.setdefault(acc, []).append((int(start), f.read_text()))

    def base(acc, pos):
        for s, seq in win.get(acc, []):
            if s <= pos < s + len(seq):
                return seq[pos - s]
        return None
    return base


def gnomad_af(keys, base):
    """gnomAD v4 joint allele frequency for VCF-style alleles (cached); None = not in gnomAD."""
    cache_f = CV / "gnomad_cache.json"
    cache = json.loads(cache_f.read_text()) if cache_f.exists() else {}
    out = {}
    for k in keys:
        acc, pos, ref, alt = k
        pos = int(pos)
        if not ref or not alt:  # indel: add the anchor base (VCF convention)
            a = base(acc, pos - 1)
            if a is None:
                out[k] = None
                continue
            vid = f"{CHR.get(acc, '?')}-{pos - 1}-{a + ref}-{a + alt}"
        else:
            vid = f"{CHR.get(acc, '?')}-{pos}-{ref}-{alt}"
        if cache.get(vid, "missing") is None or vid not in cache:
            q = json.dumps({"query": "{ variant(variantId: \"%s\", dataset: gnomad_r4) { joint { ac an } } }" % vid}).encode()
            val = None
            for attempt in range(6):
                req = urllib.request.Request("https://gnomad.broadinstitute.org/api", data=q, headers={"Content-Type": "application/json"})
                try:
                    with urllib.request.urlopen(req, timeout=60) as r:
                        j = json.loads(r.read().decode())
                    errs = " ".join(e.get("message", "") for e in j.get("errors", []) or [])
                    v = (j.get("data") or {}).get("variant")
                    if v and v.get("joint") and v["joint"]["an"]:
                        val = v["joint"]["ac"] / v["joint"]["an"]
                        break
                    if "not found" in errs.lower() or (v is not None):
                        val = 0.0
                        break
                except Exception:
                    pass
                time.sleep(5 * (attempt + 1))  # rate limit: back off
            cache[vid] = val
            time.sleep(1.0)
        out[k] = cache[vid]
    cache_f.write_text(json.dumps(cache))
    return out


def main():
    am = pd.read_csv(CV / "allele_map.csv", dtype=str, keep_default_na=False)
    amap = {(g, h): r for g, h, r in zip(am.gene, am.hgvs, am.to_dict("records"))}
    rows = load_rows()
    fid, _ = family_ids(rows)
    rows["family_id"] = fid.values
    all_patients = rows.groupby("gene").canonical_patient_id.nunique()
    rows["hgvs"] = [f"{t}:{re.sub(r'[ ]', '', c)}" if HGVS_C.match(re.sub(r"\s", "", str(c))) and
                    (t.split(".")[0] == MANE[g].split(".")[0] or t == ALT_OK.get(g)) else ""
                    for g, t, c in zip(rows.gene, rows.gnomad_transcript, rows.variant_c)]
    r = rows[rows.hgvs != ""].copy()
    # genomic allele when placed on GRCh38; otherwise the transcript HGVS stands for the variant (intronic fallback)
    r["allele"] = [("" if (g, h) not in amap or amap[(g, h)]["status"] == "Unmapped" else
                    "|".join([amap[(g, h)]["chrom"], amap[(g, h)]["pos"], amap[(g, h)]["ref"], amap[(g, h)]["alt"]])
                    if amap[(g, h)]["chrom"] else "hgvs|" + h) for g, h in zip(r.gene, r.hgvs)]
    unmapped = r[r.allele == ""]
    r = r[r.allele != ""].copy()
    for f in ["status", "significance", "review_status", "vcv"]:
        r[f] = [amap[(g, h)][f] for g, h in zip(r.gene, r.hgvs)]
    # inheritance (structured reading) for de novo families
    import comorbidity
    inh = comorbidity.raw().set_index(["gene", "canonical_patient_id"]).inheritance.to_dict()
    r["inheritance"] = [inh.get((g, c), "") for g, c in zip(r.gene, r.canonical_patient_id)]
    v = r.groupby(["gene", "allele"]).agg(
        hgvs=("hgvs", "first"), protein=("protein", "first"), consequence=("consequence", "first"),
        patients=("canonical_patient_id", "nunique"), families=("family_id", "nunique"), articles=("pmid", "nunique"),
        de_novo_families=("family_id", lambda s: s[r.loc[s.index, "inheritance"] == "De novo"].nunique()),
        status=("status", "first"), significance=("significance", "first"), review_status=("review_status", "first"),
        vcv=("vcv", "first"),
        paper_class=("variant_classification_in_paper", lambda s: "; ".join(sorted(set(x for x in s if x))))).reset_index()
    v["stars"] = v.review_status.map(stars)
    # Table 2
    t = []
    for g in GENES + ["All genes"]:
        x = v if g == "All genes" else v[v.gene == g]
        pa = r if g == "All genes" else r[r.gene == g]
        npat = all_patients.sum() if g == "All genes" else all_patients[g]
        row = {"Gene": g, "Patients analysed, n (% of gene)": f"{pa.canonical_patient_id.astype(str).add(pa.gene).nunique()} "
                                                              f"({100 * pa.canonical_patient_id.astype(str).add(pa.gene).nunique() / npat:.0f}%)",
               "Distinct variants, n": len(x)}
        for s in ORDER:
            row[s] = f"{(x.status == s).sum()} ({100 * (x.status == s).mean():.0f}%)"
        row["P/LP with ≥2-star review"] = int(((x.status == ORDER[0]) & (x.stars >= 2)).sum())
        t.append(row)
    t = pd.DataFrame(t)
    t.to_csv(TAB / "table2_variant_interpretation.csv", index=False)
    # population frequency for variants that could be offered as case-level evidence
    base = ref_base_lookup()
    pool = v[v.status.isin(["Not in ClinVar", "Uncertain significance", "Conflicting"]) & (v.families >= 2) & (v.articles >= 2)].copy()
    keys = [tuple(a.split("|")) for a in pool.allele]
    af = gnomad_af([k for k in keys if k[0] != "hgvs"], base)
    pool["gnomad_af"] = [af.get(k) for k in keys]
    pool["too_common"] = pool.gnomad_af.map(lambda a: a is not None and a > MAX_AF)
    pool_af_unknown = int(pool.gnomad_af.isna().sum())
    cand = pool[~pool.too_common].sort_values(["gene", "families"], ascending=[True, False])
    cand[["gene", "hgvs", "protein", "consequence", "families", "patients", "articles", "de_novo_families", "status", "significance",
          "review_status", "vcv", "gnomad_af", "paper_class"]].to_csv(TAB / "supp_variants_reclassification_candidates.csv", index=False)
    v.to_csv(TAB / "supp_variants_clinvar_all.csv", index=False)
    # article classification vs ClinVar (agreement among variants that have a ClinVar classification)
    pc = v.paper_class.str.lower()
    v["paper_group"] = ["VUS" if re.search(r"uncertain|\bvus\b", x) else "P/LP" if re.search(r"pathogenic", x) and "benign" not in x
                        else "Not classified" for x in pc]
    ax = pd.crosstab(v.paper_group, v.status).reindex(index=["P/LP", "VUS", "Not classified"], columns=ORDER, fill_value=0)
    ax.to_csv(TAB / "supp_variants_author_vs_clinvar.csv")
    a_plp = v[(v.paper_group == "P/LP") & v.status.isin(ORDER[:4])]
    a_vus = v[(v.paper_group == "VUS") & v.status.isin(ORDER[:4])]
    pats = r.drop_duplicates(["gene", "canonical_patient_id", "allele"])
    blb = v[v.status == "Benign/likely benign"]
    R = {"variant_interpretation": {
        "clinvar_release": (SNAP / "RELEASE").read_text().strip(), "matching": "allele-level (GRCh38)",
        "unique_variants": len(v), "patients": int(pats.canonical_patient_id.astype(str).add(pats.gene).nunique()),
        "patients_all": int(all_patients.sum()), "unmapped_variants": int(unmapped.hgvs.nunique()),
        "coverage_pct": round(100 * pats.canonical_patient_id.astype(str).add(pats.gene).nunique() / all_patients.sum(), 1),
        "coverage_by_gene": {g: round(100 * pats[pats.gene == g].canonical_patient_id.nunique() / all_patients[g], 1) for g in GENES},
        "status_n": v.status.value_counts().reindex(ORDER, fill_value=0).to_dict(),
        "status_pct": (100 * v.status.value_counts(normalize=True)).round(1).reindex(ORDER, fill_value=0).to_dict(),
        "absent_n": int((v.status == "Not in ClinVar").sum()), "absent_pct": round(100 * (v.status == "Not in ClinVar").mean(), 1),
        "absent_by_gene_pct": {g: round(100 * (v[v.gene == g].status == "Not in ClinVar").mean(), 1) for g in GENES},
        "plp_2star": int(((v.status == ORDER[0]) & (v.stars >= 2)).sum()),
        "patients_with_not_plp_variant_pct": round(100 * (pats.status != ORDER[0]).mean(), 1),
        "single_family_pct": round(100 * (v.families == 1).mean(), 1),
        "pool_n": len(pool), "pool_too_common": int(pool.too_common.sum()), "pool_af_unknown": pool_af_unknown,
        "candidates_n": len(cand), "candidates_families": int(cand.families.sum()),
        "candidates_by_status": cand.status.value_counts().to_dict(),
        "candidates_absent_n": int((cand.status == "Not in ClinVar").sum()),
        "candidates_vus_conf_n": int(cand.status.isin(["Uncertain significance", "Conflicting"]).sum()),
        "candidates_with_de_novo": int((cand.de_novo_families >= 1).sum()),
        "candidates_2plus_de_novo": int((cand.de_novo_families >= 2).sum()),
        "candidates_max": cand.sort_values("families", ascending=False).head(5)[["gene", "hgvs", "protein", "families", "de_novo_families", "status"]].to_dict(orient="records"),
        "too_common_examples": pool[pool.too_common].sort_values("families", ascending=False).head(5)[["gene", "protein", "families", "gnomad_af"]].to_dict(orient="records"),
        "blb_n": len(blb), "blb_patients": int(blb.patients.sum()),
        "author_plp_with_clinvar_class": len(a_plp),
        "author_plp_agree_pct": round(100 * (a_plp.status == ORDER[0]).mean(), 1) if len(a_plp) else None,
        "author_plp_blb_n": int((a_plp.status == "Benign/likely benign").sum()),
        "author_vus_with_clinvar_class": len(a_vus), "author_vus_now_plp": int((a_vus.status == ORDER[0]).sum()),
        "author_plp_absent_n": int(((v.paper_group == "P/LP") & (v.status == "Not in ClinVar")).sum()),
        "author_plp_n": int((v.paper_group == "P/LP").sum()),
    }}
    # copy-number variants whose extent could not be determined (kept: exclusion only when another epilepsy gene is shown)
    cnv_dir = WORKSPACE / "runs/structured_phenotypes/issue_check"
    dec = pd.concat([pd.read_csv(f, dtype=str, keep_default_na=False) for f in
                     sorted(cnv_dir.glob("sweep/decisions/*.csv")) + sorted(cnv_dir.glob("final/decisions/cnv_*.csv"))])
    dec = dec.drop_duplicates(["gene", "mention_id"], keep="last")
    live = set(zip(rows.gene, rows.mention_id.astype(str)))
    unres = dec[(dec.verdict == "unresolved_kept") & [(g, m) in live for g, m in zip(dec.gene, dec.mention_id)]]
    R["variant_interpretation"]["cnv_rows"] = int((rows.consequence == "CNV / exon-level").sum())
    R["variant_interpretation"]["cnv_extent_unknown"] = int(len(unres[unres.gene != "EEF1A2"]))
    R["variant_interpretation"]["cnv_excluded_multigene"] = int((dec.verdict.str.lower() == "excluded").sum())
    res = json.loads((PROJECT / "analysis" / "results.json").read_text())
    res.update(R)
    (PROJECT / "analysis" / "results.json").write_text(json.dumps(res, indent=1, default=str), encoding="utf-8")
    print(t.to_string(index=False))
    print(json.dumps(R, indent=1, default=str)[:3000])


if __name__ == "__main__":
    main()
