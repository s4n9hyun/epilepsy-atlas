"""Main analyses, tables and figures for the DEE atlas manuscript.

Usage: python analysis/run_analysis.py
Writes tables/*.csv, figures/*.png|pdf and analysis/results.json (numbers quoted in the manuscript).
"""
import json
import sqlite3

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from scipy import stats  # noqa: E402

from common import (family_ids, verified_labels, fill_protein, FIG, GENES, MANE, PHENO_ORDER, TAB, WORKSPACE, PROJECT, fetch_cds, is_missing, load_mentions,  # noqa: E402
                    load_patients, load_rows, pd, re)
from treatment import code_treatments  # noqa: E402
import comorbidity  # noqa: E402

plt.rcParams.update({"font.family": "Arial", "font.size": 8, "axes.spines.top": False, "axes.spines.right": False,
                     "pdf.fonttype": 42, "svg.fonttype": "none"})
GENE_COLOR = {"SCN1A": "#2C5AA0", "SCN2A": "#5B9BD5", "SCN8A": "#9DC3E6", "KCNQ2": "#E6852C", "KCNT1": "#F4B183",
              "CDKL5": "#C2408F", "PCDH19": "#E48DC0", "STXBP1": "#3A9E6F", "SYNGAP1": "#7CC49A", "GABRA1": "#B8860B",
              "FOXG1": "#8C6BB1", "SMC1A": "#B9A3D6", "EEF1A2": "#7F7F7F"}
PHENO_COLOR = {"DEE": "#B2182B", "Self-limited (familial) neonatal/infantile epilepsy": "#4393C3",
               "GEFS+ / febrile spectrum": "#92C5DE", "Generalized epilepsy": "#5AAE61", "Other epilepsy": "#BDBDBD",
               "No seizures": "#FDDBC7", "Not specified": "#F0F0F0"}
OUTCOME_COLOR = {"Benefit": "#1A9850", "No benefit": "#BDBDBD", "Worsening": "#D73027"}
R = {}
DBC_GENES = ["SCN1A", "SCN2A", "KCNQ2", "STXBP1"]  # genes with >= 80 phenotyped carriers of both missense and truncating variants


def save(fig, name):
    for ext in ("png", "pdf"):
        fig.savefig(FIG / f"{name}.{ext}", dpi=300, bbox_inches="tight")
    plt.close(fig)


def pct(a, b):
    return round(100 * a / b, 1) if b else None


def fisher(a, b, c, d):
    """2x2 [[a,b],[c,d]] odds ratio with Woolf 95% CI (Haldane correction if a zero cell)."""
    orr, pv = stats.fisher_exact([[a, b], [c, d]])
    aa, bb, cc, dd = (x + 0.5 for x in (a, b, c, d)) if 0 in (a, b, c, d) else (a, b, c, d)
    lo, hi = np.exp(np.log(aa * dd / (bb * cc)) + np.array([-1.96, 1.96]) * np.sqrt(1 / aa + 1 / bb + 1 / cc + 1 / dd))
    return {"OR": round(aa * dd / (bb * cc), 2), "CI": [round(lo, 2), round(hi, 2)], "p": float(f"{pv:.2g}")}


def wilson(k, n):
    """95% Wilson score interval, in percent."""
    if not n:
        return [None, None]
    z = 1.96
    ph = k / n
    den = 1 + z * z / n
    c = (ph + z * z / (2 * n)) / den
    h = z * np.sqrt(ph * (1 - ph) / n + z * z / (4 * n * n)) / den
    return [round(100 * (c - h), 1), round(100 * (c + h), 1)]


def boot_ci(d, ok, reps=2000, seed=1):
    """95% cluster-bootstrap interval (articles resampled) for the proportion of rows where ok is True, in percent."""
    rng = np.random.default_rng(seed)
    g = d.assign(ok=ok.values).groupby("pmid").ok.agg(["sum", "count"])
    if len(g) < 2:
        return [None, None]
    vals = []
    for _ in range(reps):
        idx = rng.integers(0, len(g), len(g))
        x = g.iloc[idx]
        vals.append(100 * x["sum"].sum() / x["count"].sum())
    return [round(float(np.percentile(vals, 2.5)), 1), round(float(np.percentile(vals, 97.5)), 1)]


def cluster_or(d):
    """Logistic regression of y on trunc with cluster-robust (article) standard errors = GEE with independence working
    correlation; point estimate equals the crude OR. Haldane-corrected crude OR is used when a cell is zero."""
    import statsmodels.api as sm
    try:
        X = sm.add_constant(d.trunc.astype(float))
        m = sm.GLM(d.y.astype(float), X, family=sm.families.Binomial()).fit(cov_type="cluster", cov_kwds={"groups": pd.factorize(d.pmid)[0]})
        b, se = m.params["trunc"], m.bse["trunc"]
        if not np.isfinite(se) or se > 10:
            raise ValueError
        return {"OR": round(float(np.exp(b)), 2), "CI": [round(float(np.exp(b - 1.96 * se)), 2), round(float(np.exp(b + 1.96 * se)), 2)],
                "articles": int(d.pmid.nunique())}
    except Exception:
        return {"OR": None, "CI": [0.01, 300.0], "articles": int(d.pmid.nunique())}


# ------------------------------------------------------------------ data
def load_analysis_data():
    global rows, pts, mentions, protein_len
    rows = load_rows()
    _fid, _prob = family_ids(rows)
    rows["family_id"] = _fid.values
    pts = load_patients(rows)
    _pf, _pp = family_ids(pts)
    pts["family_id"], pts["is_proband"] = _pf.values, _pp.values
    mentions = load_mentions()
    protein_len = {g: len(fetch_cds(MANE[g])) // 3 - 1 for g in GENES}

    used = {(g, x) for g, s in zip(rows.gene, rows.supporting_mention_ids) for x in re.split(r"[;,| ]+", str(s)) if x}
    mentions["in_canonical"] = [(g, i) in used for g, i in zip(mentions.gene, mentions.mention_id)]



# ------------------------------------------------------------------ literature flow and Table 1
def article_counts():
    out = {}
    for g in GENES:
        con = sqlite3.connect(WORKSPACE / "databases" / f"{g}.sqlite")
        screened = con.execute("select count(*) from pubmed_articles").fetchone()[0]
        processed = con.execute("select count(*) from processed_articles where status='success'").fetchone()[0]
        with_indiv = con.execute("select count(distinct pmid) from patient_mentions where mention_source_type in "
                                 "('original_patient','manual_exclude','cited_patient','previously_reported_patient','reanalyzed_previous_patient','unclear')").fetchone()[0]
        out[g] = {"screened": screened, "processed": processed, "with_individual_patients": with_indiv,
                  "contributing": int(rows.loc[rows.gene == g, "pmid"].nunique())}
    return out


INDIV = ["original_patient", "manual_exclude", "cited_patient", "previously_reported_patient", "reanalyzed_previous_patient", "unclear"]



def merged_split(g):
    """Original descriptions merged into another atlas patient: same article (within) vs another article (cross)."""
    m = mentions[(mentions.gene == g) & (mentions.mention_source_type == "original_patient")].set_index("mention_id").pmid.astype(str).to_dict()
    within = cross = 0
    for (cp, pm), grp in rows[rows.gene == g].groupby(["canonical_patient_id", "pmid"]):
        ids = {x for s_ in grp.supporting_mention_ids for x in re.split(r"[;,| ]+", str(s_)) if x and x in m}
        if len(ids) > 1:
            pms = [m[x] for x in ids]
            first = str(pm)
            extra = len(ids) - 1
            c = sum(1 for x in pms if x != first)
            cross += min(c, extra)
            within += extra - min(c, extra)
    return within, cross


def flow():
    art = article_counts()
    flow = []
    for g in GENES:
        m = mentions[mentions.gene == g]
        ind = m[m.mention_source_type.isin(INDIV)]
        orig = m[m.mention_source_type.isin(["original_patient", "manual_exclude"])]  # first reports, including those excluded on review
        rereport = ~ind.mention_source_type.isin(["original_patient", "manual_exclude", "unclear"])
        n_pat = int(pts[pts.gene == g].shape[0])
        flow.append({"gene": g, **art[g], "individual_descriptions": len(ind), "aggregate_only": int((m.mention_source_type == "aggregate_only").sum()),
                     "flagged_rereport_at_extraction": int(rereport.sum()), "unclear_originality": int((ind.mention_source_type == "unclear").sum()),
                     "original_descriptions": len(orig), "original_excluded_insufficient_or_manual": int((~orig.in_canonical).sum()),
                     "original_in_canonical": int(orig.in_canonical.sum()), "unique_patients": n_pat,
                     "merged_cross_paper": int(orig.in_canonical.sum()) - int(rows[rows.gene == g].shape[0])
                     + int(rows[rows.gene == g].shape[0]) - n_pat,
                     "merged_within_paper": merged_split(g)[0], "merged_other_paper": merged_split(g)[1],
                     "variant_rows": int(rows[rows.gene == g].shape[0])})
    df = pd.DataFrame(flow)
    tot = df.drop(columns="gene").sum()
    df = pd.concat([df, pd.DataFrame([{"gene": "Total", **tot.to_dict()}])], ignore_index=True)
    df.to_csv(TAB / "flow_counts.csv", index=False)
    R["flow"] = df.set_index("gene").to_dict(orient="index")
    return df


def table1():
    """Per-gene summary (genes in rows). DEE among patients with a phenotype; consequence per patient-variant row; ClinVar status
    per patient-variant row from the allele-level match to the 2026-09-29 release (small variants only; reviewer comments R4-M2/M8)."""
    am_f = WORKSPACE / "runs/review_fixes/clinvar/allele_map.csv"
    am = pd.read_csv(am_f, dtype=str, keep_default_na=False) if am_f.exists() else pd.DataFrame(columns=["gene", "hgvs", "status"])
    cvs = dict(zip(zip(am.gene, am.hgvs), am.status))
    rows_cv = [cvs.get((g, f"{t}:{re.sub(r'[ ]', '', str(c))}"), "") for g, t, c in zip(rows.gene, rows.gnomad_transcript, rows.variant_c)]
    recs = []
    for g in GENES + ["Total"]:
        p = pts if g == "Total" else pts[pts.gene == g]
        m = (rows.gene == g) if g != "Total" else pd.Series(True, index=rows.index)
        r = rows[m]
        cv = pd.Series(rows_cv, index=rows.index)[m]
        known_sex = p[p.sex_h != "NS"]
        onset = p.onset_months.dropna()
        ph = p[p.phenotype_h != "Not specified"]
        nd = r[~r.variant_c.map(is_missing)].groupby("gene").variant_c.apply(lambda x: x.str.replace(r"\s", "", regex=True).str.lower().nunique()).sum()
        def pc(k, n):
            return f"{k} ({100 * k / n:.1f})" if n else ""
        rec = {
            "Gene": g, "Articles": r.pmid.nunique(), "Patients": len(p), "Patient-variant rows": len(r), "Distinct variants": int(nd),
            "Female, n/N (%)": f"{(known_sex.sex_h == 'Female').sum()}/{len(known_sex)} ({pct((known_sex.sex_h == 'Female').sum(), len(known_sex))})",
            "Onset, median (IQR), mo; n": f"{onset.median():.1f} ({onset.quantile(.25):.1f}–{onset.quantile(.75):.1f}); {len(onset)}" if len(onset) else "",
            "DEE label, n/N (%)": f"{(ph.phenotype_h == 'DEE').sum()}/{len(ph)} ({pct((ph.phenotype_h == 'DEE').sum(), len(ph))})",
            "Missense, n (%)": pc((r.consequence == "Missense").sum(), len(r)),
            "Truncating, n (%)": pc((r.consequence == "Truncating").sum(), len(r)),
            "CNV/exon-level, n (%)": pc((r.consequence == "CNV / exon-level").sum(), len(r)),
            "ClinVar P/LP, n (%)": pc((cv == "Pathogenic/likely pathogenic").sum(), len(r)),
            "ClinVar VUS/conflicting, n (%)": pc(cv.isin(["Uncertain significance", "Conflicting"]).sum(), len(r)),
            "ClinVar B/LB, n (%)": pc((cv == "Benign/likely benign").sum(), len(r)),
        }
        recs.append(rec)
    t = pd.DataFrame(recs)
    t.to_csv(TAB / "table1_gene_summary.csv", index=False)
    comp = []
    fields = {"Sex": "sex_h", "Onset age (parsed)": "onset_months", "Seizure presence": "seizure_presence",
              "Phenotype category": "phenotype_review_category", "Seizure types": "seizure_types",
              "Age at last follow-up": "last_followup_age", "ASM response": "response_to_asm", "ASMs used": "asm_used",
              "Targeted treatment response": "targeted_treatment_response", "Protein change (reported or translated)": "protein"}
    for g in GENES + ["Total"]:
        p = pts if g == "Total" else pts[pts.gene == g]
        rec = {"Gene": g, "N": len(p)}
        for label, col in fields.items():
            if col == "onset_months":
                n = p[col].notna().sum()
            elif col == "sex_h":
                n = (p[col] != "NS").sum()
            else:
                n = (~p[col].map(is_missing)).sum()
            rec[label] = f"{pct(n, len(p))}"
        comp.append(rec)
    pd.DataFrame(comp).to_csv(TAB / "supp_field_completeness.csv", index=False)
    R["table1"] = t.set_index("Gene").to_dict(orient="index")
    R["field_completeness_pct"] = pd.DataFrame(comp).set_index("Gene").to_dict(orient="index")
    R["onset_quartiles_months"] = {g: [round(float(q), 3) for q in pts[pts.gene == g].onset_months.dropna().quantile([.25, .5, .75])]
                                   for g in GENES}
    return t


# ------------------------------------------------------------------ Figure 1: flow diagram
def fig1(df):
    tot = df.set_index("gene").loc["Total"]
    fig, ax = plt.subplots(figsize=(7.2, 4.6))
    ax.axis("off")
    boxes = [
        (0.5, 0.93, f"PubMed records retrieved ('<gene> AND epilepsy', plus articles added by clinicians)\nn = {tot.screened:,} gene–article records, {len(GENES)} genes"),
        (0.5, 0.77, f"Full text and supplements processed by LLM extraction\nn = {tot.processed:,} records"),
        (0.5, 0.61, f"Individual-level patient descriptions extracted\nn = {tot.individual_descriptions:,} descriptions (from {tot.with_individual_patients:,} records)"),
        (0.5, 0.40, f"Original reports (first article reporting the genotype)\nn = {tot.original_descriptions:,} descriptions"),
        (0.5, 0.20, f"Original reports retained after source verification by AI review agents\nn = {tot.original_in_canonical:,} descriptions"),
        (0.5, 0.04, f"Unique patients after deduplication\nn = {tot.unique_patients:,} ({tot.variant_rows:,} patient-variant rows)"),
    ]
    for x, y, s in boxes:
        ax.text(x, y, s, ha="center", va="center", fontsize=7.5, bbox=dict(boxstyle="round,pad=0.4", fc="#F7F7F7", ec="#555555", lw=0.8))
    for (x1, y1, _), (x2, y2, _) in zip(boxes[:-1], boxes[1:]):
        ax.annotate("", xy=(x2, y2 + 0.055), xytext=(x1, y1 - 0.055), arrowprops=dict(arrowstyle="->", lw=0.8))
    side = [
        (0.61, 0.69, f"Cohort-level summaries only: {tot.aggregate_only:,}"),
        (0.505, 0.50, f"Curated re-report records: {tot.flagged_rereport_at_extraction:,}\n"
                      f"Originality unclear: {tot.unclear_originality:,}"),
        (0.30, 0.30, f"Excluded on source review: {tot.original_excluded_insufficient_or_manual:,}\n(no variant, no patient-level clinical data, unaffected carrier,\nbenign or non-causal variant, not in the article)"),
        (0.12, 0.12, f"Additional records consolidated to patient level: {tot.original_in_canonical - tot.unique_patients:,}\n(cross-article {tot.merged_other_paper:,}; within-article {tot.merged_within_paper:,}, including multiple variants)"),
    ]
    for x, y, s in side:
        ax.text(0.99, y, s, ha="right", va="center", fontsize=6.5, color="#7F0000")
    save(fig, "figure1_flow")


# ------------------------------------------------------------------ Figure 2: redundancy in the literature
def naive_variant_key(c, p):
    # protein change first (unifies cDNA notations and legacy numbering), cDNA when no protein change is given
    m = re.search(r"p\.\(?([A-Z][a-z]{2}\d+(?:[A-Z][a-z]{2}|Ter))\)?$", str(p).strip())
    if m:
        return m.group(1).lower()
    k = re.sub(r"\s", "", str(c)).lower()
    return k[k.index("c."):] if "c." in k else ""


def fig2(df):
    d = df[df.gene != "Total"].set_index("gene")
    ind = d.individual_descriptions
    unique = d.unique_patients
    rere = d.flagged_rereport_at_extraction + (d.original_in_canonical - d.unique_patients)
    excl = d.original_excluded_insufficient_or_manual + d.unclear_originality
    R["redundancy"] = {"total_descriptions": int(ind.sum()), "rereport_total": int(rere.sum()),
                       "rereport_pct": pct(rere.sum(), ind.sum()), "excluded_total": int(excl.sum()),
                       "by_gene_rereport_pct": {g: pct(rere[g], ind[g]) for g in GENES}}

    # recurrent variants: naive (all individual descriptions) vs deduplicated
    m = mentions[mentions.mention_source_type.isin(INDIV)].copy()
    # same key for both counts: protein translated from the MANE cDNA when possible, otherwise the reported protein/cDNA
    def key(g, c, p):
        mane_p = fill_protein(g, str(c).split(":")[-1].strip(), "") if "c." in str(c) else ""
        return naive_variant_key(c, mane_p if naive_variant_key("", mane_p) else p)
    cmap = rows.assign(k=[key(g, c, p) for g, c, p in zip(rows.gene, rows.variant_c, rows.variant_p)])
    m["k"] = [key(g, cg if "c." in str(cg) else c, p) for g, c, cg, p in zip(m.gene, m.variant_c_original, m.variant_c_gnomad, m.variant_p)]
    rec = []
    for g in GENES:
        nv = m[(m.gene == g) & (m.k != "")].k.value_counts()
        dv = cmap[(cmap.gene == g) & (cmap.k != "")].drop_duplicates(["canonical_patient_id", "k"]).k.value_counts()
        for k in dv.index[:5]:
            rec.append({"gene": g, "variant": k, "naive": int(nv.get(k, 0)), "dedup": int(dv[k]),
                        "protein": rows[(rows.gene == g) & (cmap.k == k)].protein.iloc[0]})
    rv = pd.DataFrame(rec)
    rv["inflation"] = (rv.naive / rv.dedup).round(2)
    rv.to_csv(TAB / "supp_recurrent_variants_naive_vs_dedup.csv", index=False)
    R["recurrent_variants"] = rv.to_dict(orient="records")

    # phenotype proportion (DEE) naive vs dedup
    m["ph"] = [__import__("common").harmonize_phenotype(g, c, s) for g, c, s in zip(m.gene, m.phenotype_review_category, m.seizure_presence)]
    ph = []
    for g in GENES:
        mm = m[(m.gene == g) & (m.ph != "Not specified")]
        pp = pts[(pts.gene == g) & (pts.phenotype_h != "Not specified")]
        ph.append({"gene": g, "naive_DEE_pct": pct((mm.ph == "DEE").sum(), len(mm)), "dedup_DEE_pct": pct((pp.phenotype_h == "DEE").sum(), len(pp)),
                   "naive_n": len(mm), "dedup_n": len(pp)})
    ph = pd.DataFrame(ph)
    ph.to_csv(TAB / "supp_phenotype_naive_vs_dedup.csv", index=False)
    R["phenotype_naive_vs_dedup"] = ph.to_dict(orient="records")

    fig, axes = plt.subplots(1, 3, figsize=(7.2, 2.6), gridspec_kw={"width_ratios": [1.1, 1.3, 1]})
    ax = axes[0]
    y = np.arange(len(GENES))
    for i, (vals, lab, col) in enumerate([(unique, "Unique patients", "#2166AC"), (rere, "Re-report/additional records", "#F4A582"),
                                          (excl, "Excluded/unclear", "#D9D9D9")]):
        left = sum([unique, rere, excl][:i]) if i else 0
        ax.barh(y, vals[GENES] / ind[GENES] * 100, left=(left[GENES] / ind[GENES] * 100) if i else 0, color=col, label=lab, height=0.7)
    ax.set_yticks(y, GENES, fontsize=6.5)
    ax.invert_yaxis()
    ax.set_xlabel("% of individual-level descriptions")
    for i, g in enumerate(GENES):
        ax.text(101, i, f"{ind[g]:,}", va="center", fontsize=6)
    ax.legend(fontsize=6, frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.22), ncol=1)
    ax.set_title("A  Patient descriptions in the literature", loc="left", fontsize=8, fontweight="bold")

    ax = axes[1]
    top = rv.sort_values("naive", ascending=False).groupby("gene").head(2).sort_values("naive", ascending=True).tail(10)
    lab = [f"{g} {re.sub(r'([acgt])>([acgt])', lambda m: m.group(0).upper(), v) if is_missing(p) or '_' in str(p) or '?' in str(p) else re.sub(r'^p\.', '', str(p))}"
           for g, p, v in zip(top.gene, top.protein, top.variant)]
    yy = np.arange(len(top))
    ax.barh(yy + 0.2, top.naive, height=0.4, color="#F4A582", label="Naive count")
    ax.barh(yy - 0.2, top.dedup, height=0.4, color="#2166AC", label="Unique patients")
    ax.set_yticks(yy, lab, fontsize=6)
    ax.set_xlabel("Patients carrying variant")
    ax.legend(fontsize=6, frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.22), ncol=2)
    ax.set_title("B  Recurrent variants", loc="left", fontsize=8, fontweight="bold")

    ax = axes[2]
    ph5 = ph[ph.gene != "SYNGAP1"]  # SYNGAP1 has no DEE category in the phenotype schema
    x = np.arange(len(ph5))
    ax.bar(x - 0.2, ph5.naive_DEE_pct, width=0.4, color="#F4A582", label="Naive")
    ax.bar(x + 0.2, ph5.dedup_DEE_pct, width=0.4, color="#2166AC", label="Deduplicated")
    ax.set_xticks(x, ph5.gene, rotation=60, ha="right", fontsize=6)
    ax.set_ylabel("% DEE among phenotyped")
    ax.legend(fontsize=6, frameon=False)
    ax.set_title("C  DEE proportion", loc="left", fontsize=8, fontweight="bold")
    fig.tight_layout()
    save(fig, "figS1_redundancy")


# ------------------------------------------------------------------ Figure 3: cross-gene landscape
def fig3():
    fig = plt.figure(figsize=(7.2, 9.6))
    gs = fig.add_gridspec(3, 2, height_ratios=[1.15, 1.25, 1.9], hspace=0.75, wspace=0.35)

    ax = fig.add_subplot(gs[0, 0])
    data = [np.log10(pts[pts.gene == g].onset_months.dropna().clip(lower=0.03)) for g in GENES]
    parts = ax.violinplot(data, showextrema=False, widths=0.8)
    for body, g in zip(parts["bodies"], GENES):
        body.set_facecolor(GENE_COLOR[g])
        body.set_alpha(0.6)
    for i, d in enumerate(data, 1):
        ax.plot([i - 0.25, i + 0.25], [np.median(d)] * 2, color="k", lw=1)
    ticks = [0.1, 1, 3, 12, 60, 240]
    ax.set_yticks(np.log10(ticks), ["3 d", "1 m", "3 m", "1 y", "5 y", "20 y"])
    ax.set_xticks(range(1, len(GENES) + 1), [f"{g} ({len(d)})" for g, d in zip(GENES, data)], fontsize=5.5, rotation=60, ha="right")
    ax.set_ylabel("Age at seizure onset")
    ax.set_title("A  Age at seizure onset", loc="left", fontsize=8, fontweight="bold")
    R["onset_kruskal_p"] = float(f"{stats.kruskal(*data).pvalue:.2g}")
    R["onset_median_months"] = {g: round(float(10 ** np.median(d)), 2) for g, d in zip(GENES, data)}

    ax = fig.add_subplot(gs[0, 1])
    ct = pd.crosstab(pts.gene, pts.phenotype_h).reindex(index=GENES, columns=PHENO_ORDER, fill_value=0)
    share = ct.div(ct.sum(axis=1), axis=0) * 100
    left = np.zeros(len(GENES))
    for c in PHENO_ORDER:
        ax.barh(range(len(GENES)), share[c], left=left, color=PHENO_COLOR[c], label=c, height=0.7, edgecolor="white", lw=0.3)
        left += share[c].values
    ax.set_yticks(range(len(GENES)), GENES, fontsize=6.5)
    ax.invert_yaxis()
    ax.set_xlabel("% of patients")
    ax.legend(fontsize=5.5, frameon=False, loc="upper center", bbox_to_anchor=(0.45, -0.28), ncol=2)
    ax.set_title("B  Epilepsy phenotype", loc="left", fontsize=8, fontweight="bold")
    ct.to_csv(TAB / "supp_phenotype_by_gene.csv")

    ax = fig.add_subplot(gs[1, 0])
    cons = ["Missense", "Truncating", "In-frame indel", "Other splice/intronic/synonymous", "CNV / exon-level", "Unclassified"]
    ccol = ["#4575B4", "#D73027", "#FDAE61", "#FEE090", "#762A83", "#E0E0E0"]
    ct2 = pd.crosstab(rows.gene, rows.consequence).reindex(index=GENES, columns=cons, fill_value=0)
    sh2 = ct2.div(ct2.sum(axis=1), axis=0) * 100
    left = np.zeros(len(GENES))
    for c, col in zip(cons, ccol):
        ax.barh(range(len(GENES)), sh2[c], left=left, color=col, label=c, height=0.7, edgecolor="white", lw=0.3)
        left += sh2[c].values
    ax.set_yticks(range(len(GENES)), GENES, fontsize=6.5)
    ax.invert_yaxis()
    ax.set_xlabel("% of patient-variant rows")
    ax.legend(fontsize=5.5, frameon=False, loc="upper center", bbox_to_anchor=(0.45, -0.28), ncol=2)
    ax.set_title("C  Variant consequence", loc="left", fontsize=8, fontweight="bold")
    ct2.to_csv(TAB / "supp_consequence_by_gene.csv")

    # d: truncating vs missense, two outcomes, all genes with >=10 carriers of each class (reviewer comments R1-1/2, R2-M4)
    ax = fig.add_subplot(gs[1, 1])
    dev = comorbidity.raw().set_index(["gene", "canonical_patient_id"]).dev_impairment.to_dict()
    pp = pts.assign(dev=[dev.get((g, c), "NS") for g, c in zip(pts.gene, pts.canonical_patient_id)])
    rec, plot = [], []
    for g in GENES:
        d = pp[pp.consequence.isin(["Missense", "Truncating"])]
        d = d[d.gene == g]
        outcomes = {
            "Reported DEE label": d[d.phenotype_h != "Not specified"].assign(y=lambda x: (x.phenotype_h == "DEE").astype(int)),
            "Developmental impairment (patients with epilepsy)": d[(~d.phenotype_h.isin(["No seizures", "Not specified"]))
                                                                   & d.dev.isin(["Yes", "No"])].assign(y=lambda x: (x.dev == "Yes").astype(int)),
        }
        for oname, s in outcomes.items():
            if g == "SYNGAP1" and oname == "Reported DEE label":
                pass
            t, m = s[s.consequence == "Truncating"], s[s.consequence == "Missense"]
            if len(t) < 10 or len(m) < 10:
                continue
            fz = fisher(int(t.y.sum()), len(t) - int(t.y.sum()), int(m.y.sum()), len(m) - int(m.y.sum()))
            cr = cluster_or(s.assign(trunc=(s.consequence == "Truncating").astype(int)))
            rec.append({"gene": g, "outcome": oname, "truncating": f"{int(t.y.sum())}/{len(t)} ({pct(t.y.sum(), len(t))}%)",
                        "missense": f"{int(m.y.sum())}/{len(m)} ({pct(m.y.sum(), len(m))}%)", "OR": fz["OR"], "CI": fz["CI"],
                        "p": fz["p"], "cluster_robust_CI": cr["CI"], "articles": cr["articles"],
                        "t_pct": pct(t.y.sum(), len(t)), "m_pct": pct(m.y.sum(), len(m)), "n_t": len(t), "n_m": len(m)})
    cd = pd.DataFrame(rec)
    cd.to_csv(TAB / "supp_dee_by_consequence.csv", index=False)
    R["dee_by_consequence"] = {f"{r.gene}|{r.outcome}": {"truncating": r.truncating, "missense": r.missense, "OR": r.OR,
                                                         "CI": r.CI, "p": r.p, "cluster_CI": r.cluster_robust_CI}
                               for r in cd.itertuples()}
    genes_d = [g for g in GENES if g in set(cd.gene)]
    for j, (oname, col, mk) in enumerate([("Reported DEE label", "#B2182B", "o"),
                                          ("Developmental impairment (patients with epilepsy)", "#2166AC", "s")]):
        s = cd[cd.outcome == oname].set_index("gene")
        labelled = set() if j == 0 else labelled
        for i, g in enumerate(genes_d):
            if g not in s.index:
                continue
            r = s.loc[g]
            y = i + (j - 0.5) * 0.32
            a, b_, c_, d_ = r.n_t * r.t_pct / 100, r.n_t * (100 - r.t_pct) / 100, r.n_m * r.m_pct / 100, r.n_m * (100 - r.m_pct) / 100
            zero = min(round(a), round(b_), round(c_), round(d_)) == 0
            if zero:  # Haldane-corrected crude OR with Woolf CI; open marker
                a, b_, c_, d_ = (round(v) + 0.5 for v in (a, b_, c_, d_))
                orr = a * d_ / (b_ * c_)
                se = np.sqrt(1 / a + 1 / b_ + 1 / c_ + 1 / d_)
                lo, hi = orr * np.exp(-1.96 * se), orr * np.exp(1.96 * se)
            else:
                orr = a * d_ / (b_ * c_)
                lo, hi = r.cluster_robust_CI if r.cluster_robust_CI != [0.01, 300.0] else r.CI
            L, H = 0.02, 100
            ax.plot([max(lo, L), min(hi, H)], [y, y], color=col, lw=1)
            if hi > H:
                ax.plot(H, y, ">", color=col, ms=3)
            if lo < L:
                ax.plot(L, y, "<", color=col, ms=3)
            ax.plot(min(max(orr, L), H), y, mk, color=col, ms=4, mfc="white" if zero else col,
                    label=None if oname in labelled else (labelled.add(oname) or oname))
    ax.axvline(1, color="#888888", lw=0.6, ls="--")
    ax.set_xscale("log")
    ax.set_xlim(0.02, 100)
    ax.set_xticks([0.1, 1, 10, 100], ["0.1", "1", "10", "100"])
    ax.set_yticks(range(len(genes_d)), genes_d, fontsize=7)
    ax.invert_yaxis()
    ax.set_xlabel("Odds ratio, truncating vs missense (95% CI, article-clustered)", fontsize=7)
    ax.legend(fontsize=6.5, frameon=False, loc="upper center", bbox_to_anchor=(0.45, -0.22), ncol=1)
    ax.set_title("D  Truncating vs missense variants", loc="left", fontsize=8, fontweight="bold")

    # e: missense position density along the protein
    ax = fig.add_subplot(gs[2, :])
    mis = rows[(rows.consequence == "Missense") & rows.aa_pos.notna()].drop_duplicates(["gene", "canonical_patient_id", "aa_pos"])
    for i, g in enumerate(GENES):
        pos = mis[mis.gene == g].aa_pos / protein_len[g] * 100
        h, edges = np.histogram(pos, bins=100, range=(0, 100))
        h = h / max(h.max(), 1)
        ax.fill_between(edges[:-1], i - h * 0.85, i, step="post", color=GENE_COLOR[g], alpha=0.85)
        mg = mis[mis.gene == g].drop_duplicates(["canonical_patient_id", "protein"])
        top = mg.drop_duplicates(["family_id", "protein"]).protein.value_counts().head(2)  # one count per family
        R.setdefault("missense_top", {})[g] = {k: int(v) for k, v in top.items()}
        ax.text(101, i - 0.35, ", ".join(f"{re.sub(r'^p.', '', k)} ({v})" for k, v in top.items()), fontsize=5.5, va="center")
    ax.set_yticks(np.arange(len(GENES)) - 0.35, [f"{g} ({protein_len[g]} aa)" for g in GENES], fontsize=6.5)
    ax.set_ylim(len(GENES) - 0.9, -1.1)
    ax.set_xlim(0, 125)
    ax.set_xticks([0, 25, 50, 75, 100])
    ax.set_xlabel("Missense position (% of protein length; MANE transcript)")
    ax.set_title("E  Distribution of missense variants along each protein (unique patients; top variants, number of families)", loc="left",
                 fontsize=8, fontweight="bold")
    save(fig, "figure2_landscape")

    # sex
    sx = pts[pts.sex_h != "NS"]
    R["female_pct"] = {g: pct((sx[sx.gene == g].sex_h == "Female").sum(), (sx.gene == g).sum()) for g in GENES}
    k1 = pts[pts.gene == "KCNT1"].phenotype_review_category.str.lower()
    R["kcnt1_categories_pct"] = {"EIMFS": pct(k1.str.contains("migrating").sum(), len(k1)),
                                 "SHE/ADNFLE": pct(k1.str.contains("frontal lobe").sum(), len(k1))}
    R["phenotype_share"] = share.round(1).to_dict(orient="index")
    R["consequence_share"] = sh2.round(1).to_dict(orient="index")


# ------------------------------------------------------------------ Figure 4: precision therapy
def fig4():
    t = code_treatments(pts)
    t.to_csv(TAB / "supp_treatment_coded.csv", index=False)
    # primary: narrative statements only (benefit, no benefit and worsening are all observable there);
    # sensitivity: also count drugs listed only as "responsive ASM" (can only add Benefit)
    known_all = t[t.outcome != "Used, response not reported"]
    known = known_all[known_all.evidence == "narrative"]
    # Source-verified pairs replace the automated label (SCN1A/SCN2A/SCN8A sodium channel blockers, KCNT1 quinidine);
    # pairs not confirmed as an individual response drop out.
    ver = verified_labels()
    live = set(zip(pts.gene, pts.canonical_patient_id))  # patients merged away since verification drop out
    ver = ver[[(g, c) in live for g, c in zip(ver.gene, ver.canonical_patient_id)]]
    key = set(zip(ver.gene, ver.canonical_patient_id, ver.treatment))
    known = known[[(g, c, tr) not in key for g, c, tr in zip(known.gene, known.canonical_patient_id, known.treatment)]]
    add = ver[ver.source_label.isin(list(OUTCOME_COLOR))][["gene", "canonical_patient_id", "pmid", "treatment", "source_label"]]
    known = pd.concat([known, add.rename(columns={"source_label": "outcome"}).assign(evidence="source_verified")], ignore_index=True)
    auto = t.set_index(["gene", "canonical_patient_id", "treatment"]).outcome
    ver["auto"] = [auto.get((g, c, tr), "") for g, c, tr in zip(ver.gene, ver.canonical_patient_id, ver.treatment)]
    R["source_verification"] = {f"{tr}|{g}": {"pairs": len(v), "labels": v.source_label.value_counts().to_dict(),
                                               "db_text_wrong": int((v.db_text_correct.str.lower() == "no").sum()),
                                               "auto_agreement_pct": round(100 * (v.auto == v.source_label).mean(), 1)}
                                for (tr, g), v in ver.groupby(["treatment", "gene"])}
    v6 = ver[ver.gene.isin(["SCN1A", "SCN2A"])]
    R["scb_source_verification"] = {"pairs": len(v6), "labels": v6.source_label.value_counts().to_dict(),
                                    "db_text_wrong": int((v6.db_text_correct.str.lower() == "no").sum()),
                                    "auto_agreement_pct": round(100 * (v6.auto == v6.source_label).mean(), 1)}
    R["treatment_sensitivity_with_responsive_field"] = {
        f"{tr}|{g}": {"n": len(d), "benefit_pct": pct((d.outcome == "Benefit").sum(), len(d)),
                      "worsening_pct": pct((d.outcome == "Worsening").sum(), len(d))}
        for (tr, g), d in known_all.groupby(["treatment", "gene"]) if len(d) >= 5}
    summ = pd.crosstab([t.treatment, t.gene], t.outcome)
    summ.to_csv(TAB / "table2_treatment_response_by_gene.csv")
    R["treatment_counts"] = {f"{a}|{b}": v for (a, b), v in summ.to_dict(orient="index").items()}

    meta = pts.set_index(["gene", "canonical_patient_id"])
    scb = known[known.treatment == "Sodium channel blocker"].copy()
    scb["onset"] = [meta.loc[(g, c)].onset_months for g, c in zip(scb.gene, scb.canonical_patient_id)]
    scb["ph"] = [meta.loc[(g, c)].phenotype_h for g, c in zip(scb.gene, scb.canonical_patient_id)]
    scb["cons"] = [meta.loc[(g, c)].consequence for g, c in zip(scb.gene, scb.canonical_patient_id)]

    named = {
        "SCN1A": scb[scb.gene == "SCN1A"],
        "SCN2A onset <3 m": scb[(scb.gene == "SCN2A") & (scb.onset < 3)],
        "SCN2A onset ≥3 m": scb[(scb.gene == "SCN2A") & (scb.onset >= 3)],
        "SCN8A DEE": scb[(scb.gene == "SCN8A") & (scb.ph == "DEE")],
        "SCN8A self-limited": scb[(scb.gene == "SCN8A") & (scb.ph.str.startswith("Self-limited"))],
        "SCN8A other/not specified": scb[(scb.gene == "SCN8A") & ~scb.ph.isin(["DEE"]) & ~scb.ph.str.startswith("Self-limited")],
        "KCNQ2 self-limited": scb[(scb.gene == "KCNQ2") & (scb.ph.str.startswith("Self-limited"))],
        "KCNQ2 DEE": scb[(scb.gene == "KCNQ2") & (scb.ph == "DEE")],
        "KCNT1": scb[scb.gene == "KCNT1"],
        "PCDH19": scb[scb.gene == "PCDH19"],
        "STXBP1": scb[scb.gene == "STXBP1"],
        "CDKL5": scb[scb.gene == "CDKL5"],
        "SYNGAP1": scb[scb.gene == "SYNGAP1"],
    }
    R["scb"] = {name: {"n": len(s)} for name, s in named.items()}
    groups = [(name, s) for name, s in named.items() if len(s) >= 10]  # figure: groups with >= 10 patients
    for name, s in groups:
        n = len(s)
        R["scb"][name] = {"n": n, **{o: f"{(s.outcome == o).sum()} ({pct((s.outcome == o).sum(), n)}%)" for o in OUTCOME_COLOR}}
    # confidence intervals (Wilson and article-cluster bootstrap) for every shown group (reviewer comments R1-3, R2-M5/M8)
    R["scb_ci"] = {}
    for name, s_ in groups:
        R["scb_ci"][name] = {o: {"wilson": wilson(int((s_.outcome == o).sum()), len(s_)), "cluster": boot_ci(s_, s_.outcome == o)}
                             for o in OUTCOME_COLOR}
    # small-study pattern: outcome by number of patients the article contributes to the group
    def by_size(s_, label):
        sz = s_.groupby("pmid").pmid.transform("size")
        out = {}
        for lab, m in [("1-3", sz <= 3), ("4-10", (sz > 3) & (sz <= 10)), (">10", sz > 10)]:
            x = s_[m]
            out[lab] = {"n": len(x), "articles": int(x.pmid.nunique()), **{o: int((x.outcome == o).sum()) for o in OUTCOME_COLOR}}
        return out
    R["small_study"] = {k: by_size(v, k) for k, v in {"SCN1A SCB": named["SCN1A"], "SCN2A SCB": scb[scb.gene == "SCN2A"],
                                                      "SCN8A SCB": scb[scb.gene == "SCN8A"]}.items()}
    fc = rows[(rows.gene == "FOXG1") & (rows.consequence == "CNV / exon-level")].cnv_type.value_counts()
    R["foxg1_cnv_main"] = (f"duplications ({fc.get('gain', 0)} of {fc.sum()})" if fc.get("gain", 0) > fc.get("loss", 0)
                           else f"deletions ({fc.get('loss', 0)} of {fc.sum()})")
    pd.DataFrame([{"Group": k, "Patients per article": sz, "Patients, n": d["n"], "Articles, n": d["articles"],
                   "Benefit, n (%)": f"{d['Benefit']} ({100 * d['Benefit'] / d['n']:.0f}%)" if d["n"] else "",
                   "No benefit, n": d["No benefit"], "Worsening, n (%)": f"{d['Worsening']} ({100 * d['Worsening'] / d['n']:.0f}%)" if d["n"] else ""}
                  for k, v in R["small_study"].items() for sz, d in v.items()]).to_csv(TAB / "supp_small_study.csv", index=False)
    # SCN8A with and without the caregiver registry (PMID 41034159)
    s8r = scb[(scb.gene == "SCN8A") & (scb.pmid != "41034159")]
    R["scb_scn8a_no_registry"] = {"n": len(s8r), **{o: f"{(s8r.outcome == o).sum()} ({pct((s8r.outcome == o).sum(), len(s8r))}%)" for o in OUTCOME_COLOR}}
    reg = scb[(scb.gene == "SCN8A") & (scb.pmid == "41034159")]
    R["scb_scn8a_registry"] = {"n": len(reg), **{o: int((reg.outcome == o).sum()) for o in OUTCOME_COLOR}}
    # SCN2A onset rule without the cohort that proposed it (Wolff 2017, PMID 28379373)
    s2 = scb[scb.gene == "SCN2A"].dropna(subset=["onset"]).assign(trunc=lambda x: (x.onset < 3).astype(int),
                                                                     y=lambda x: (x.outcome == "Benefit").astype(int))
    R["scb_scn2a_onset_or"] = {"all": cluster_or(s2), "without_wolff2017": cluster_or(s2[s2.pmid != "28379373"]),
                               "wolff_pairs": int((s2.pmid == "28379373").sum()), "n": len(s2)}
    # SCN1A: worsening by drug (source-verified pairs list the drugs)
    vs = pd.read_csv(PROJECT / "analysis" / "scb_source_verified.csv", dtype=str, keep_default_na=False)
    vs = vs[(vs.gene == "SCN1A") & vs.source_label.isin(list(OUTCOME_COLOR))]
    vs = vs[[(g, c) in set(zip(pts.gene, pts.canonical_patient_id)) for g, c in zip(vs.gene, vs.canonical_patient_id)]]
    R["scn1a_scb_by_drug"] = {}
    for drug, rx in [("lamotrigine", "lamotrigine"), ("carbamazepine/oxcarbazepine", "carbamazepine|oxcarbazepine"),
                     ("phenytoin", "phenytoin"), ("lacosamide", "lacosamide")]:
        x = vs[vs.drugs.str.lower().str.contains(rx)]
        R["scn1a_scb_by_drug"][drug] = {"n": len(x), **{o: int((x.source_label == o).sum()) for o in OUTCOME_COLOR}}
    early, late = named["SCN2A onset <3 m"], named["SCN2A onset ≥3 m"]
    R["scb_scn2a_early_vs_late_benefit"] = fisher((early.outcome == "Benefit").sum(), (early.outcome != "Benefit").sum(),
                                                  (late.outcome == "Benefit").sum(), (late.outcome != "Benefit").sum())
    R["scb_scn2a_early_vs_late_worsening"] = fisher((late.outcome == "Worsening").sum(), (late.outcome != "Worsening").sum(),
                                                    (early.outcome == "Worsening").sum(), (early.outcome != "Worsening").sum())
    s8 = scb[scb.gene == "SCN8A"]
    R["scb_scn8a_registry_n"] = int((s8.pmid == "41034159").sum())
    R["scb_scn8a_all"] = {"n": len(s8), **{o: f"{(s8.outcome == o).sum()} ({pct((s8.outcome == o).sum(), len(s8))}%)" for o in OUTCOME_COLOR}}
    R["scb_scn8a_vs_scn1a_benefit"] = fisher((s8.outcome == "Benefit").sum(), (s8.outcome != "Benefit").sum(),
                                             (named["SCN1A"].outcome == "Benefit").sum(), (named["SCN1A"].outcome != "Benefit").sum())
    s1 = scb[scb.gene == "SCN1A"]
    R["scb_scn1a_by_phenotype"] = {ph: {"n": int((s1.ph == ph).sum()), "worsening": int(((s1.ph == ph) & (s1.outcome == "Worsening")).sum()),
                                        "benefit": int(((s1.ph == ph) & (s1.outcome == "Benefit")).sum())}
                                   for ph in ["DEE", "GEFS+ / febrile spectrum", "Other epilepsy"]}
    a = s1[s1.ph == "DEE"]
    b = s1[s1.ph != "DEE"]
    R["scb_scn1a_dee_vs_nondee_worsening"] = fisher((a.outcome == "Worsening").sum(), (a.outcome != "Worsening").sum(),
                                                    (b.outcome == "Worsening").sum(), (b.outcome != "Worsening").sum())

    fig, axes = plt.subplots(1, 2, figsize=(7.2, 4.2), gridspec_kw={"width_ratios": [1.2, 1]})
    ax = axes[0]
    verified_genes = {"SCN1A", "SCN2A", "SCN8A"}
    for i, (name, s) in enumerate(groups):
        left = 0
        auto = name.split()[0] not in verified_genes
        for o, col in OUTCOME_COLOR.items():
            v = pct((s.outcome == o).sum(), len(s)) or 0
            ax.barh(i, v, left=left, color=col, height=0.7, label=o if i == 0 else None, hatch="///" if auto else None,
                    edgecolor="white" if auto else None, lw=0)
            left += v
        lo, hi = wilson(int((s.outcome == "Benefit").sum()), len(s))
        ax.plot([lo, hi], [i, i], color="black", lw=0.8)
        ax.text(101, i, f"n={len(s)}", va="center", fontsize=6.5)
    ax.set_yticks(range(len(groups)), [g[0] + (" *" if g[0].split()[0] not in verified_genes else "") for g in groups], fontsize=7)
    ax.invert_yaxis()
    ax.set_xlabel("% of patients with an individually reported response")
    ax.set_xlim(0, 112)
    from matplotlib.patches import Patch
    h, l = ax.get_legend_handles_labels()
    h += [Patch(facecolor="white", edgecolor="#555555", hatch="///"), plt.Line2D([], [], color="black", lw=1)]
    l += ["* Automated coding (not source-verified)", "Wilson 95% CI of benefit"]
    ax.legend(h, l, fontsize=6, frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.18), ncol=3)
    ax.set_title("A  Sodium channel blockers", loc="left", fontsize=8, fontweight="bold")

    ax = axes[1]
    combos = [("Stiripentol", "SCN1A"), ("Fenfluramine", "SCN1A"), ("Cannabidiol", "SCN1A"), ("Ketogenic diet", "SCN1A"),
              ("Ketogenic diet", "SCN2A"), ("Ketogenic diet", "CDKL5"), ("Cannabidiol", "CDKL5"), ("Fenfluramine", "CDKL5"),
              ("Vigabatrin/ACTH/steroids", "CDKL5"), ("Retigabine", "KCNQ2"), ("Quinidine", "KCNT1"), ("Ketogenic diet", "KCNT1"),
              ("Cannabidiol", "KCNT1"), ("Ketogenic diet", "SCN8A"), ("Vigabatrin/ACTH/steroids", "SCN8A"),
              ("Vigabatrin/ACTH/steroids", "STXBP1"), ("Ketogenic diet", "STXBP1"), ("Stiripentol", "PCDH19"),
              ("Vigabatrin/ACTH/steroids", "PCDH19")]
    R["other_treatments"] = {}
    combos = [c for c in combos if ((known.treatment == c[0]) & (known.gene == c[1])).sum() >= 10]
    for i, (tr, g) in enumerate(combos):
        s = known[(known.treatment == tr) & (known.gene == g)]
        left = 0
        auto = not (tr == "Quinidine" and g == "KCNT1")
        for o, col in OUTCOME_COLOR.items():
            v = pct((s.outcome == o).sum(), len(s)) or 0
            ax.barh(i, v, left=left, color=col, height=0.7, hatch="///" if auto else None, edgecolor="white" if auto else None, lw=0)
            left += v
        lo, hi = wilson(int((s.outcome == "Benefit").sum()), len(s))
        ax.plot([lo, hi], [i, i], color="black", lw=0.8)
        ax.text(101, i, f"n={len(s)}", va="center", fontsize=6.5)
        R["other_treatments"][f"{tr}|{g}"] = {"n": len(s), **{o: int((s.outcome == o).sum()) for o in OUTCOME_COLOR}}
    ax.set_yticks(range(len(combos)), [f"{g}: {tr}" + ("" if (tr == "Quinidine" and g == "KCNT1") else " *") for tr, g in combos], fontsize=7)
    ax.invert_yaxis()
    ax.set_xlim(0, 118)
    ax.set_xlabel("% of patients with an individually reported response")
    ax.set_title("B  Other precision / syndrome-directed therapies", loc="left", fontsize=8, fontweight="bold")
    fig.tight_layout()
    save(fig, "figure4_treatment")
    # primary-analysis table behind Figure 4 (replaces the automated crosstab as Table S9)
    rows_t = []
    for name, s_ in [(n_, x) for n_, x in named.items() if len(x)] + [(f"{g}: {tr}", known[(known.treatment == tr) & (known.gene == g)])
                                                                    for tr, g in sorted(set(zip(known.treatment, known.gene)))
                                                                    if tr != "Sodium channel blocker"]:
        if not len(s_):
            continue
        n_ = len(s_)
        verified = (name.split()[0] in verified_genes and "Sodium" not in name and ":" not in name) or name == "KCNT1: Quinidine"
        rec_ = {"Group": name, "n": n_, "Articles": s_.pmid.nunique(), "Labels": "source-verified" if verified else "automated"}
        for o in OUTCOME_COLOR:
            k_ = int((s_.outcome == o).sum())
            lo_, hi_ = wilson(k_, n_)
            rec_[o] = f"{k_} ({100 * k_ / n_:.1f}%; {lo_}–{hi_})"
        rows_t.append(rec_)
    pd.DataFrame(rows_t).to_csv(TAB / "supp_S9_treatment_primary.csv", index=False)
    sens = known_all.groupby(["treatment", "gene"]).outcome.value_counts().unstack(fill_value=0)
    sens.to_csv(TAB / "supp_S9b_treatment_with_responsive_field.csv")


if __name__ == "__main__":
    load_analysis_data()
    for d in (FIG, TAB):
        d.mkdir(exist_ok=True)
    df = flow()
    table1()
    fig1(df)
    fig2(df)
    fig3()
    fig4()
    (PROJECT / "analysis" / "results.json").write_text(json.dumps(R, indent=1, default=lambda o: o.item() if hasattr(o, "item") else str(o)),
                                                       encoding="utf-8")
    print(json.dumps(R, indent=1, default=str)[:6000])
