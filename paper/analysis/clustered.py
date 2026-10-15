"""Sensitivity analysis: key comparisons re-estimated with GEE logistic models clustered by source article.

Patients from the same article are not independent (shared ascertainment and reporting style). Each comparison is
refit as logit(outcome) ~ exposure with an exchangeable working correlation within PMID; robust (sandwich) 95% CI.
Writes tables/supp_S10_clustered_sensitivity.csv and merges results into analysis/results.json.
"""
import json
import warnings

import numpy as np
import statsmodels.api as sm
from scipy import stats

from common import PROJECT, TAB, load_patients, pd, verified_labels
from treatment import code_treatments

warnings.filterwarnings("ignore")


def wilson_ci(k, n):
    z = 1.96
    if not n:
        return [None, None]
    ph = k / n
    den = 1 + z * z / n
    c = (ph + z * z / (2 * n)) / den
    h = z * np.sqrt(ph * (1 - ph) / n + z * z / (4 * n * n)) / den
    return [round(c - h, 3), round(c + h, 3)]


def gee(df, y, x):
    """Crude OR (Woolf CI, Haldane if a zero cell), cluster-robust logistic (GEE independence; patients as unit, article clusters),
    GEE exchangeable (reported only when the estimated within-article correlation is admissible, 0 <= rho < 1) and the
    within-article conditional-logit OR (reviewer comments R2-M4/M8)."""
    d = df[[y, x, "pmid"]].dropna().astype({y: int, x: int}).sort_values("pmid")
    X = sm.add_constant(d[[x]])
    ind = sm.GEE(d[y], X, groups=d["pmid"], family=sm.families.Binomial(), cov_struct=sm.cov_struct.Independence()).fit()
    exch, rho = None, None
    try:
        e = sm.GEE(d[y], X, groups=d["pmid"], family=sm.families.Binomial(), cov_struct=sm.cov_struct.Exchangeable()).fit()
        rho = float(e.model.cov_struct.dep_params)
        if np.isfinite(e.bse[x]) and 0 <= rho < 1:
            exch = e
    except Exception:
        pass
    a = ((d[x] == 1) & (d[y] == 1)).sum(), ((d[x] == 1) & (d[y] == 0)).sum()
    c = ((d[x] == 0) & (d[y] == 1)).sum(), ((d[x] == 0) & (d[y] == 0)).sum()
    orr, p = stats.fisher_exact([[a[0], a[1]], [c[0], c[1]]])
    aa, bb, cc, dd = (v + 0.5 for v in (a[0], a[1], c[0], c[1])) if 0 in (a[0], a[1], c[0], c[1]) else (a[0], a[1], c[0], c[1])
    lor, se0 = np.log(aa * dd / (bb * cc)), np.sqrt(1 / aa + 1 / bb + 1 / cc + 1 / dd)
    clogit = None
    try:
        from statsmodels.discrete.conditional_models import ConditionalLogit
        inf = d.groupby("pmid").filter(lambda g: g[y].nunique() > 1 and g[x].nunique() > 1)
        if len(inf) >= 10:
            cl = ConditionalLogit(inf[y], inf[[x]], groups=inf["pmid"]).fit(disp=0)
            clogit = {"OR": round(float(np.exp(cl.params[x])), 2),
                      "CI": [round(float(np.exp(cl.params[x] - 1.96 * cl.bse[x])), 2), round(float(np.exp(cl.params[x] + 1.96 * cl.bse[x])), 2)],
                      "articles": int(inf.pmid.nunique())}
    except Exception:
        pass
    b, se = ind.params[x], ind.bse[x]
    out = {"n": len(d), "articles": d.pmid.nunique(), "crude_OR": round(float(np.exp(lor)), 2),
           "crude_CI": [round(float(np.exp(lor - 1.96 * se0)), 2), round(float(np.exp(lor + 1.96 * se0)), 2)], "crude_p": float(f"{p:.2g}"),
           "GEE_OR": round(float(np.exp(b)), 2), "GEE_CI": [round(float(np.exp(b - 1.96 * se)), 2), round(float(np.exp(b + 1.96 * se)), 2)],
           "GEE_p": float(f"{ind.pvalues[x]:.2g}"), "working_correlation": "independence (cluster-robust)",
           "exch_OR": round(float(np.exp(exch.params[x])), 2) if exch is not None else None,
           "exch_CI": [round(float(np.exp(exch.params[x] - 1.96 * exch.bse[x])), 2), round(float(np.exp(exch.params[x] + 1.96 * exch.bse[x])), 2)] if exch is not None else None,
           "rho": round(rho, 2) if rho is not None else None, "within_article": clogit}
    return out


def main():
    pts = load_patients()
    out = {}
    ph = pts[pts.phenotype_h != "Not specified"].copy()
    ph["dee"] = (ph.phenotype_h == "DEE").astype(int)
    for g in ["SCN1A", "SCN2A", "SCN8A", "KCNQ2", "STXBP1"]:
        d = ph[(ph.gene == g) & ph.consequence.isin(["Truncating", "Missense"])].copy()
        d["trunc"] = (d.consequence == "Truncating").astype(int)
        out[f"{g}: DEE, truncating vs missense"] = gee(d, "dee", "trunc")
    t = code_treatments(pts)
    t = t[(t.treatment == "Sodium channel blocker") & (t.evidence == "narrative") & (t.outcome != "Used, response not reported")]
    v = verified_labels()  # source-verified labels (same as Figure 4)
    v = v[(v.treatment == "Sodium channel blocker") & v.source_label.isin(["Benefit", "No benefit", "Worsening"])]
    live = set(zip(pts.gene, pts.canonical_patient_id))
    v = v[[(g, c) in live for g, c in zip(v.gene, v.canonical_patient_id)]]
    t = pd.concat([t[~t.gene.isin(set(v.gene))], v[["gene", "canonical_patient_id", "pmid"]].assign(outcome=v.source_label.values)],
                  ignore_index=True)
    meta = pts.set_index(["gene", "canonical_patient_id"])
    t["onset"] = [meta.at[(g, c), "onset_months"] for g, c in zip(t.gene, t.canonical_patient_id)]
    t["ph"] = [meta.at[(g, c), "phenotype_h"] for g, c in zip(t.gene, t.canonical_patient_id)]
    s2 = t[(t.gene == "SCN2A") & t.onset.notna()].copy()
    s2["early"] = (s2.onset < 3).astype(int)
    s2["benefit"] = (s2.outcome == "Benefit").astype(int)
    s2["worse"] = (s2.outcome == "Worsening").astype(int)
    s2["late"] = 1 - s2.early
    out["SCN2A: SCB benefit, onset <3 m vs >=3 m"] = gee(s2, "benefit", "early")
    out["SCN2A: SCB worsening, onset >=3 m vs <3 m"] = gee(s2, "worse", "late")
    s1 = t[t.gene == "SCN1A"].copy()
    s1["dee"] = (s1.ph == "DEE").astype(int)
    s1["worse"] = (s1.outcome == "Worsening").astype(int)
    out["SCN1A: SCB worsening, DEE vs non-DEE"] = gee(s1, "worse", "dee")

    rows = []
    for k, v in out.items():
        w = v["within_article"] or {}
        rows.append({"Comparison": k, "n (articles)": f"{v['n']:,} ({v['articles']})",
                     "Crude OR (95% CI)": f"{v['crude_OR']} ({v['crude_CI'][0]}–{v['crude_CI'][1]})", "Fisher P": v["crude_p"],
                     "Cluster-robust OR (95% CI)": f"{v['GEE_OR']} ({v['GEE_CI'][0]}–{v['GEE_CI'][1]})", "Cluster-robust P": v["GEE_p"],
                     "Exchangeable GEE OR (95% CI); rho": (f"{v['exch_OR']} ({v['exch_CI'][0]}–{v['exch_CI'][1]}); {v['rho']}" if v["exch_OR"] is not None else f"not admissible (rho {v['rho']})"),
                     "Within-article OR (95% CI); articles": (f"{w['OR']} ({w['CI'][0]}–{w['CI'][1]}); {w['articles']}" if w else "—")})
    pd.DataFrame(rows).to_csv(TAB / "supp_S10_clustered_sensitivity.csv", index=False)
    res = json.loads((PROJECT / "analysis" / "results.json").read_text())
    res["clustered_sensitivity"] = out
    (PROJECT / "analysis" / "results.json").write_text(json.dumps(res, indent=1, default=str), encoding="utf-8")
    print(pd.DataFrame(rows).to_string())


if __name__ == "__main__":
    main()
