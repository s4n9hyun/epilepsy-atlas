"""Comorbidities and inheritance by gene, read from the source articles (runs/structured_phenotypes), and Figure 3
(A seizure types, B comorbidities, C inheritance).

Every unique patient was looked up in the original article and the fields were coded from the source (Yes / No / NS).
Comorbidity percentages are reported frequencies: patients in whom the feature was reported, over all patients of the gene
(features are often not mentioned, so these are lower bounds); developmental impairment is also given among patients with
a statement either way. Inheritance is summarised among patients whose parents were tested (de novo, maternal, paternal or
inherited from an unspecified parent).
Writes tables/supp_comorbidity_by_gene.csv, tables/supp_inheritance_by_gene.csv, figures/figure3_clinical_profile.*,
and merges results into results.json.
"""
import json

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from common import FIG, GENES, PROJECT, TAB, WORKSPACE, family_ids, load_patients, pd, re  # noqa: E402

plt.rcParams.update({"font.family": "Arial", "font.size": 7.5, "axes.spines.top": False, "axes.spines.right": False,
                     "pdf.fonttype": 42, "svg.fonttype": "none"})
SRC = WORKSPACE / "runs/structured_phenotypes/results"
FEATURES = [("dev_impairment", "Developmental impairment"), ("sev_severe", "Severe/profound impairment"), ("autism", "Autism"),
            ("movement_disorder", "Movement disorder"), ("microcephaly", "Microcephaly"), ("cvi", "Cortical visual impairment"),
            ("hypotonia", "Hypotonia"), ("deceased", "Deceased")]
INH = ["De novo", "Maternal", "Paternal", "Inherited, parent not specified"]
INH_COLOR = {"De novo": "#3B6FB6", "Maternal": "#E07B39", "Paternal": "#F2B880", "Inherited, parent not specified": "#C9C9C9"}


FIX = WORKSPACE / "runs/structured_phenotypes/issue_check/structured_fix"
VAL = WORKSPACE / "runs/structured_phenotypes/validation"
ADD = WORKSPACE / "runs/structured_phenotypes/addendum/results"  # patients added or re-linked after the first reading
VFIELDS = ["dev_impairment", "dev_severity", "autism", "movement_disorder", "microcephaly", "cvi", "hypotonia", "deceased",
           "inheritance", "parental_mosaic", "patient_mosaic"]


def raw():
    """Readings per patient; corrections made during the issue check (structured_fix) replace the first reading."""
    fr = [pd.read_csv(f, dtype=str, keep_default_na=False) for f in sorted(SRC.glob("*.csv"))]
    fr += [pd.read_csv(f, dtype=str, keep_default_na=False) for f in sorted(FIX.glob("*.csv"))] if FIX.exists() else []
    fr += [pd.read_csv(f, dtype=str, keep_default_na=False) for f in sorted(ADD.glob("*.csv"))] if ADD.exists() else []
    return pd.concat(fr, ignore_index=True).drop_duplicates(["gene", "canonical_patient_id"], keep="last")


def validation(s):
    """Blinded second reading of random patients versus the dataset values."""
    files = sorted((VAL / "results").glob("*.csv"))
    if not files:
        return None
    v = pd.concat([pd.read_csv(f, dtype=str, keep_default_na=False) for f in files])
    m = v.merge(s, on=["gene", "canonical_patient_id"], suffixes=("_v", ""))
    m = m[(m.found_v == "1") & (m.found == "1")]
    rows = []
    for f in VFIELDS:
        a, b = m[f + "_v"].replace({"Features": "Yes"} if f == "autism" else {}), m[f].replace({"Features": "Yes"} if f == "autism" else {})
        both = (a != "NS") | (b != "NS")
        pos = {"Yes"} if f not in ("inheritance", "parental_mosaic", "dev_severity") else (
            {"De novo", "Maternal", "Paternal"} if f == "inheritance" else {"Maternal", "Paternal"} if f == "parental_mosaic" else {"Severe", "Profound"})
        pa, pb = a.isin(pos), b.isin(pos)
        ppos = round(100 * 2 * int((pa & pb & (a == b)).sum()) / (int(pa.sum()) + int(pb.sum())), 1) if pa.sum() + pb.sum() else None
        rows.append({"Field": f, "Patients": len(m), "Agreement, %": round(100 * (a == b).mean(), 1),
                     "Agreement where either reading recorded a value, %": round(100 * (a[both] == b[both]).mean(), 1) if both.any() else None,
                     "Positive agreement, %": ppos})
    t = pd.DataFrame(rows)
    t.to_csv(TAB / "supp_comorbidity_validation.csv", index=False)
    vals = sum(len(m) for _ in VFIELDS)
    agree = sum(int(((m[f + "_v"]) == m[f]).sum()) for f in VFIELDS)
    return {"patients": len(m), "agreement_pct": round(100 * agree / vals, 1),
            "agreement_nonmissing_pct": round(float(t["Agreement where either reading recorded a value, %"].mean()), 1),
            "by_field": dict(zip(t.Field, t["Agreement, %"])),
            "positive_agreement_by_field": dict(zip(t.Field, t["Positive agreement, %"])),
            "positive_agreement_min": float(t["Positive agreement, %"].min()), "positive_agreement_max": float(t["Positive agreement, %"].max())}


def load():
    s = raw()
    pts = load_patients()
    fid, prob = family_ids(pts)
    pts = pts.assign(family_id=fid.values, is_proband=prob.values)[["gene", "canonical_patient_id", "pmid", "patient_id", "family_id", "is_proband"]]
    s = pts.merge(s.drop(columns=[c for c in ["pmid", "patient_id"] if c in s.columns]), on=["gene", "canonical_patient_id"], how="left")
    s = s.fillna("NS")
    s["sev_severe"] = np.where(s.dev_severity.isin(["Severe", "Profound"]), "Yes", np.where(s.dev_severity.isin(["Mild", "Moderate"]), "No", "NS"))
    s["autism"] = s.autism.replace({"Features": "Yes"})
    return s


def main():
    s = load()
    covered = (s.found != "NS").mean()
    rows = []
    for g in GENES:
        d = s[s.gene == g]
        r = {"Gene": g, "Patients, n": len(d)}
        for f, lab in FEATURES:
            yes, no = (d[f] == "Yes").sum(), (d[f] == "No").sum()
            r[f"{lab}, n"] = int(yes)
            r[f"{lab}, % of all"] = round(100 * yes / len(d), 3)
            r[f"{lab}, % of stated"] = round(100 * yes / (yes + no), 3) if yes + no else None
            r[f"{lab}, stated n"] = int(yes + no)
        rows.append(r)
    com = pd.DataFrame(rows)
    com.to_csv(TAB / "supp_comorbidity_by_gene.csv", index=False)

    # inheritance per proband (one per family; reviewer comments R3-6, R2-M8.9); per patient as secondary view
    def inh_table(d):
        k = d[d.inheritance.isin(INH)]
        t = pd.crosstab(k.gene, k.inheritance).reindex(index=GENES, columns=INH, fill_value=0)
        return t, t.sum(axis=1)
    prob = s[s.is_proband]
    inh, tested = inh_table(prob)
    inh_p, tested_p = inh_table(s)
    inh_pct = inh.div(tested.replace(0, np.nan), axis=0) * 100
    parent_lab = s.patient_id.astype(str).str.contains(r"mother|father|parent", case=False)
    # mosaic transmitting parents: count families (not every child), and do not count transmitting parents as mosaic patients
    pm = s[s.parental_mosaic.isin(["Maternal", "Paternal"])].drop_duplicates("family_id").groupby("gene").size().reindex(GENES, fill_value=0)
    inherited = inh[["Maternal", "Paternal", "Inherited, parent not specified"]].sum(axis=1)
    out = inh.copy()
    out.insert(0, "Probands with tested parents, n", tested)
    out["De novo, % of probands"] = inh_pct["De novo"].round(1)
    out["De novo, % of patients (all relatives counted)"] = (100 * inh_p["De novo"] / tested_p.replace(0, np.nan)).round(1)
    out["Families with a mosaic transmitting parent, n"] = pm
    out["Mosaic patients (excluding transmitting parents), n"] = s[(s.patient_mosaic == "Yes") & ~parent_lab].groupby("gene").size().reindex(GENES, fill_value=0)
    out["Probands with unknown/not tested inheritance, n"] = prob[~prob.inheritance.isin(INH)].groupby("gene").size().reindex(GENES, fill_value=0)
    out.to_csv(TAB / "supp_inheritance_by_gene.csv")
    # deaths: causes and age at death (reviewer comments R1-7, R4 minor 17)
    def age_months(v):
        t = str(v).lower()
        m = re.search(r"(\d+(?:\.\d+)?)\s*(y|year|m|month|w|week|d|day)", t)
        if not m:
            return None
        x = float(m.group(1))
        return x * {"y": 12, "m": 1, "w": 7 / 30.44, "d": 1 / 30.44}[m.group(2)[0]]
    dd = s[s.deceased == "Yes"].assign(cause=lambda x: x.cause_of_death.str.split(":").str[0].replace({"NS": "Unknown"}),
                                       age=lambda x: x.age_at_death.map(age_months))
    drow = []
    for g in GENES + ["All genes"]:
        x = dd if g == "All genes" else dd[dd.gene == g]
        a = x.age.dropna()
        drow.append({"Gene": g, "Deaths, n": len(x), **{c: int((x.cause == c).sum()) for c in ["SUDEP", "Status epilepticus", "Respiratory/infection", "Other", "Unknown"]},
                     "Age at death known, n": len(a), "Age at death, median (IQR), months": f"{a.median():.0f} ({a.quantile(.25):.0f}–{a.quantile(.75):.0f})" if len(a) else ""})
    pd.DataFrame(drow).to_csv(TAB / "supp_deaths_by_gene.csv", index=False)
    death_age = {c: (round(float(dd[dd.cause == c].age.median()), 0) if dd[dd.cause == c].age.notna().any() else None)
                 for c in ["SUDEP", "Status epilepticus"]}

    # ---------------------------------------------------------------- Figure 3
    sz = pd.read_csv(TAB / "supp_seizure_types_by_gene.csv", index_col=0)
    raw_sz = pd.read_csv(TAB / "_raw_seizure_types_pct.csv", index_col=0)  # unrounded values for the heat map
    sz_cols = [c for c in sz.columns if not c.endswith(", n")]
    fig = plt.figure(figsize=(7.2, 8.6))
    gs = fig.add_gridspec(2, 2, height_ratios=[1, 1], width_ratios=[1.35, 1], hspace=0.42, wspace=0.55)
    axA = fig.add_subplot(gs[0, :])
    axB = fig.add_subplot(gs[1, 0])
    axC = fig.add_subplot(gs[1, 1])

    def heat(ax, mat, xl, yl, vmax, cmap):
        im = ax.imshow(mat, cmap=cmap, vmin=0, vmax=vmax, aspect="auto")
        ax.set_xticks(range(len(xl)))
        ax.set_xticklabels(xl, rotation=40, ha="right")
        ax.set_yticks(range(len(yl)))
        ax.set_yticklabels(yl)
        for i in range(mat.shape[0]):
            for j in range(mat.shape[1]):
                v = mat[i, j]
                if not np.isnan(v):
                    ax.text(j, i, f"{v:.0f}", ha="center", va="center", fontsize=6, color="white" if v >= vmax * 0.6 else "black")
        for sp in ax.spines.values():
            sp.set_visible(False)
        ax.tick_params(length=0)
        return im

    nsz = sz["Patients with a classified seizure type, n"]
    imA = heat(axA, raw_sz.loc[GENES, sz_cols].values, [c.replace("Bilateral tonic-clonic", "Tonic-clonic") for c in sz_cols],
               [f"{g} (n={nsz[g]})" for g in GENES], 80, "Reds")
    fig.colorbar(imA, ax=axA, fraction=0.025, pad=0.01).set_label("% of patients with a classified seizure type")
    axA.set_title("A  Seizure types", loc="left", fontweight="bold", fontsize=8.5)

    bcols = [f"{lab}, % of all" for _, lab in FEATURES]
    imB = heat(axB, com.set_index("Gene").loc[GENES, bcols].values.astype(float), [lab for _, lab in FEATURES],
               [f"{g} (n={len(s[s.gene == g])})" for g in GENES], 80, "Purples")
    fig.colorbar(imB, ax=axB, fraction=0.04, pad=0.02).set_label("% of all patients (reported)")
    axB.set_title("B  Comorbidities and outcome", loc="left", fontweight="bold", fontsize=8.5)

    left = np.zeros(len(GENES))
    y = np.arange(len(GENES))
    for c in INH:
        vals = inh_pct[c].fillna(0).values
        axC.barh(y, vals, left=left, color=INH_COLOR[c], label=c.replace("Inherited, parent not specified", "Inherited, parent n.s."),
                 edgecolor="white", linewidth=0.4)
        left += vals
    for i, g in enumerate(GENES):
        axC.text(101, i, f"{tested[g]}" + (f" ({pm[g]} pm)" if pm[g] else ""), va="center", fontsize=6)
    axC.set_yticks(y)
    axC.set_yticklabels(GENES)
    axC.invert_yaxis()
    axC.set_xlim(0, 100)
    axC.set_xlabel("% of probands with tested parents (one per family)")
    axC.legend(loc="upper center", bbox_to_anchor=(0.45, -0.17), ncol=2, frameon=False, fontsize=6.5)
    axC.set_title("C  Inheritance", loc="left", fontweight="bold", fontsize=8.5)
    for ext in ("png", "pdf"):
        fig.savefig(FIG / f"figure3_clinical_profile.{ext}", dpi=300, bbox_inches="tight")
    plt.close(fig)

    c = com.set_index("Gene")
    val = validation(s)
    R = {"comorbidity": {
        "patients_covered_pct": round(100 * covered, 1), "patients_found_pct": round(100 * (s.found == "1").mean(), 1),
        "pct_all": {g: {lab: c.loc[g, f"{lab}, % of all"] for _, lab in FEATURES} for g in GENES},
        "pct_stated": {g: {lab: c.loc[g, f"{lab}, % of stated"] for _, lab in FEATURES} for g in GENES},
        "stated_n": {g: {lab: int(c.loc[g, f"{lab}, stated n"]) for _, lab in FEATURES} for g in GENES},
        "dev_stated_pct_all": round(100 * s.dev_impairment.isin(["Yes", "No"]).mean(), 1),
        "deaths": int((s.deceased == "Yes").sum()), "death_age_median_months": death_age,
        "death_unknown_cause": int(((s.deceased == "Yes") & (~s.cause_of_death.str.split(":").str[0].isin(["SUDEP", "Status epilepticus", "Respiratory/infection", "Other"]))).sum()),
        "kcnt1_death_stated": int(s[(s.gene == "KCNT1")].deceased.isin(["Yes", "No"]).sum()), "death_causes": s[s.deceased == "Yes"].cause_of_death.str.split(":").str[0].value_counts().to_dict(),
        "sudep": int(s.cause_of_death.str.startswith("SUDEP").sum()),
    }, "inheritance": {
        "tested_n": int(tested.sum()), "probands": int(len(prob)), "tested_pct": round(100 * tested.sum() / len(prob), 1),
        "de_novo_pct_per_patient": round(100 * inh_p["De novo"].sum() / tested_p.sum(), 1),
        "de_novo_pct_per_patient_by_gene": (100 * inh_p["De novo"] / tested_p.replace(0, np.nan)).round(1).to_dict(),
        "de_novo_pct": round(100 * inh["De novo"].sum() / tested.sum(), 1),
        "de_novo_pct_by_gene": inh_pct["De novo"].round(1).to_dict(), "tested_by_gene": tested.to_dict(),
        "inherited_by_gene": inherited.to_dict(), "maternal_by_gene": inh["Maternal"].to_dict(), "paternal_by_gene": inh["Paternal"].to_dict(),
        "parental_mosaic_n": int(pm.sum()), "parental_mosaic_by_gene": pm.to_dict(),
        "parental_mosaic_pct_of_inherited": round(100 * pm.sum() / inherited.sum(), 1) if inherited.sum() else None,
        "patient_mosaic_n": int(((s.patient_mosaic == "Yes") & ~parent_lab).sum()),
        "pcdh19_paternal_families": int(inh.loc["PCDH19", "Paternal"]), "pcdh19_maternal_families": int(inh.loc["PCDH19", "Maternal"]),
    }, "comorbidity_validation": val}
    res = json.loads((PROJECT / "analysis" / "results.json").read_text())
    res.update(R)
    (PROJECT / "analysis" / "results.json").write_text(json.dumps(res, indent=1, default=str), encoding="utf-8")
    print(com[["Gene", "Patients, n"] + bcols].to_string(index=False))
    print(out.to_string())
    print(json.dumps({k: v for k, v in R["comorbidity"].items() if k not in ("pct_all", "pct_stated", "stated_n")}, indent=0))


if __name__ == "__main__":
    main()
