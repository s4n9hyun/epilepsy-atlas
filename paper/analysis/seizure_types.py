"""Reported seizure types by gene (Figure 3).

Free-text seizure_types of each unique patient is mapped to ILAE-based categories by keyword rules (a patient can have several).
Denominator: patients with at least one classifiable seizure type. A type not mentioned is not evidence of absence, so the
percentages are reported frequencies, not prevalences.
Validated on a blinded source check of random patients (runs/paper/seizure_type_validation).
Writes figures/figure3_seizure_types.*, tables/supp_seizure_types_by_gene.csv, tables/supp_seizure_types_validation.csv and merges results into results.json.
"""
import json

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from scipy import stats  # noqa: E402

from common import FIG, GENES, PROJECT, TAB, WORKSPACE, is_missing, load_patients, load_rows, pd, re  # noqa: E402

plt.rcParams.update({"font.family": "Arial", "font.size": 8, "axes.spines.top": False, "axes.spines.right": False,
                     "pdf.fonttype": 42, "svg.fonttype": "none"})

# order = column order in the figure; each rule is applied to the lower-cased text
TYPES = [
    # "tonic/clonic/myoclonic spasm" describes neonatal seizure semiology, not epileptic spasms (reviewer comment R1-6)
    ("Epileptic spasms", r"epileptic spasm|infantile spasm|west syndrome|salaam|hypsarrhythmi|\biess\b|"
                         r"(?<!tonic )(?<!tonic-)(?<!clonic )(?<!clonic-)(?<!myoclonic )(?<!myoclonic-)\bspasms?\b"),
    ("Tonic", r"(?<![-a-z])tonic\b(?![\s-]*clonic)"),
    ("Focal", r"focal|partial|hemiclon|hemiconvuls|unilateral|migrating|dyscognitive|hypermotor|automatism|versive"),
    ("Hemiclonic", r"hemiclon|hemi-clon|hemiconvuls|unilateral (?:clonic|tonic-clonic|seizure|convuls)|alternating"),
    ("Bilateral tonic-clonic", r"tonic[\s-]*clonic|\bgtcs?\b|grand mal|secondar\w* generali|convuls"),
    ("Clonic", r"(?<!tonic-)(?<!tonic )(?<!hemi)(?<!hemi-)(?<!focal )(?<!unilateral )\bclonic\b"),
    ("Myoclonic", r"myoclon"),
    ("Absence", r"absence|petit mal"),
    ("Atonic", r"atonic|drop attack|astatic"),
    ("Febrile seizures", r"(?<!a)febrile|fever|hyperthermi"),
    ("Status epilepticus", r"status epilepticus|\bc?se\b|prolonged (?:convulsive |febrile )?seizure|\bncse\b"),
]
LABEL = {"Bilateral tonic-clonic": "Tonic-clonic"}
# febrile and status epilepticus are often recorded with the clinical course rather than in the seizure-type field
COURSE_TYPES = {"Febrile seizures", "Status epilepticus"}
NEG = r"\bno\b|\bnot\b|without|absent|denied|family|mother|father|sibling|sister|brother|relative|parent"


def clauses(text):
    t = str(text).lower().replace("–", "-").replace("—", "-")
    return [c for c in re.split(r"[;.]", t) if c.strip()]


def classify(text, course=""):
    t = str(text).lower()
    if is_missing(t) or re.fullmatch(r"(not specified.*|seizures?|epilepsy|unknown.*)", t.strip()):
        return set()
    out = set()
    for name, rx in TYPES:
        for c in clauses(text):
            if re.search(rx, c) and not (name in COURSE_TYPES and re.search(NEG, c)):
                out.add(name)
        if name in COURSE_TYPES and not is_missing(course):
            if any(re.search(rx, c) and not re.search(NEG, c) for c in clauses(course)):
                out.add(name)
    return out


VAL = WORKSPACE / "runs/paper/seizure_type_validation"
VCOL = {"Epileptic spasms": "spasms", "Tonic": "tonic", "Focal": "focal", "Hemiclonic": "hemiclonic",
        "Bilateral tonic-clonic": "tonic_clonic", "Clonic": "clonic", "Myoclonic": "myoclonic", "Absence": "absence",
        "Atonic": "atonic", "Febrile seizures": "febrile", "Status epilepticus": "status"}


def validation(pts, sub="test", out="supp_seizure_types_validation.csv"):
    """Blinded source check of random patients (20 per gene): categories read from the source vs categories derived from
    the dataset by the rules above. 'dev' = sample used once to refine the rules; 'test' = fresh sample drawn afterwards."""
    d = VAL / sub
    key = pd.read_csv(d / "key.csv", dtype=str, keep_default_na=False)[["gene", "canonical_patient_id"]]
    res = pd.concat([pd.read_csv(f, dtype=str, keep_default_na=False) for f in sorted(d.glob("results/*.csv"))])
    m = key.merge(res, on=["gene", "canonical_patient_id"])
    cur = pts.set_index(["gene", "canonical_patient_id"])
    m = m[[(g, c) in cur.index for g, c in zip(m.gene, m.canonical_patient_id)]]  # patients excluded from the dataset since sampling
    m["our_categories"] = ["; ".join(sorted(classify(cur.loc[(g, c)].seizure_types, cur.loc[(g, c)].other_phenotypes)))
                           for g, c in zip(m.gene, m.canonical_patient_id)]
    ok = m[(m.found == "1") & (m.types_reported == "1")]
    rows, tp_all, fp_all, fn_all = [], 0, 0, 0
    for name, col in VCOL.items():
        ours = ok.our_categories.str.split("; ").map(lambda x, n=name: n in x)
        src = ok[col].str.strip() == "1"
        tp, fp, fn = int((ours & src).sum()), int((ours & ~src).sum()), int((~ours & src).sum())
        tp_all, fp_all, fn_all = tp_all + tp, fp_all + fp, fn_all + fn
        rows.append({"Seizure type": LABEL.get(name, name), "Reported at source, n": tp + fn, "Assigned in dataset, n": tp + fp,
                     "Sensitivity, %": round(100 * tp / (tp + fn), 1) if tp + fn else None,
                     "Precision, %": round(100 * tp / (tp + fp), 1) if tp + fp else None})
    rows.append({"Seizure type": "All types", "Reported at source, n": tp_all + fn_all, "Assigned in dataset, n": tp_all + fp_all,
                 "Sensitivity, %": round(100 * tp_all / (tp_all + fn_all), 1), "Precision, %": round(100 * tp_all / (tp_all + fp_all), 1)})
    t = pd.DataFrame(rows)
    t.to_csv(TAB / out, index=False)
    return {"sampled": len(m), "found_with_types": len(ok), "sensitivity": float(t.iloc[-1]["Sensitivity, %"]),
            "precision": float(t.iloc[-1]["Precision, %"]),
            "sens_min": float(t.iloc[:-1]["Sensitivity, %"].min()), "prec_min": float(t.iloc[:-1]["Precision, %"].min())}


def main():
    pts = load_patients()
    pts["types"] = [classify(t, o) for t, o in zip(pts.seizure_types, pts.other_phenotypes)]
    # patients with seizures: seizure presence recorded as yes, or a seizure type described (reviewer comment R4-M5b)
    pts = pts[pts.seizure_presence.str.lower().str.startswith("yes") | (pts.types.map(len) > 0)].copy()
    has = pts[pts.types.map(len) > 0]
    names = [n for n, _ in TYPES]
    tab = pd.DataFrame({n: has.groupby("gene").types.apply(lambda s, n=n: 100 * s.map(lambda x: n in x).mean()) for n in names})
    tab = tab.reindex(GENES)
    npat = has.groupby("gene").size().reindex(GENES)
    out = tab.round(1)
    out.insert(0, "Patients with a classified seizure type, n", npat)
    out.insert(1, "Patients with seizures, n", pts.groupby("gene").size().reindex(GENES))
    out.to_csv(TAB / "supp_seizure_types_by_gene.csv")
    tab.to_csv(TAB / "_raw_seizure_types_pct.csv")  # unrounded, for Figure 3A

    # figure: heatmap
    fig, ax = plt.subplots(figsize=(7.0, 4.2))
    im = ax.imshow(tab.values, cmap="Reds", vmin=0, vmax=80, aspect="auto")
    ax.set_xticks(range(len(names)))
    ax.set_xticklabels([LABEL.get(n, n) for n in names], rotation=40, ha="right")
    ax.set_yticks(range(len(GENES)))
    ax.set_yticklabels([f"{g}  (n={npat[g]})" for g in GENES])
    for i in range(len(GENES)):
        for j in range(len(names)):
            v = tab.values[i, j]
            ax.text(j, i, f"{v:.0f}", ha="center", va="center", fontsize=6.5, color="white" if v >= 50 else "black")
    for s in ax.spines.values():
        s.set_visible(False)
    ax.tick_params(length=0)
    cb = fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02)
    cb.set_label("Patients with the seizure type reported, %")
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(FIG / f"panel_seizure_types.{ext}", dpi=300, bbox_inches="tight")
    plt.close(fig)

    # selected contrasts (chi-square / Fisher: gene vs all other genes)
    def contrast(gene, typ):
        a = int(has[has.gene == gene].types.map(lambda x: typ in x).sum())
        b = int((has.gene == gene).sum()) - a
        c = int(has[has.gene != gene].types.map(lambda x: typ in x).sum())
        d = int((has.gene != gene).sum()) - c
        return {"gene_pct": round(100 * a / (a + b), 1), "others_pct": round(100 * c / (c + d), 1),
                "p": float(f"{stats.fisher_exact([[a, b], [c, d]])[1]:.2g}")}

    fx = has[has.gene == "FOXG1"]
    gains = set(load_rows()[lambda d: (d.gene == "FOXG1") & (d.cnv_type == "gain")].canonical_patient_id)
    fxn = fx[~fx.canonical_patient_id.isin(gains)]
    foxg1_sp = {"all": round(100 * fx.types.map(lambda x: "Epileptic spasms" in x).mean(), 1),
                "without_gains": round(100 * fxn.types.map(lambda x: "Epileptic spasms" in x).mean(), 1),
                "gain_carriers_with_spasms": f"{int(fx[fx.canonical_patient_id.isin(gains)].types.map(lambda x: 'Epileptic spasms' in x).sum())}/{int(fx.canonical_patient_id.isin(gains).sum())}"}
    kq = has[has.gene == "KCNQ2"]
    kq_ph = {ph: round(100 * kq[kq.phenotype_h == ph].types.map(lambda x: "Tonic" in x).mean(), 1)
             for ph in ["DEE", "Self-limited (familial) neonatal/infantile epilepsy"] if (kq.phenotype_h == ph).any()}
    top = {g: [n for n in tab.loc[g].sort_values(ascending=False).index[:3]] for g in GENES}
    R = {"seizure_types": {
        "patients_with_seizures": len(pts), "patients_classified": len(has),
        "classified_pct": round(100 * len(has) / len(pts), 1), "foxg1_spasms": foxg1_sp, "kcnq2_tonic_by_phenotype": kq_ph,
        "pct": tab.round(3).to_dict(orient="index"), "n": npat.to_dict(), "top3": top,
        "contrasts": {f"{g}|{t}": contrast(g, t) for g, t in [
            ("SCN1A", "Febrile seizures"), ("SCN1A", "Myoclonic"), ("SCN1A", "Hemiclonic"), ("SCN1A", "Status epilepticus"),
            ("CDKL5", "Epileptic spasms"), ("STXBP1", "Epileptic spasms"), ("FOXG1", "Epileptic spasms"),
            ("KCNQ2", "Tonic"), ("SCN8A", "Tonic"), ("KCNT1", "Focal"), ("PCDH19", "Focal"), ("PCDH19", "Febrile seizures"),
            ("SYNGAP1", "Absence"), ("SYNGAP1", "Myoclonic"), ("SYNGAP1", "Atonic"), ("GABRA1", "Absence"), ("GABRA1", "Myoclonic")]},
    }}
    allp = load_patients()
    R["seizure_types"]["validation_dev"] = validation(allp, "dev", "supp_seizure_types_validation_dev.csv")
    if list((VAL / "test" / "results").glob("*.csv")):
        R["seizure_types"]["validation"] = validation(allp, "test")
    res = json.loads((PROJECT / "analysis" / "results.json").read_text())
    res.update(R)
    (PROJECT / "analysis" / "results.json").write_text(json.dumps(res, indent=1, default=str), encoding="utf-8")
    print(out.to_string())
    print(json.dumps(R["seizure_types"]["contrasts"], indent=0))


if __name__ == "__main__":
    main()
