"""Korean patients in the literature: identification, within-Korea overlap, comparison with other reports (Figure 5).

Korean papers were verified at the source (korean_papers_verified.csv; all_korean = every counted patient recruited in Korea).
Within-Korea duplicates: candidates (same gene and variant, different Korean papers sharing an institution) were verified
at the source. Primary set: one patient per verified same-patient group (earliest PMID), using patient-level and
cohort-level (nested cohort) evidence; strict set: patient-level evidence only; sensitivity: all reported.
"""
import json

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from scipy import stats  # noqa: E402

from common import FIG, GENES, PROJECT, TAB, load_patients, load_rows, pd, re  # noqa: E402

plt.rcParams.update({"font.family": "Arial", "font.size": 8, "axes.spines.top": False, "axes.spines.right": False,
                     "pdf.fonttype": 42})

INSTITUTIONS = {
    "Severance": ["41573391", "32139178", "29141279", "29455050", "29390993", "29295803", "33000761", "39227335", "32917465",
                  "38073125", "35998430", "31675620", "42314303"],
    "SNUCH/SNUH": ["35087721", "28664031", "39787995", "21868258", "23916143", "25524840", "31572294", "38563110", "31401500",
                   "38946282", "40540813", "33349918"],
    "Samsung": ["32613771", "33067208", "35087721", "28664031", "35190816"],
    "Asan": ["20119593", "35087721", "28664031", "25524840"],
    "Eulji": ["28442529"],
    "Korea Univ.": ["41630268"],
    "Catholic Univ.": ["29263050"],
    "Gachon": ["31675620"],
}
MAIN_SITE = {}
for inst, pm in INSTITUTIONS.items():
    for p in pm:
        MAIN_SITE.setdefault(p, inst)


def vkey(c, p):
    k = re.sub(r"\s", "", str(c)).lower()
    if k.startswith("c."):
        return k
    m = re.search(r"p\.\(?([A-Za-z]{3}\d+[A-Za-z]{3})", str(p))
    return m.group(1).lower() if m else k


def main():
    R = {}
    ver = pd.read_csv(PROJECT / "analysis" / "korean_papers_verified.csv", dtype=str, encoding="utf-8-sig")
    kor = ver[ver.korean_status == "all_korean"]
    kpairs = set(zip(kor.gene, kor.pmid))
    # partly Korean articles: only the listed canonical patients (e.g. one Samsung patient in an international STXBP1 cohort)
    part = ver[ver.korean_status == "partly"]
    kpts = {(g, c.strip()) for g, ids in zip(part.gene, part.korean_patient_ids) for c in str(ids).split(";") if c.strip()}
    rows = load_rows()
    pts = load_patients(rows)
    rows["korean"] = [(g, p) in kpairs or (g, c) in kpts for g, p, c in zip(rows.gene, rows.pmid, rows.canonical_patient_id)]
    pts["korean"] = [(g, p) in kpairs or (g, c) in kpts for g, p, c in zip(pts.gene, pts.pmid, pts.canonical_patient_id)]
    kp = pts[pts.korean].copy()
    kp["vk"] = [vkey(c, p) for c, p in zip(kp.variant_c, kp.protein)]
    kp["sites"] = kp.pmid.map(lambda p: {i for i, pm in INSTITUTIONS.items() if p in pm})

    # within-Korea duplicates, verified at the source (runs/paper/source_qc/korean_clusters_decisions.csv):
    # groups of mention_ids judged to be the same patient, on patient-level evidence or cohort-level nesting.
    # Canonical patient IDs are unique only within a gene, so every key includes the gene.
    kp["uid"] = kp.gene + ":" + kp.canonical_patient_id
    kp["mids"] = kp.supporting_mention_ids.map(lambda s: {x for x in re.split(r"[;,| ]+", str(s)) if x})
    groups = pd.read_csv(PROJECT / "analysis" / "korean_same_patient_groups.csv", dtype=str)

    def collapse(bases):
        cl = dict(zip(kp.uid, kp.uid))
        for _, g in groups[groups.basis.isin(bases)].iterrows():
            mids = set(g.mention_ids.split(";"))
            hit = kp[(kp.gene == g.gene) & kp.mids.map(lambda m: bool(m & mids))]
            if len(hit) > 1:
                root = cl[hit.uid.iloc[0]]
                for u in hit.uid:
                    old = cl[u]
                    for k, v in cl.items():
                        if v == old:
                            cl[k] = root
        c = kp.uid.map(cl)
        return kp.assign(cluster=c, pm=kp.pmid.astype(int)).sort_values("pm").drop_duplicates("cluster", keep="first")

    primary = collapse(["patient-level", "cohort-level"])
    strict = collapse(["patient-level"])
    dup = kp[~kp.uid.isin(primary.uid)]
    kp[["gene", "canonical_patient_id", "pmid", "patient_id", "variant_c", "protein", "phenotype_h", "sex_h", "age_of_onset"]].assign(
        kept_in_primary=kp.uid.isin(primary.uid), kept_in_strict=kp.uid.isin(strict.uid)).to_csv(
        TAB / "supp_korean_patients.csv", index=False)
    R["korean_strict_collapse"] = strict.groupby("gene").size().reindex(GENES).fillna(0).astype(int).to_dict()

    by_gene = pd.DataFrame({"all_reported": kp.groupby("gene").size(), "after_overlap_collapse": primary.groupby("gene").size(),
                            "papers": kp.groupby("gene").pmid.nunique()}).reindex(GENES).fillna(0).astype(int)
    by_gene["global_patients"] = pts.groupby("gene").size().reindex(GENES)
    by_gene["korean_share_pct"] = (100 * by_gene.after_overlap_collapse / by_gene.global_patients).round(1)
    by_gene.to_csv(TAB / "table3_korean_by_gene.csv")
    R["korean_by_gene"] = by_gene.to_dict(orient="index")
    R["korean_totals"] = {"papers": int(kor.pmid.nunique()), "partly_papers": int(part.pmid.nunique()), "all_reported": int(len(kp)), "after_collapse": int(len(primary)),
                          "after_patient_level_only": int(len(strict)), "removed_as_duplicates": int(len(dup)),
                          "same_patient_groups": int(len(groups)), "groups_patient_level": int((groups.basis == "patient-level").sum()),
                          "groups_cohort_level": int((groups.basis == "cohort-level").sum()),
                          "candidate_papers": int(ver.pmid.nunique()), "not_korean_papers": int(ver.loc[ver.korean_status == "not_korean", "pmid"].nunique())}
    R["korean_site_counts"] = primary.pmid.map(MAIN_SITE).value_counts().to_dict()

    # comparisons Korean (primary) vs non-Korean, per gene with >= 20 Korean patients, plus pooled
    comp = []
    other = pts[~pts.korean]
    for g in GENES:
        k = primary[primary.gene == g]
        o = other[other.gene == g]
        if len(k) < 15:
            continue
        rec = {"gene": g, "n_korean": len(k), "n_other": len(o)}
        for label, fk, fo in [("DEE", k.phenotype_h == "DEE", o.phenotype_h == "DEE"),
                              ("Truncating", k.consequence == "Truncating", o.consequence == "Truncating"),
                              ("Female", k.sex_h == "Female", o.sex_h == "Female")]:
            if label == "Female":
                kk, oo = k[k.sex_h != "NS"], o[o.sex_h != "NS"]
                fk, fo = kk.sex_h == "Female", oo.sex_h == "Female"
            else:
                kk, oo = (k[k.phenotype_h != "Not specified"], o[o.phenotype_h != "Not specified"]) if label == "DEE" else (k, o)
                fk = kk.phenotype_h == "DEE" if label == "DEE" else fk
                fo = oo.phenotype_h == "DEE" if label == "DEE" else fo
            a, b, c, d = int(fk.sum()), int(len(fk) - fk.sum()), int(fo.sum()), int(len(fo) - fo.sum())
            rec[f"{label}_korean"] = f"{a}/{a + b} ({100 * a / max(a + b, 1):.1f}%)"
            rec[f"{label}_other"] = f"{c}/{c + d} ({100 * c / max(c + d, 1):.1f}%)"
            rec[f"{label}_p"] = float(f"{stats.fisher_exact([[a, b], [c, d]])[1]:.2g}") if a + b and c + d else None
        ko, oo = k.onset_months.dropna(), o.onset_months.dropna()
        rec["onset_korean"] = f"{ko.median():.1f} ({ko.quantile(.25):.1f}-{ko.quantile(.75):.1f}); n={len(ko)}"
        rec["onset_other"] = f"{oo.median():.1f} ({oo.quantile(.25):.1f}-{oo.quantile(.75):.1f}); n={len(oo)}"
        rec["onset_p"] = float(f"{stats.mannwhitneyu(ko, oo).pvalue:.2g}") if len(ko) > 4 else None
        rec["field_completeness_korean_onset_pct"] = round(100 * len(ko) / len(k), 1)
        rec["field_completeness_other_onset_pct"] = round(100 * len(oo) / len(o), 1)
        comp.append(rec)
    comp = pd.DataFrame(comp)
    comp.to_csv(TAB / "supp_korean_vs_other.csv", index=False)
    R["korean_vs_other"] = comp.to_dict(orient="records")

    # SCN1A recurrent variants in Korea vs elsewhere
    k1 = primary[primary.gene == "SCN1A"].protein.value_counts().head(8)
    R["korean_scn1a_top_variants"] = k1.to_dict()

    # sensitivity: all Korean-reported patients
    sens = []
    for g in ["SCN1A"]:
        k = kp[kp.gene == g]
        o = other[other.gene == g]
        kk, oo = k[k.phenotype_h != "Not specified"], o[o.phenotype_h != "Not specified"]
        sens.append({"gene": g, "n": len(k), "DEE_pct": round(100 * (kk.phenotype_h == "DEE").mean(), 1),
                     "DEE_other_pct": round(100 * (oo.phenotype_h == "DEE").mean(), 1)})
    R["korean_sensitivity_all_reported"] = sens

    # ------------------------------------------------ figure
    fig, axes = plt.subplots(1, 3, figsize=(7.2, 2.7), gridspec_kw={"width_ratios": [1.1, 1.2, 1.1]})
    ax = axes[0]
    sites = ["Severance", "SNUCH/SNUH", "Samsung", "Asan", "Eulji", "Korea Univ.", "Catholic Univ.", "Gachon"]
    scol = ["#1F78B4", "#33A02C", "#6A3D9A", "#FF7F00", "#B15928", "#A6CEE3", "#FB9A99", "#CAB2D6"]
    primary = primary.assign(site=primary.pmid.map(MAIN_SITE))
    ct = pd.crosstab(primary.gene, primary.site).reindex(index=GENES, columns=sites, fill_value=0)
    left = np.zeros(len(GENES))
    for s, c in zip(sites, scol):
        ax.barh(range(len(GENES)), ct[s], left=left, color=c, label=s, height=0.7)
        left += ct[s].values
    for i, g in enumerate(GENES):
        ax.text(left[i] + 3, i, f"{int(left[i])} (all: {int(by_gene.loc[g, 'all_reported'])})", va="center", fontsize=6)
    ax.set_yticks(range(len(GENES)), GENES)
    ax.invert_yaxis()
    ax.set_xlim(0, max(left) * 1.45)
    ax.set_xlabel("Korean patients\n(after overlap collapse)")
    ax.legend(fontsize=5.5, frameon=False, loc="lower right")
    ax.set_title("a  Korean reports by institution", loc="left", fontsize=8, fontweight="bold")

    ax = axes[1]
    cats = ["DEE", "GEFS+ / febrile spectrum", "Other epilepsy", "No seizures"]
    k = primary[(primary.gene == "SCN1A") & (primary.phenotype_h != "Not specified")]
    o = other[(other.gene == "SCN1A") & (other.phenotype_h != "Not specified")]
    x = np.arange(len(cats))
    ax.bar(x - 0.2, [100 * (k.phenotype_h == c).mean() for c in cats], 0.4, color="#C51B7D", label=f"Korea (n={len(k)})")
    ax.bar(x + 0.2, [100 * (o.phenotype_h == c).mean() for c in cats], 0.4, color="#7F7F7F", label=f"Other (n={len(o)})")
    ax.set_xticks(x, ["DEE (Dravet)", "GEFS+", "Other epilepsy", "No seizures"], rotation=30, ha="right")
    ax.set_ylabel("% of SCN1A patients")
    ax.legend(fontsize=6, frameon=False)
    ax.set_title("b  SCN1A phenotype", loc="left", fontsize=8, fontweight="bold")

    ax = axes[2]
    gs = [g for g in GENES if (primary.gene == g).sum() >= 15]
    for i, g in enumerate(gs):
        for j, (df, col) in enumerate([(primary[primary.gene == g], "#C51B7D"), (other[other.gene == g], "#7F7F7F")]):
            v = df.onset_months.dropna()
            if len(v) < 3:
                continue
            q = np.log10(v.clip(lower=0.03).quantile([.25, .5, .75]))
            ax.plot([i + (j - 0.5) * 0.3] * 2, [q.iloc[0], q.iloc[2]], color=col, lw=1.5)
            ax.plot(i + (j - 0.5) * 0.3, q.iloc[1], "o", color=col, ms=4, label=("Korea" if j == 0 else "Other") if i == 0 else None)
    ticks = [0.1, 1, 3, 12, 60]
    ax.set_yticks(np.log10(ticks), ["3 d", "1 m", "3 m", "1 y", "5 y"])
    ax.set_xticks(range(len(gs)), gs, rotation=45, ha="right", fontsize=6.5)
    ax.set_ylabel("Onset age, median (IQR)")
    ax.legend(fontsize=6, frameon=False)
    ax.set_title("c  Age at onset", loc="left", fontsize=8, fontweight="bold")
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(FIG / f"figureS1_korean.{ext}", dpi=300, bbox_inches="tight")
    res = json.loads((PROJECT / "analysis" / "results.json").read_text())
    res.update(R)
    (PROJECT / "analysis" / "results.json").write_text(json.dumps(res, indent=1, default=str), encoding="utf-8")
    print(json.dumps(R, indent=1, default=str))


if __name__ == "__main__":
    main()
