"""Fill manuscript placeholders from results.json and tables, check Epilepsia limits, and write the build copies.

manuscript/tex/main.tex (manuscript + Supporting Information) -> manuscript/build/main.tex ({{...}} filled, Table 1 inserted)
build.sh then compiles paper/epilepsy_main.pdf (tectonic) and the two Word files for Epilepsia submission (pandoc).
"""
import json
import re
import subprocess
from pathlib import Path

import pandas as pd

PROJECT = Path(__file__).resolve().parents[1]
TEX = PROJECT / "manuscript" / "tex"
TAB = PROJECT / "tables"
OUT = PROJECT / "manuscript" / "build"
CHANGE_LOG = Path(__import__("os").environ.get("EPILEPSY_WORKSPACE", Path(__file__).resolve().parents[2])) / "outputs/current/curation_change_log.csv"
R = json.loads((PROJECT / "analysis" / "results.json").read_text())
GENES = ["SCN1A", "SCN2A", "SCN8A", "KCNQ2", "KCNT1", "CDKL5", "PCDH19", "STXBP1", "SYNGAP1", "GABRA1", "FOXG1", "SMC1A", "EEF1A2"]
X_LINKED = ["CDKL5", "PCDH19", "SMC1A"]
DBC_GENES = ["SCN1A", "SCN2A", "KCNQ2", "STXBP1"]


def pnum(p):
    p = float(p)
    return "P < 0.001" if p < 0.001 else f"P = {p:.3f}" if p < 0.01 else f"P = {p:.2f}"


def pct_of(s):
    return re.search(r"\(([\d.]+)%\)", s).group(1)


def esc(v):
    s = str(v)
    for a, b in [("\\", r"\textbackslash{}"), ("&", r"\&"), ("%", r"\%"), ("_", r"\_"), ("#", r"\#"), ("$", r"\$")]:
        s = s.replace(a, b)
    return s


def latex_table(df, title, note="", size=r"\small"):
    cols = list(df.columns)
    spec = "l" * len(cols)
    lines = [r"\begingroup" + size, r"\setlength{\tabcolsep}{4pt}", rf"\begin{{longtable}}{{{spec}}}",
             rf"\caption*{{{title}}}\\", r"\toprule", " & ".join(rf"\textbf{{{esc(c)}}}" for c in cols) + r" \\", r"\midrule", r"\endhead"]
    lines += [" & ".join(esc(v) for v in row) + r" \\" for row in df.itertuples(index=False)]
    lines += [r"\bottomrule", r"\end{longtable}"]
    if note:
        lines.append(rf"\par\footnotesize {esc(note)}")
    lines.append(r"\endgroup")
    return "\n".join(lines)


def table1():
    """Genes as rows (reviewer comment R4 minor: legibility)."""
    t = pd.read_csv(TAB / "table1_gene_summary.csv", dtype=str, keep_default_na=False)
    pc = lambda x: (re.search(r"\(([\d.]+)\)", x).group(1) if re.search(r"\(([\d.]+)\)", x) else "")
    t = pd.DataFrame({"Gene": t.Gene, "Articles": t.Articles, "Patients": t.Patients,
                      "Female, %": t["Female, n/N (%)"].map(pc), "Onset, median (IQR), mo": t["Onset, median (IQR), mo; n"].str.split(";").str[0],
                      "DEE, %": t["DEE label, n/N (%)"].map(pc), "Missense, %": t["Missense, n (%)"].map(pc),
                      "Truncating, %": t["Truncating, n (%)"].map(pc), "CNV, %": t["CNV/exon-level, n (%)"].map(pc),
                      "ClinVar P/LP, %": t["ClinVar P/LP, n (%)"].map(pc), "ClinVar VUS/conflicting, %": t["ClinVar VUS/conflicting, n (%)"].map(pc),
                      "ClinVar B/LB, %": t["ClinVar B/LB, n (%)"].map(pc)})
    note = ("DEE: reported DEE label among patients with a phenotype category; female among patients with known sex; onset age parsed "
            "from free text (numerators, denominators and the number with a parsed onset: [[T:table1counts]]). Consequence and ClinVar columns are per patient–variant row; ClinVar from the allele-level "
            "match to the release of " + R["variant_interpretation"]["clinvar_release"] + " (small variants only; rows without a match or with "
            "copy-number variants are in neither column). "
            "Abbreviations: B/LB, benign/likely benign; CNV, copy-number variant; DEE, developmental and epileptic encephalopathy; "
            "IQR, interquartile range; P/LP, pathogenic/likely pathogenic; VUS, variant of uncertain significance.")
    cols = list(t.columns)
    rows = [" & ".join((rf"\textit{{{esc(v)}}}" if i == 0 and v in GENES else esc(v)) for i, v in enumerate(r)) + r" \\" for r in t.itertuples(index=False)]
    rows.insert(len(rows) - 1, r"\midrule")
    head = " & ".join(rf"\textbf{{{esc(c)}}}" for c in cols) + r" \\"
    return "\n".join([r"\begin{table*}[t]", r"\caption{Dataset summary by gene.}\label{tab1}", r"\centering\footnotesize",
                      r"\resizebox{\textwidth}{!}{\begin{tabular}{l" + "r" * (len(cols) - 1) + "}", r"\toprule", head, r"\midrule", *rows,
                      r"\bottomrule", r"\end{tabular}}", r"\begin{tablenotes}", rf"\item {esc(note)}", r"\end{tablenotes}", r"\end{table*}"])


def values():
    scb = R["scb"]["SCN1A"]
    ver = R["scb_source_verification"]
    dee = R["scb_scn1a_dee_vs_nondee_worsening"]
    by = R["scb_scn1a_by_phenotype"]
    dee_n, dee_w = by["DEE"]["n"], by["DEE"]["worsening"]
    oth_n = sum(v["n"] for k, v in by.items() if k != "DEE")
    oth_w = sum(v["worsening"] for k, v in by.items() if k != "DEE")
    cl = R["clustered_sensitivity"]["SCN1A: SCB worsening, DEE vs non-DEE"]
    sig = dee["p"] < 0.05
    sentence = (f"Worsening was {'more' if sig else 'not significantly more'} frequent in DEE (mostly Dravet syndrome) "
                f"than in other SCN1A phenotypes ({dee_w}/{dee_n} vs {oth_w}/{oth_n}; OR {dee['OR']:.2f}, "
                f"{dee['CI'][0]:.2f}--{dee['CI'][1]:.2f}; {pnum(dee['p'])}"
                + (f"; clustered OR {cl['GEE_OR']:.2f}, {cl['GEE_CI'][0]:.2f}--{cl['GEE_CI'][1]:.2f}" if sig else "") + ").")
    return {
        "S1_N": str(scb["n"]),
        "S1_W_PCT": pct_of(scb["Worsening"]),
        "S1_B_PCT": pct_of(scb["Benefit"]),
        "S1_NB_PCT": pct_of(scb["No benefit"]),
        "S1_DEE_SENTENCE": sentence,
        "VER_PAIRS": str(ver["pairs"]),
        "VER_AGREE": str(ver["auto_agreement_pct"]),
        "VER_WRONG": str(ver["db_text_wrong"]),
        "VER_WRONG_PCT": f"{100 * ver['db_text_wrong'] / ver['pairs']:.1f}",
        **recurrent_values(),
        "VER_LABELS": esc(", ".join(f"{k} {v}" for k, v in ver["labels"].items())),
    }


def recurrent_values():
    rv = pd.DataFrame(R["recurrent_variants"])
    rv["gap"] = rv.naive - rv.dedup
    top = rv[rv.protein.str.match(r"p\.[A-Z]")].sort_values(["gap", "naive"], ascending=False).head(6)
    parts = [f"{g} {p} ({n} vs {d}" + (" unique patients)" if i == 0 else ")")
             for i, (g, p, n, d) in enumerate(zip(top.gene, top.protein, top.naive, top.dedup))]
    sentence = ("Counting every description as a patient (naive counting) overstated the number of patients carrying " + ", ".join(parts[:-1])
                + " and " + parts[-1] + ".")
    mt = R["missense_top"]
    fmt = lambda g, k, first=False: f"{g} " + and_join([f"{p} ({n}" + (" families)" if first and i == 0 else ")")
                                                     for i, (p, n) in enumerate(list(mt[g].items())[:k])])
    return {"RECURRENT_SENTENCE": sentence, "MAX_FOLD": f"{rv.inflation.max():.1f}",
            "MIS_SCN8A": fmt("SCN8A", 2, True), "MIS_KCNT1": fmt("KCNT1", 2), "MIS_KCNQ2": fmt("KCNQ2", 2)}


def table2():
    t = pd.read_csv(TAB / "table2_variant_interpretation.csv", dtype=str, keep_default_na=False)
    t = t.rename(columns={"Distinct variants, n": "Variants, n", "Pathogenic/likely pathogenic": "P/LP", "Uncertain significance": "VUS",
                          "Benign/likely benign": "B/LB", "P/LP with ≥2-star review": "P/LP, ≥2 stars"})
    note = (f"Distinct small variants (GRCh38 alleles; MANE Select and exons absent from MANE; copy-number and exon-level variants excluded), "
            f"matched at the allele level to the ClinVar release of {R['variant_interpretation']['clinvar_release']}. Patients analysed: "
            "patients with such a variant (% of all patients of the gene). Percentages of variants in the gene. Other: ClinVar records without "
            "a germline classification. ≥2 stars: multiple submitters without conflicts, expert panel or practice guideline. "
            "Abbreviations: B/LB, benign/likely benign; P/LP, pathogenic/likely pathogenic; VUS, uncertain significance.")
    cols = list(t.columns)
    rows = [" & ".join((rf"\textit{{{esc(v)}}}" if i == 0 and v in GENES else esc(v)) for i, v in enumerate(r)) + r" \\" for r in t.itertuples(index=False)]
    rows.insert(len(rows) - 1, r"\midrule")
    return "\n".join([r"\begin{table*}[t]", r"\caption{ClinVar status of the variants reported in the atlas.}\label{tab2}", r"\centering\footnotesize",
                      r"\resizebox{\textwidth}{!}{\begin{tabular}{l" + "r" * (len(cols) - 1) + "}", r"\toprule",
                      " & ".join(rf"\textbf{{{esc(c)}}}" for c in cols) + r" \\", r"\midrule", *rows, r"\bottomrule",
                      r"\end{tabular}}", r"\begin{tablenotes}", rf"\item {esc(note)}", r"\end{tablenotes}", r"\end{table*}"])


def new_values():
    """Seizure types (Figure 3) and ClinVar status (Table 2)."""
    st, vi = R["seizure_types"], R["variant_interpretation"]
    pc, con = st["pct"], st["contrasts"]
    val = st.get("validation") or st["validation_dev"]
    sp = lambda g, t: f"{pc[g][t]:.0f}"
    v = dict(ST_N_SZ=n(st["patients_with_seizures"]), ST_N_CL=n(st["patients_classified"]), ST_PCT=f1(st["classified_pct"]),
             ST_VAL_N=n(val["found_with_types"]), ST_VAL_SENS=f1(val["sensitivity"]), ST_VAL_PREC=f1(val["precision"]))
    for g, t, k in [("CDKL5", "Epileptic spasms", "CDKL5_SP"), ("STXBP1", "Epileptic spasms", "STXBP1_SP"), ("FOXG1", "Epileptic spasms", "FOXG1_SP"),
                    ("SCN2A", "Epileptic spasms", "SCN2A_SP"), ("SCN1A", "Epileptic spasms", "SCN1A_SP"), ("PCDH19", "Epileptic spasms", "PCDH19_SP"),
                    ("KCNQ2", "Tonic", "KCNQ2_TO"), ("CDKL5", "Tonic", "CDKL5_TO"), ("KCNT1", "Focal", "KCNT1_FO"), ("SMC1A", "Focal", "SMC1A_FO"),
                    ("PCDH19", "Focal", "PCDH19_FO"), ("PCDH19", "Febrile seizures", "PCDH19_FE"), ("GABRA1", "Febrile seizures", "GABRA1_FE"),
                    ("EEF1A2", "Myoclonic", "EEF1A2_MY"), ("EEF1A2", "Absence", "EEF1A2_AB")]:
        v["ST_" + k] = sp(g, t)
    for key, k in [("SCN1A|Febrile seizures", "SCN1A_FE"), ("SCN1A|Status epilepticus", "SCN1A_SE"), ("SCN1A|Hemiclonic", "SCN1A_HC"),
                   ("SCN1A|Myoclonic", "SCN1A_MY"), ("SYNGAP1|Myoclonic", "SYNGAP1_MY"), ("SYNGAP1|Absence", "SYNGAP1_AB"),
                   ("SYNGAP1|Atonic", "SYNGAP1_AT")]:
        v["ST_" + k] = f"{con[key]['gene_pct']:.0f}"
        v["ST_" + k + "_O"] = f"{con[key]['others_pct']:.0f}"
    claim(max(pc, key=lambda g: pc[g]["Epileptic spasms"]) == "CDKL5", "spasms highest in CDKL5")
    claim(max(pc, key=lambda g: pc[g]["Myoclonic"]) == "SYNGAP1" and max(pc, key=lambda g: pc[g]["Absence"]) == "SYNGAP1",
          "SYNGAP1 highest myoclonic and absence")
    claim(max(pc, key=lambda g: pc[g]["Focal"]) == "KCNT1", "focal highest in KCNT1")
    claim(max(pc, key=lambda g: pc[g]["Hemiclonic"]) == "SCN1A", "hemiclonic highest in SCN1A")
    sp_pct = vi["status_pct"]
    ab = vi["absent_by_gene_pct"]
    lo, hi = min(ab, key=ab.get), max(ab, key=ab.get)
    frac_word = lambda x: ("half" if 45 <= x < 55 else "two in five" if 35 <= x < 45 else "one third" if 28 <= x < 35
                           else "one quarter" if 20 <= x < 28 else f"{x:.0f}%")
    v.update(VI_RELEASE=vi["clinvar_release"], VI_N=n(vi["unique_variants"]), VI_PAT=n(vi["patients"]), VI_COV=f1(vi["coverage_pct"]),
             VI_SINGLE_PCT=f1(vi["single_family_pct"]), VI_PLP_PCT=f1(sp_pct["Pathogenic/likely pathogenic"]), VI_PLP2=n(vi["plp_2star"]),
             VI_VC_PCT=f1(sp_pct["Uncertain significance"] + sp_pct["Conflicting"]), VI_ABS_PCT=f1(vi["absent_pct"]),
             VI_ABS_N=n(vi["absent_n"]), VI_ABS_LO=f"{ab[lo]:.0f}", VI_ABS_LO_G=lo, VI_ABS_HI=f"{ab[hi]:.0f}", VI_ABS_HI_G=hi,
             VI_ABS_WORD=frac_word(vi["absent_pct"]),
             VI_PAT_NOTPLP_PCT=f1(vi["patients_with_not_plp_variant_pct"]), VI_POOL=n(vi["pool_n"]), VI_AF_UNK=n(vi.get("pool_af_unknown", 0)), VI_TOO_COMMON=n(vi["pool_too_common"]),
             VI_CAND_N=n(vi["candidates_n"]), VI_CAND_FAM=n(vi["candidates_families"]), VI_CAND_ABS=n(vi["candidates_absent_n"]),
             VI_CAND_VC=n(vi["candidates_vus_conf_n"]), VI_CAND_DN=n(vi["candidates_with_de_novo"]).capitalize() if False else n(vi["candidates_with_de_novo"]),
             VI_CAND_DN2=n(vi["candidates_2plus_de_novo"]),
             VI_BLB_N=n(vi["blb_n"]), VI_BLB_PAT=n(vi["blb_patients"]),
             CNV_UNK=n(vi["cnv_extent_unknown"]), CNV_N=n(vi["cnv_rows"]), CNV_EXCL=n(vi["cnv_excluded_multigene"]),
             VI_A_PLP_CV=f1(vi["author_plp_agree_pct"] or 0), VI_A_VUS_PLP=n(vi["author_vus_now_plp"]), VI_A_VUS_N=n(vi["author_vus_with_clinvar_class"]))
    return v


def cm_values():
    """Comorbidities and inheritance (Figure 3B-C)."""
    cm, inh, val = R["comorbidity"], R["inheritance"], R.get("comorbidity_validation") or {}
    pa = cm["pct_all"]
    g = lambda gene, lab: f"{pa[gene][lab]:.0f}"
    dev = {k: v["Developmental impairment"] for k, v in pa.items()}
    dn = inh["de_novo_pct_by_gene"]
    v = dict(CM_DEV_STATED=f1(cm["dev_stated_pct_all"]), CM_DEV_LO=f"{dev['SCN1A']:.0f}", CM_DEV_HI=f"{dev['SYNGAP1']:.0f}",
             CM_SEV_STXBP1=g("STXBP1", "Severe/profound impairment"), CM_SEV_CDKL5=g("CDKL5", "Severe/profound impairment"),
             CM_AUT_SYNGAP1=g("SYNGAP1", "Autism"), CM_AUT_EEF1A2=g("EEF1A2", "Autism"), CM_AUT_PCDH19=g("PCDH19", "Autism"),
             CM_MOV_FOXG1=g("FOXG1", "Movement disorder"), CM_MOV_STXBP1=g("STXBP1", "Movement disorder"),
             CM_MIC_FOXG1=g("FOXG1", "Microcephaly"), CM_MIC_SMC1A=g("SMC1A", "Microcephaly"),
             CM_CVI_CDKL5=g("CDKL5", "Cortical visual impairment"), CM_DEC_KCNT1=g("KCNT1", "Deceased"),
             CM_DEATHS=n(cm["deaths"]), CM_SUDEP=n(cm["sudep"]), CM_SE=n(cm["death_causes"].get("Status epilepticus", 0)),
             IN_TESTED=n(inh["tested_n"]), IN_TESTED_PCT=f1(inh["tested_pct"]), IN_DN_PCT=f1(inh["de_novo_pct"]),
             IN_PAT_PCDH19=n(inh["paternal_by_gene"]["PCDH19"]), IN_MAT_PCDH19=n(inh["maternal_by_gene"]["PCDH19"]),
             IN_PM_N=n(inh["parental_mosaic_n"]), IN_PM_SCN1A=n(inh["parental_mosaic_by_gene"]["SCN1A"]), IN_PTM_N=n(inh["patient_mosaic_n"]),
             CM_VAL_N=n(val.get("patients", 0)), CM_VAL_AGREE=f1(val.get("agreement_pct", 0)),
             CM_VAL_PMIN=f1(val.get("positive_agreement_min", 0)), CM_VAL_PMAX=f"{val.get('positive_agreement_max', 0):.0f}")
    for gene in ["CDKL5", "STXBP1", "SYNGAP1", "KCNQ2", "PCDH19", "SCN1A"]:
        v["IN_DN_" + gene] = f"{dn[gene]:.0f}"
    ps_, sn = cm["pct_stated"], cm["stated_n"]
    gs = lambda gene, lab: f"{ps_[gene][lab]:.0f}" if ps_[gene][lab] is not None else "--"
    deaths = pd.read_csv(TAB / "supp_deaths_by_gene.csv", dtype=str).set_index("Gene")
    v.update(CM_DEV_LO_S=gs("SCN1A", "Developmental impairment"), CM_DEV_HI_S=gs("SYNGAP1", "Developmental impairment"),
             CM_AUT_SYNGAP1_S=gs("SYNGAP1", "Autism"), CM_MOV_FOXG1_S=gs("FOXG1", "Movement disorder"),
             CM_CVI_CDKL5_S=gs("CDKL5", "Cortical visual impairment"), CM_CVI_CDKL5_N=n(sn["CDKL5"]["Cortical visual impairment"]),
             CM_DEC_KCNT1_S=gs("KCNT1", "Deceased"), CM_DEC_KCNT1_N=n(sn["KCNT1"]["Deceased"]),
             CM_DUNK=n(cm["death_unknown_cause"]), CM_DAGE=deaths.loc["All genes", "Age at death, median (IQR), months"].replace("–", "--"),
             IN_PROBANDS=n(inh["probands"]))
    claim(min(dev, key=dev.get) == "SCN1A" and max(dev, key=dev.get) == "SYNGAP1", "developmental impairment lowest SCN1A, highest SYNGAP1")
    top = lambda lab: sorted(pa, key=lambda x: -pa[x][lab])
    claim(set(top("Severe/profound impairment")[:2]) == {"STXBP1", "CDKL5"}, "severe/profound most often STXBP1 and CDKL5")
    claim(set(top("Autism")[:3]) == {"SYNGAP1", "EEF1A2", "PCDH19"}, "autism most often SYNGAP1, EEF1A2, PCDH19")
    claim(set(top("Movement disorder")[:2]) == {"FOXG1", "STXBP1"}, "movement disorders most often FOXG1, STXBP1")
    claim(set(top("Microcephaly")[:2]) == {"FOXG1", "SMC1A"}, "microcephaly most often FOXG1, SMC1A")
    claim(top("Cortical visual impairment")[0] == "CDKL5", "CVI most often CDKL5")
    claim(top("Deceased")[0] == "KCNT1", "deaths most often KCNT1")
    claim(set(sorted(dn, key=lambda x: -dn[x])[:3]) == {"CDKL5", "STXBP1", "SYNGAP1"}, "de novo highest CDKL5, STXBP1, SYNGAP1")
    claim(inh["paternal_by_gene"]["PCDH19"] > inh["maternal_by_gene"]["PCDH19"], "PCDH19 paternal > maternal")
    claim(sum(1 for x in inh["parental_mosaic_by_gene"].values() if x == 0) == 1, "mosaic parents in all but one gene")
    return v


CLAIMS = []  # wording that depends on the data (direction/significance); printed as CHECK when no longer true


def claim(ok, msg):
    if not ok:
        CLAIMS.append(f"text claim no longer supported by the data: {msg}")


def ptex(p):
    p = float(p)
    return r"P \textless{} 0.001" if p < 0.001 else pnum(p)


def n(x):
    return f"{int(x):,}"


def f1(x):
    return f"{float(x):.1f}"


def ci(c, d=2):
    return f"{c[0]:.{d}f}--{c[1]:.{d}f}"


def frac(s):
    """'213/259 (82.2%)' -> (213, 259, '82.2')"""
    m = re.match(r"(\d+)/(\d+) \(([\d.]+)%\)", s)
    return int(m.group(1)), int(m.group(2)), m.group(3)


def lead(s):
    """'82 (85.4%)' or '5.0 (3.0-6.5); n=89' -> leading number as string"""
    return re.match(r"[\d.]+", str(s)).group(0)


def and_join(items):
    return ", ".join(items[:-1]) + " and " + items[-1] if len(items) > 1 else items[0]


NUMWORDS = ("zero one two three four five six seven eight nine ten eleven twelve thirteen fourteen fifteen sixteen seventeen "
         "eighteen nineteen").split()
TENS = "twenty thirty forty fifty sixty seventy eighty ninety".split()


def num_words(k):
    k = int(k)
    if k < 20:
        return NUMWORDS[k]
    if k < 100:
        return TENS[k // 10 - 2] + ("-" + NUMWORDS[k % 10] if k % 10 else "")
    return f"{k:,}"


def text_values():
    """All other data-derived numbers in main.tex (flow, redundancy, landscape, treatment, Korean)."""
    v = {}
    F = R["flow"]
    T = F["Total"]
    t1 = R["table1"]
    # ---- dataset / flow
    v.update(SCREENED=n(T["screened"]), PROCESSED=n(T["processed"]), WITH_INDIV=n(T["with_individual_patients"]),
             DESC_N=n(T["individual_descriptions"]), ORIG_N=n(T["original_descriptions"]),
             REREP_EXTR=n(T["flagged_rereport_at_extraction"]), UNCLEAR_N=n(T["unclear_originality"]),
             AGG_N=n(T["aggregate_only"]), EXCL_N=n(T["original_excluded_insufficient_or_manual"]),
             CANON_N=n(T["original_in_canonical"]), N_PATIENTS=n(T["unique_patients"]), N_ROWS=n(T["variant_rows"]),
             N_ARTICLES=n(t1["Total"]["Articles"]), MERGED_XP=n(T["merged_other_paper"]))
    claim(11500 <= T["unique_patients"] < 12000, f"title says 'nearly 12,000' patients, now {T['unique_patients']}")
    v["N_BY_GENE"] = and_join([f"{g} {n(F[g]['unique_patients'])}" for g in GENES])
    s1 = pd.read_csv(TAB / "supp_data_S1_source_articles.csv", dtype=str, encoding="utf-8-sig")
    claim(len(s1) == int(t1["Total"]["Articles"]), f"Data S1 lists {len(s1)} articles but Table 1 counts {t1['Total']['Articles']}")
    fc = R["field_completeness_pct"]["Total"]
    for k, lab in [("PHENO", "Phenotype category"), ("SZ", "Seizure presence"), ("PROTEIN", "Protein change (reported or translated)"),
                   ("SEX", "Sex"), ("ONSET", "Onset age (parsed)"), ("FU", "Age at last follow-up"), ("ASM", "ASM response"),
                   ("TT", "Targeted treatment response")]:
        v[f"FC_{k}"] = f1(fc[lab])
    # ---- redundancy
    rd = R["redundancy"]
    v.update(RED_N=n(rd["rereport_total"]), RED_PCT=f1(rd["rereport_pct"]))
    claim(100 / 3 - 3.5 <= rd["rereport_pct"] <= 100 / 3 + 3.5, f"'one third / one in three' re-described, now {rd['rereport_pct']}%")
    bg = rd["by_gene_rereport_pct"]
    claim(max(bg, key=bg.get) == "EEF1A2", "redundancy highest for EEF1A2")
    lo = min(GENES, key=lambda x: bg[x])
    v.update(RED_MIN=f1(bg[lo]), RED_MIN_G=lo, RED_EEF1A2=f1(bg["EEF1A2"]), RED_SMC1A=f1(bg["SMC1A"]), RED_GABRA1=f1(bg["GABRA1"]))
    claim(sorted(GENES, key=lambda x: -bg[x])[:3] == ["EEF1A2", "SMC1A", "GABRA1"], "redundancy highest for EEF1A2, SMC1A, GABRA1")
    rv = pd.DataFrame(R["recurrent_variants"])
    fold = rv.inflation.max()
    v["MAX_FOLD_WORD"] = num_words(round(fold)) + "fold"
    t2 = rv.sort_values("inflation", ascending=False).head(2)
    v["TOP_INFLATION"] = and_join([f"{g} {p} ({nv} descriptions, {d} patients)" for g, p, nv, d in zip(t2.gene, t2.protein, t2.naive, t2.dedup)])
    claim(rv.loc[rv.inflation.idxmax(), "gene"] in ("EEF1A2", "SMC1A", "GABRA1"), "largest inflation in a review-dominated rare gene")
    ph = {r["gene"]: r for r in R["phenotype_naive_vs_dedup"]}
    big = [g for g in GENES if F[g]["unique_patients"] > 500 and g != "SYNGAP1"]
    v["DEE_PP_MAX4"] = f1(max(abs(ph[g]["naive_DEE_pct"] - ph[g]["dedup_DEE_pct"]) for g in big))
    v.update(EEF1A2_NAIVE=f1(ph["EEF1A2"]["naive_DEE_pct"]), EEF1A2_DEDUP=f1(ph["EEF1A2"]["dedup_DEE_pct"]),
             EEF1A2_DEE_PP=f"{ph['EEF1A2']['naive_DEE_pct'] - ph['EEF1A2']['dedup_DEE_pct']:.0f}")
    # ---- landscape: onset
    v["P_KW"] = ptex(R["onset_kruskal_p"])
    claim(R["onset_kruskal_p"] < 0.05, "onset differed between genes (Kruskal-Wallis)")
    oq = R["onset_quartiles_months"]

    def med(g):
        m = oq[g][1]
        return f"{m * 30.44:.0f} days" if m < 1 else f"{m:.0f} months" if m >= 12 else f"{m:.1f} months"
    early = sorted((g for g in GENES if oq[g][1] < 1), key=lambda g: oq[g][1])
    late = sorted((g for g in GENES if oq[g][1] >= 9), key=lambda g: oq[g][1])
    v["ONSET_SENTENCE"] = ("Median onset was in the neonatal period for " + and_join([f"{g} ({med(g)})" for g in early])
                           + ", in infancy for most other genes, and latest for " + and_join([f"{g} ({med(g)})" for g in late])
                           + "; SCN2A onset was bimodal, with neonatal and infantile peaks.")
    # ---- landscape: phenotype, sex, consequence
    ps, cs, fem = R["phenotype_share"], R["consequence_share"], R["female_pct"]
    # Table 1 uses only patients with a recorded phenotype category.
    dee_pct = {}
    for g in GENES:
        cell = t1[g]["DEE label, n/N (%)"]
        counts = re.match(r"([\d,]+)/([\d,]+)", cell)
        numerator, denominator = (int(x.replace(",", "")) for x in counts.groups())
        dee_pct[g] = 100 * numerator / denominator if denominator else 0
        assert f1(dee_pct[g]) == re.search(r"\(([\d.]+)\)", cell).group(1), f"DEE Table 1 mismatch: {g}"
    dee = sorted((x for x in GENES if dee_pct[x] > 0), key=lambda x: -dee_pct[x])
    hi = [x for x in dee if dee_pct[x] >= 45]
    lo = [x for x in dee if dee_pct[x] < 20]
    v["DEE_SHARE_LIST"] = (and_join([rf"{f1(dee_pct[x])}\% of {x}" + (" (predominantly Dravet syndrome)" if x == "SCN1A" else "") for x in hi])
                           + " patients but only " + and_join([rf"{f1(dee_pct[x])}\% of {x}" for x in lo]) + " patients")
    for g in GENES:
        v[f"PH_{g}_SL"] = f1(ps[g]["Self-limited (familial) neonatal/infantile epilepsy"])
        v[f"PH_{g}_NOSZ"] = f1(ps[g]["No seizures"])
        for k, lab in [("MIS", "Missense"), ("TRUNC", "Truncating"), ("CNV", "CNV / exon-level")]:
            v[f"CONS_{g}_{k}"] = f1(cs[g][lab])
    v["PH_SYNGAP1_GEN"] = f1(ps["SYNGAP1"]["Generalized epilepsy"])
    v.update(FEM_CDKL5=f1(fem["CDKL5"]), FEM_PCDH19=f1(fem["PCDH19"]), FEM_SMC1A=f1(fem["SMC1A"]))
    kc = R["kcnt1_categories_pct"]
    v.update(KCNT1_EIMFS_PCT=f1(kc["EIMFS"]), KCNT1_SHE_PCT=f1(kc["SHE/ADNFLE"]), PH_SCN8A_SL=f1(ps["SCN8A"]["Self-limited (familial) neonatal/infantile epilepsy"]))
    claim(ps["PCDH19"]["Other epilepsy"] == max(ps["PCDH19"].values()), "PCDH19 mostly clustering (other) epilepsy")
    for g in ["SYNGAP1", "CDKL5", "SMC1A"]:
        claim(cs[g]["Truncating"] > cs[g]["Missense"], f"truncating variants predominate in {g}")
    for g in ["EEF1A2", "SCN2A", "KCNQ2", "SCN8A", "KCNT1", "GABRA1"]:
        claim(cs[g]["Missense"] > 50, f"missense variants predominate in {g}")
    claim(cs["FOXG1"]["CNV / exon-level"] >= 20, "copy-number events frequent in FOXG1")
    # ---- DEE by consequence (+ clustered)
    cl = R["clustered_sensitivity"]
    for g in DBC_GENES:
        d, c = R["dee_by_consequence"][f"{g}|Reported DEE label"], cl[f"{g}: DEE, truncating vs missense"]
        v.update({f"DBC_{g}_T": frac(d["truncating"])[2], f"DBC_{g}_M": frac(d["missense"])[2], f"DBC_{g}_OR": f"{d['OR']:.2f}",
                  f"DBC_{g}_CI": ci(d["CI"]), f"GEE_{g}_OR": f"{c['GEE_OR']:.2f}", f"GEE_{g}_CI": ci(c["GEE_CI"])})
        claim(d["p"] < 0.05 and c["GEE_p"] < 0.05, f"{g} truncating vs missense DEE association significant, crude and clustered")
        claim((d["OR"] > 1) == (g in ("SCN1A", "STXBP1")) and (c["GEE_OR"] > 1) == (g in ("SCN1A", "STXBP1")),
              f"{g} direction of truncating-DEE association")
    # ---- treatment
    sc = R["scb"]
    e, l = sc["SCN2A onset <3 m"], sc["SCN2A onset ≥3 m"]
    b, w = R["scb_scn2a_early_vs_late_benefit"], R["scb_scn2a_early_vs_late_worsening"]
    cb, cw = cl["SCN2A: SCB benefit, onset <3 m vs >=3 m"], cl["SCN2A: SCB worsening, onset >=3 m vs <3 m"]
    v.update(S2_E_N=n(e["n"]), S2_L_N=n(l["n"]), S2_E_B_PCT=pct_of(e["Benefit"]), S2_L_B_PCT=pct_of(l["Benefit"]),
             S2_E_W_PCT=pct_of(e["Worsening"]), S2_L_W_PCT=pct_of(l["Worsening"]),
             S2_B_OR=f"{b['OR']:.1f}", S2_B_CI=ci(b["CI"], 1), S2_W_OR=f"{w['OR']:.1f}", S2_W_CI=ci(w["CI"], 1),
             S2_GEE_B_OR=f"{cb['GEE_OR']:.1f}", S2_GEE_B_CI=ci(cb["GEE_CI"], 1),
             S2_GEE_W_OR=f"{cw['GEE_OR']:.1f}", S2_GEE_W_CI=ci(cw["GEE_CI"], 1))
    claim(b["p"] < 0.05 and w["p"] < 0.05 and b["OR"] > 1 and w["OR"] > 1, "SCN2A SCB benefit (early) and worsening (late) differences")
    claim(cb["GEE_p"] < 0.05 and cw["GEE_p"] < 0.05, "SCN2A SCB differences 'held after clustering by article'")
    for k, name in [("SL", "KCNQ2 self-limited"), ("DEE", "KCNQ2 DEE")]:
        v[f"KQ_{k}_B"], v[f"KQ_{k}_N"] = lead(sc[name]["Benefit"]), n(sc[name]["n"])
        claim(float(pct_of(sc[name]["Benefit"])) >= (100 if k == "SL" else 70), f"'all'/'most' {name} patients benefited from SCBs")
    s8, sd, sl = R["scb_scn8a_all"], sc["SCN8A DEE"], sc["SCN8A self-limited"]
    so = sc["SCN8A other/not specified"]
    assert sum(x["n"] for x in (sd, sl, so)) == s8["n"]
    for outcome in ("Benefit", "No benefit", "Worsening"):
        assert sum(int(lead(x[outcome])) for x in (sd, sl, so)) == int(lead(s8[outcome]))
    v.update(S8_B=lead(s8["Benefit"]), S8_OTHER_B=lead(so["Benefit"]), S8_OTHER_N=n(so["n"]))
    v["S8_REG_N"] = num_words(R["scb_scn8a_registry_n"]).capitalize()
    v.update(S8_N=n(s8["n"]), S8_B_PCT=pct_of(s8["Benefit"]), S8_W_PCT=pct_of(s8["Worsening"]), S8_NB_PCT=pct_of(s8["No benefit"]),
             S8_DEE_B=lead(sd["Benefit"]), S8_DEE_N=n(sd["n"]), S8_SL_B=lead(sl["Benefit"]), S8_SL_N=n(sl["n"]))
    claim(float(pct_of(s8["Benefit"])) > 60 and float(pct_of(sd["Benefit"])) > 60 and float(pct_of(sl["Benefit"])) > 60,
          "SCN8A SCB benefit in most patients, in DEE and self-limited")
    sv = R["source_verification"]
    allv = list(sv.values())
    v.update(VER_ALL_PAIRS=n(sum(x["pairs"] for x in allv)),
             VER_ALL_AGREE=f1(sum(x["pairs"] * x["auto_agreement_pct"] for x in allv) / sum(x["pairs"] for x in allv)),
             VER_ALL_WRONG=n(sum(x["db_text_wrong"] for x in allv)))
    for k, key in [("VER8", "Sodium channel blocker|SCN8A"), ("VERQ", "Quinidine|KCNT1")]:
        v.update({f"{k}_PAIRS": n(sv[key]["pairs"]), f"{k}_AGREE": f1(sv[key]["auto_agreement_pct"]),
                  f"{k}_LABELS": esc(", ".join(f"{a} {b}" for a, b in sv[key]["labels"].items()))})
    ot = R["other_treatments"]
    q = ot["Quinidine|KCNT1"]
    v.update(QND_B=str(q["Benefit"]), QND_N=str(q["n"]), QND_NB=str(q["No benefit"]), QND_W=str(q["Worsening"]),
             QND_PCT=f"{100 * q['Benefit'] / q['n']:.0f}")
    claim(q["Benefit"] / q["n"] < 0.5, "quinidine benefit in 'fewer than half' of KCNT1 patients")
    for k, key in [("STP", "Stiripentol|SCN1A"), ("FFA", "Fenfluramine|SCN1A"), ("CBD", "Cannabidiol|SCN1A"),
                   ("KD", "Ketogenic diet|SCN1A"), ("KD_CDKL5", "Ketogenic diet|CDKL5"), ("KD_KCNT1", "Ketogenic diet|KCNT1"),
                   ("VGB_STXBP1", "Vigabatrin/ACTH/steroids|STXBP1")]:
        v[f"OT_{k}_B"], v[f"OT_{k}_N"], v[f"OT_{k}_W"] = str(ot[key]["Benefit"]), str(ot[key]["n"]), str(ot[key]["Worsening"])
    kd = ot["Ketogenic diet|CDKL5"]
    v["OT_KD_CDKL5_PCT"] = f"{100 * kd['Benefit'] / kd['n']:.1f}"
    claim(max((x for x in ot if x.endswith("|CDKL5")), key=lambda x: ot[x]["n"]) == "Ketogenic diet|CDKL5",
          "ketogenic diet the CDKL5 therapy most often reported with an individual response")
    ea = R["extraction_accuracy"]
    bg = ea["rows_added_pct_by_gene"]
    v.update(ACC_CHANGED_PCT=f1(100 - ea["field_unchanged_pct"]), ACC_UNCHANGED_PCT=f1(ea["field_unchanged_pct"]), ACC_VALUES=n(ea["field_values"]), ACC_ADDED_PCT=f1(ea["rows_added_pct"]),
             ACC_ADDED=n(ea["rows_added_on_review"]), ACC_FINAL=n(ea["final_rows"]), ACC_LLM_ROWS=n(ea["rows_llm"]),
             ACC_ADDED_MIN=f1(min(bg.values())), ACC_ADDED_MAX=f1(max(bg.values())),
             ACC_ARTICLES=n(ea["subset_articles"]), ACC_FOUND=n(ea["subset_patients_found"]), ACC_PATIENTS=n(ea["subset_patients_checked"]),
             ACC_NOTINPAPER=n(ea["subset_not_in_paper"]))
    claim(15 <= ea["rows_added_pct"] <= 25, "'about one in five' final records missed by the model")
    return v


def wilson_s(k, nn):
    import math
    if not nn:
        return "--"
    z, ph = 1.96, k / nn
    den = 1 + z * z / nn
    c = (ph + z * z / (2 * nn)) / den
    h = z * math.sqrt(ph * (1 - ph) / nn + z * z / (4 * nn * nn)) / den
    return f"{100 * (c - h):.1f}--{100 * (c + h):.1f}"


def review_values():
    """Numbers added in the revision after review (Methods, Results, Supporting tables)."""
    v = {}
    q, vq, crf, ea = R["review_fixes"], R["variant_qc"], R["crf_check"], R["extraction_accuracy"]
    sw = q["software"]
    v.update(ABS_ONLY_ART=num_words(q["abstract_only_articles"]) if q["abstract_only_articles"] < 10 else n(q["abstract_only_articles"]),
             ABS_ONLY_PAT=n(q["abstract_only_patients"]),
             CRF_N=n(crf["values"]), CRF_OK=n(crf["confirmed"]), CRF_NC=n(crf["not_confirmed"]), CRF_NC_PCT=f"{crf['not_confirmed_pct']:.0f}",
             ACC_NE_PCT=f1(ea["nonempty_changed_pct"]), ACC_NE_N=n(ea["nonempty_values"]),
             VQ_ROWS=n(vq["rows_checked"]), VQ_REL=n(vq["relabelled"]), VQ_NOT=n(vq["notation"]), VQ_INC=n(vq["inconsistent"]),
             VQ_INC_PCT=f1(vq["inconsistent_pct"]), BLB_EXCL=n(q["blb_excluded"]), INC_CHECKED=n(q["inc_checked"]), INC_EXCL=n(q["inc_excluded"]),
             DUP_PAIRS=n(q["dup_pairs"]), DUP_SAME=n(q["dup_same"]), DUP_DIFF=n(q["dup_different"]), DUP_UNC=n(q["dup_uncertain"]),
             N_FAMILIES=n(q["families"]), SW_PY=sw["python"], SW_PD=sw["pandas"], SW_NP=sw["numpy"], SW_SP=sw["scipy"],
             SW_SM=sw["statsmodels"], SW_MPL=sw["matplotlib"])
    # FOXG1 CNV type
    cnv = pd.read_csv(TAB / "supp_consequence_by_gene.csv", dtype=str) if (TAB / "supp_consequence_by_gene.csv").exists() else None
    v["FOXG1_CNV_MAIN"] = R.get("foxg1_cnv_main", "deletions and duplications")
    # DEE label / developmental impairment by consequence
    d = R["dee_by_consequence"]
    for g in ["SCN1A", "SCN2A", "KCNQ2", "STXBP1"]:
        x = d[f"{g}|Reported DEE label"]
        a, b_ = frac(x["truncating"])[:2]
        c, dd = frac(x["missense"])[:2]
        v.update({f"DBC_{g}_T": frac(x["truncating"])[2], f"DBC_{g}_M": frac(x["missense"])[2]})
        claim((a / b_ > c / dd) == (g in ("SCN1A", "STXBP1")), f"{g} DEE label direction")
    x = d["SCN2A|Developmental impairment (patients with epilepsy)"]
    a, na, pa = frac(x["truncating"])
    c, nc, pc = frac(x["missense"])
    import math
    if 0 in (a, na - a, c, nc - c):
        aa, bb, cc, dd = a + .5, na - a + .5, c + .5, nc - c + .5
    else:
        aa, bb, cc, dd = a, na - a, c, nc - c
    orr = aa * dd / (bb * cc)
    lo, hi = x["cluster_CI"]
    v.update(DEV_SCN2A_T=pa, DEV_SCN2A_M=pc, DEV_SCN2A_OR=f"{orr:.1f}", DEV_SCN2A_CI=f"{lo:.1f}--{hi:.1f}")
    claim(orr > 1 and lo > 1, "SCN2A truncating carriers more often developmental impairment (CI above 1)")
    # seizure types extras
    st = R["seizure_types"]
    fx = st["foxg1_spasms"]
    v.update(ST_FOXG1_GAIN=fx["gain_carriers_with_spasms"].replace("/", " of "), ST_FOXG1_NOGAIN=f"{fx['without_gains']:.0f}",
             ST_KCNQ2_TO_DEE=f"{st['kcnq2_tonic_by_phenotype']['DEE']:.0f}",
             ST_KCNQ2_TO_SL=f"{st['kcnq2_tonic_by_phenotype']['Self-limited (familial) neonatal/infantile epilepsy']:.0f}")
    # treatment CIs and sensitivity
    sc = R["scb"]
    k1 = sc["SCN1A"]
    v["S1_W_CI"] = wilson_s(int(lead(k1["Worsening"])), k1["n"])
    lt = R["scn1a_scb_by_drug"]["lamotrigine"]
    v.update(S1_LTG_W=n(lt["Worsening"]), S1_LTG_N=n(lt["n"]))
    s8 = R["scb_scn8a_all"]
    v["S8_B_CI"] = wilson_s(int(lead(s8["Benefit"])), s8["n"])
    reg, nreg = R["scb_scn8a_registry"], R["scb_scn8a_no_registry"]
    v.update(S8_REG_B=n(reg["Benefit"]), S8_REG_W=n(reg["Worsening"]), S8_NOREG_PCT=pct_of(nreg["Benefit"]))
    e, l = sc["SCN2A onset <3 m"], sc["SCN2A onset ≥3 m"]
    v.update(S2_E_B_CI=wilson_s(int(lead(e["Benefit"])), e["n"]), S2_L_B_CI=wilson_s(int(lead(l["Benefit"])), l["n"]))
    so = R["scb_scn2a_onset_or"]
    v.update(S2_WOLFF_N=n(so["wolff_pairs"]), S2_NOWOLFF_OR=f"{so['without_wolff2017']['OR']:.1f}",
             S2_NOWOLFF_CI=ci(so["without_wolff2017"]["CI"], 1))
    qd = R["other_treatments"]["Quinidine|KCNT1"]
    v["QND_CI"] = wilson_s(qd["Benefit"], qd["n"])
    v["MERGED_WP"] = n(R["flow"]["Total"].get("merged_within_paper", 0))
    ss = R["small_study"]
    claim(all(d[">10"]["Worsening"] / d[">10"]["n"] > d["1-3"]["Worsening"] / d["1-3"]["n"] for d in ss.values()),
          "worsening more often in large series than in 1-3 patient articles, all three genes")
    return v


def fill(text):
    for k, v in {**values(), **text_values(), **new_values(), **cm_values(), **review_values()}.items():
        text = text.replace("{{" + k + "}}", v)
    return text


def fill_supp(text, frozen_tables=False):
    if frozen_tables:
        s8 = pd.read_csv(TAB / "supp_S8_change_log_summary.csv")
    else:
        log = pd.read_csv(CHANGE_LOG, dtype=str, keep_default_na=False)
        s8 = log.groupby(["change_type", "field"]).size().rename("n").reset_index().sort_values("n", ascending=False)
        s8.to_csv(TAB / "supp_S8_change_log_summary.csv", index=False)
    s10 = pd.read_csv(TAB / "supp_S10_clustered_sensitivity.csv", dtype=str, keep_default_na=False)
    s10["Fisher P"] = [pnum(p) for p in s10["Fisher P"]]
    s10["Cluster-robust P"] = [pnum(p) for p in s10["Cluster-robust P"]]
    text = fill(text)
    text = text.replace("{{S8_TOTAL}}", f"{int(s8.n.sum()):,}")
    text = text.replace("{{S7A}}", latex_table(pd.read_csv(TAB / "supp_S7a_extraction_accuracy.csv", dtype=str), ""))
    text = text.replace("{{S7B}}", latex_table(pd.read_csv(TAB / "supp_S7b_adjudicated_subset.csv", dtype=str), "", size=r"\footnotesize"))
    text = text.replace("{{S8}}", latex_table(s8.rename(columns={"change_type": "Change type", "field": "Field"}), ""))
    text = text.replace("{{S10}}", latex_table(s10, "", size=r"\footnotesize"))
    text = text.replace("{{S7C}}", latex_table(pd.read_csv(TAB / "supp_S7c_extraction_nonempty.csv", dtype=str), "", size=r"\footnotesize"))
    text = text.replace("{{CRF_TABLE}}", latex_table(pd.read_csv(TAB / "supp_crf_reconciliation.csv", dtype=str), ""))
    text = text.replace("{{VQ_TABLE}}", latex_table(pd.read_csv(TAB / "supp_variant_qc.csv", dtype=str), "", size=r"\footnotesize"))
    return text


SI_FILES = {
    "T:table1counts": [TAB / "table1_gene_summary.csv"],
    "T:crf": [TAB / "supp_crf_reconciliation.csv", TAB / "supp_crf_reconciliation_by_field.csv"],
    "T:accuracy": [TAB / "supp_S7a_extraction_accuracy.csv", TAB / "supp_S7c_extraction_nonempty.csv", TAB / "supp_S7b_adjudicated_subset.csv",
                   TAB / "supp_S7_llm_verification_errors.csv"],
    "T:changelog": [TAB / "supp_S8_change_log_summary.csv"], "T:variantqc": [TAB / "supp_variant_qc.csv"],
    "T:seizureval": [TAB / "supp_seizure_types_validation.csv", TAB / "supp_seizure_types_validation_dev.csv"],
    "T:comorbval": [TAB / "supp_comorbidity_validation.csv"], "T:treatval": [],
    "T:clustered": [TAB / "supp_S10_clustered_sensitivity.csv"], "T:completeness": [TAB / "supp_field_completeness.csv"],
    "T:phenotype": [TAB / "supp_phenotype_by_gene.csv", TAB / "supp_consequence_by_gene.csv"], "T:dbc": [TAB / "supp_dee_by_consequence.csv"],
    "T:seizure": [TAB / "supp_seizure_types_by_gene.csv"], "T:comorbidity": [TAB / "supp_comorbidity_by_gene.csv"],
    "T:deaths": [TAB / "supp_deaths_by_gene.csv"], "T:inheritance": [TAB / "supp_inheritance_by_gene.csv"],
    "T:candidates": [TAB / "supp_variants_reclassification_candidates.csv"], "T:authorvs": [TAB / "supp_variants_author_vs_clinvar.csv"],
    "T:clinvarall": [TAB / "supp_variants_clinvar_all.csv"],
    "T:treatment": [TAB / "supp_S9_treatment_primary.csv", TAB / "supp_S9b_treatment_with_responsive_field.csv"],
    "T:smallstudy": [TAB / "supp_small_study.csv"], "T:recurrent": [TAB / "supp_recurrent_variants_naive_vs_dedup.csv"],
    "T:deenaive": [TAB / "supp_phenotype_naive_vs_dedup.csv"],
}
SI_NUM = {}


def number_si(tex):
    """Number Supporting Information items in order of first citation in the main text: [[T:key]] -> Table Sn, [[A:key]] -> Appendix Sn,
    [[F:key]] -> Figure Sn (reviewer comment R4: SI numbered in citation order)."""
    body = tex.split(r"\section{Introduction}", 1)[1].split(r"\bmsection*{Author Contributions}", 1)[0]
    order = {"T": [], "A": [], "F": []}
    for kind, key in re.findall(r"\[\[([TAF]):(\w+)\]\]", body + tex):
        if key not in order[kind]:
            order[kind].append(key)
    word = {"T": "Table", "A": "Appendix", "F": "Figure"}
    for kind, keys in order.items():
        for i, key in enumerate(keys, start=1):
            SI_NUM[f"{kind}:{key}"] = i
            tex = tex.replace(f"[[{kind}:{key}]]", f"{word[kind]} S{i}")
    # Supporting-information list in numerical order
    a = tex.index(r"\begin{itemize}", tex.index(r"\bmsection*{Supporting information}"))
    b = tex.index(r"\end{itemize}", a)
    items = re.findall(r"\\item\s+\\textbf\{(.+?)\.\}\s(.+?)(?=\n\\item|\n\\end|$)", tex[a:b], flags=re.S)
    rank = lambda lab: ({"Appendix": 0, "Table": 1, "Figure": 2, "Data": 3}[lab.split()[0]], int(re.sub(r"\D", "", lab) or 0))
    items = sorted(items, key=lambda x: rank(x[0]))
    lst = "\\begin{itemize}\n\\tightlist\n" + "".join(f"\\item\n  \\textbf{{{lab}.}} {d.strip()}\n" for lab, d in items)
    tex = tex[:a] + lst + tex[b:]
    rows = [f"S{SI_NUM[k]} & {esc(d.strip())} \\\\" for lab, d in items for k in SI_NUM if k.startswith("T:") and lab == f"Table S{SI_NUM[k]}"]
    supp_list = "\n".join([r"\begin{longtable}{lp{0.8\linewidth}}", r"\toprule", r"\textbf{Table} & \textbf{Content} \\", r"\midrule",
                           *rows, r"Data S1 & The source articles with bibliographic details, genes and patients contributed per gene \\",
                           r"\bottomrule", r"\end{longtable}"])
    return tex.replace("{{SUPP_TABLE_LIST}}", supp_list)


def plain(tex):
    return subprocess.run(["pandoc", "-f", "latex", "-t", "plain", "--wrap=none"], input=tex,
                          capture_output=True, text=True, check=True).stdout


def braced(tex, start):
    """Content of the balanced {...} group that follows `start`."""
    i = tex.index(start) + len(start)
    depth, j = 1, i
    while depth:
        depth += {"{": 1, "}": -1}.get(tex[j], 0)
        j += 1
    return tex[i:j - 1]


def njd_parts(tex):
    body_tex = tex.split(r"\begin{document}", 1)[1]
    return {
        "title": braced(tex, r"\title{"),
        "abstract": braced(tex, r"\abstract[Abstract]{"),
        "keywords": braced(tex, r"\keywords{"),
        "keypoints": body_tex.split(r"\begin{keypointsbox}", 1)[1].split(r"\end{keypointsbox}", 1)[0],
        "body": body_tex.split(r"\end{keypointsbox}", 1)[1].split(r"\bmsection*{Author Contributions}", 1)[0],
        "back": r"\bmsection*{Author Contributions}" + body_tex.split(r"\bmsection*{Author Contributions}", 1)[1].split(r"\bibliography{", 1)[0],
        "suppinfo": body_tex.split(r"\bmsection*{Supporting information}", 1)[1].split("% SUPPLEMENTARY", 1)[0],
        "supp": body_tex.split("% SUPPLEMENTARY", 1)[1].split(r"\end{document}", 1)[0],
    }


def checks(tex):
    problems = []
    part = njd_parts(tex)
    aw = len(re.sub(r"\b(Objective|Methods|Results|Significance):", "", plain(part["abstract"])).split())
    if aw > 300:
        problems.append(f"abstract {aw} words > 300")
    kp = [plain(x).strip() for x in re.findall(r"\\item\s+(.+)", part["keypoints"])]
    for k in kp:
        if len(k) > 140:
            problems.append(f"key point {len(k)} chars > 140: {k[:50]}...")
    if not 3 <= len(kp) <= 5:
        problems.append(f"{len(kp)} key points")
    body = re.sub(r"\\begin\{(figure|table)\*\}.*?\\end\{\1\*\}", "", part["body"], flags=re.S)
    body = re.sub(r"\\(sub)?section\{[^}]*\}", "", body)
    words = len(plain(body).split())
    if words > 4000:
        problems.append(f"main text {words} words > 4000")
    left = [x for x in re.findall(r"\{\{\w+\}\}", tex) if x != "{{WORDS}}"]
    if left:
        problems.append(f"unfilled placeholders: {sorted(set(left))}")
    return aw, len(kp), words, problems


def word_files(tex, words):
    """Epilepsia wants Word: plain title page, Key Points, abstract, text, references, legends and tables at the end,
    figures uploaded separately; the Supporting Information as a second file."""
    part = njd_parts(tex)
    body = part["body"]
    figs = re.findall(r"\\begin\{figure\*\}.*?\\end\{figure\*\}", body, flags=re.S)
    legends = [(n, braced(f, r"\caption{")) for n, f in enumerate(figs, start=1)]
    tables = re.findall(r"\\begin\{table\*\}.*?\\end\{table\*\}", body, flags=re.S)
    body = re.sub(r"\\begin\{(figure|table)\*\}.*?\\end\{\1\*\}", "", body, flags=re.S)
    body = body.replace(r"\section{", r"\section*{").replace(r"\subsection{", r"\subsection*{")
    back = part["back"].replace(r"\bmsection*{", r"\section*{")
    table_tex = "\n\n".join(t.replace(r"\caption{" + braced(t, r"\caption{") + "}", r"\textbf{Table " + str(i) + ".} " + braced(t, r"\caption{") + "\n")
                              .replace(r"\begin{table*}[t]", "").replace(r"\end{table*}", "")
                              .replace(r"\begin{tablenotes}", "").replace(r"\end{tablenotes}", "").replace(r"\item ", "")
                              .replace(r"\resizebox{\textwidth}{!}{", "").replace(r"\end{tabular}}", r"\end{tabular}")
                              for i, t in enumerate(tables, start=1))
    head = (rf"\title{{{part['title']}}}" + "\n\\begin{document}\n\\maketitle\n\n"
            + r"\textbf{Running title:} Patient-level atlas of DEE genes" + "\n\n"
            + r"\textbf{Authors:} Sanghyun Lee\textsuperscript{1} (s4n9hyun@snu.ac.kr), Byung Chan Lim\textsuperscript{2,*}, Kwangsoo Kim\textsuperscript{3,*}" + "\n\n"
            + r"\textsuperscript{1}Interdisciplinary Program in Bioinformatics, Seoul National University, Seoul, Republic of Korea; \textsuperscript{2}Department of Pediatrics, Seoul National University Children's Hospital, Seoul, Republic of Korea; \textsuperscript{3}Transdisciplinary Department of Medicine and Advanced Technology, Seoul National University Hospital, Seoul, Republic of Korea. *Co-corresponding authors." + "\n\n"
            + r"\textbf{Corresponding authors:} Byung Chan Lim, MD, PhD, Department of Pediatrics, Seoul National University Children's Hospital, Seoul, Republic of Korea [TODO: address, email, ORCID]; Kwangsoo Kim, Transdisciplinary Department of Medicine and Advanced Technology, Seoul National University Hospital, Seoul, Republic of Korea [TODO: address, email, ORCID]" + "\n\n"
            + rf"\textbf{{Keywords:}} {part['keywords'].replace(', ', '; ')}" + "\n\n"
            + rf"\textbf{{Counts:}} Words (main text): {words:,} · References: {(TEX / 'references.bib').read_text().count('@article')} · Figures: {len(legends)} · Tables: {len(tables)}" + "\n\n"
            + r"\section*{Key Points}" + "\n" + part["keypoints"] + "\n\n"
            + r"\section*{Abstract}" + "\n" + part["abstract"].replace(r"\par\smallskip", "\n") + "\n\n")
    main_tex = head + body + back + "\n\\end{document}\n"
    md = subprocess.run(["pandoc", "-f", "latex", "-t", "markdown", "--wrap=none"], input=main_tex, capture_output=True, text=True, check=True).stdout
    tail_tex = ("\\begin{document}\n" + r"\section*{Supporting Information}" + part["suppinfo"]
                + "\n\n" + r"\section*{Figure Legends}" + "\n\n"
                + "\n\n".join(rf"\textbf{{Figure {n}.}} {c.strip()}" for n, c in legends)
                + "\n\n" + r"\section*{Tables}" + "\n\n" + table_tex + "\n\\end{document}\n")
    tail_md = subprocess.run(["pandoc", "-f", "latex", "-t", "markdown", "--wrap=none"], input=tail_tex, capture_output=True, text=True, check=True).stdout
    md = md + "\n\n# References\n\n::: {#refs}\n:::\n\n" + tail_md
    (OUT / "word_main.md").write_text(md)
    subprocess.run(["pandoc", "word_main.md", "--citeproc", "--bibliography", "references.bib", "--csl", "vancouver.csl",
                    "-o", "manuscript_epilepsia.docx"], cwd=OUT, check=True)
    supp_tex = (r"\documentclass{article}\begin{document}" + "\n" + r"\section*{Supporting Information}" + "\n"
                + part["supp"].replace(r"\onecolumn", "").replace(r"\clearpage", "") + "\n\\end{document}\n")
    supp_tex = supp_tex.replace("figures/figS1_redundancy.pdf", "figures/figS1_redundancy.png")  # Word cannot show PDF images
    (OUT / "word_supp.tex").write_text(supp_tex)
    subprocess.run(["pandoc", "word_supp.tex", "-o", "supporting_information_epilepsia.docx"], cwd=OUT, check=True)
    (OUT / "word_main.md").unlink()
    (OUT / "word_supp.tex").unlink()


def supporting_workbook():
    out = OUT / "Supporting_Tables_Data_S1.xlsx"
    with pd.ExcelWriter(out) as w:
        for key, files in sorted(SI_FILES.items(), key=lambda kv: SI_NUM.get(kv[0], 999)):
            if key not in SI_NUM:
                continue
            files = [f for f in files if f.exists()]
            for j, f in enumerate(files):
                name = f"Table S{SI_NUM[key]}" + ("abcd"[j] if len(files) > 1 else "")
                pd.read_csv(f, dtype=str, keep_default_na=False, encoding="utf-8-sig").to_excel(w, sheet_name=name, index=False)
        pd.read_csv(TAB / "supp_data_S1_source_articles.csv", dtype=str, keep_default_na=False, encoding="utf-8-sig").to_excel(w, sheet_name="Data S1", index=False)


def main(frozen_tables=False):
    OUT.mkdir(exist_ok=True)
    tex = number_si(fill_supp((TEX / "main.tex").read_text(), frozen_tables=frozen_tables).replace("{{TABLE1}}", table1()).replace("{{TABLE2}}", table2()))
    aw, nkp, words, problems = checks(tex)
    tex = tex.replace("{{WORDS}}", f"{words:,}")
    left = sorted(set(re.findall(r"\{\{\w+\}\}", tex)))
    if left:
        raise SystemExit(f"render_manuscript: unfilled placeholders {left}")
    (OUT / "main.tex").write_text(tex)
    for f in ("references.bib", "vancouver.csl"):
        (OUT / f).write_text((TEX / f).read_text())
    word_files(tex, words)
    supporting_workbook()
    print(f"abstract {aw} words; key points {nkp}; main text {words} words")
    for p in problems + CLAIMS:
        print("CHECK:", p)


if __name__ == "__main__":
    main()
