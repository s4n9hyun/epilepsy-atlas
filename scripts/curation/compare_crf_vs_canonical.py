"""교수님 curated CRF 통합본과 현재 canonical patient-variant CSV를 행 단위로 비교한다.

- SCN1A/SCN2A/KCNQ2/CDKL5: CRF Patient ID가 canonical_patient_id(CP번호)이므로 CP ID + 변이로 매칭
- SYNGAP1/EEF1A2: CRF가 자체 ID를 쓰므로 PMID + 변이로 매칭
- 매칭된 행은 필드별로 same / changed / filled_by_crf / blanked_by_crf 로 분류
"""
import glob
import re
from collections import defaultdict
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "runs/crf_0904/crf_vs_current/validation_outputs"
CRF = ROOT / "feedback/20260923/유전자별_Phenotype_CRF_통합_0904.xlsx"
GENES = ["SCN1A", "SCN2A", "KCNQ2", "CDKL5", "SYNGAP1", "EEF1A2"]
CP_GENES = {"SCN1A", "SCN2A", "KCNQ2", "CDKL5"}
MISSING = {"", "not specified", "n/a", "na", "none", "nan", "not used", "unknown", "-"}

AA = {"A": "Ala", "R": "Arg", "N": "Asn", "D": "Asp", "C": "Cys", "Q": "Gln", "E": "Glu",
      "G": "Gly", "H": "His", "I": "Ile", "L": "Leu", "K": "Lys", "M": "Met", "F": "Phe",
      "P": "Pro", "S": "Ser", "T": "Thr", "W": "Trp", "Y": "Tyr", "V": "Val"}

# CRF 컬럼 -> canonical 컬럼
FIELDS = {
    "variant_c": ("Variant cDNA change", "variant_c"),
    "variant_p": ("Variant protein change", "variant_p"),
    "clinvar_pathogenicity": ("ClinVar pathogenicity", "clinvar_pathogenicity"),
    "clinvar_id": ("ClinVar ID", "ClinVar ID"),
    "sex": ("Sex", "sex"),
    "age_of_onset": ("Age at seizure onset", "age_of_onset"),
    "last_followup_age": ("Age at last follow-up", "last_followup_age"),
    "phenotype": ("Phenotype", "phenotype_review_category"),
    "targeted_response": ("Response to targeted ASM", "targeted_treatment_response"),
}


def missing(v):
    return str(v).strip().lower() in MISSING or pd.isna(v)


def norm_text(v):
    if missing(v):
        return ""
    s = str(v).strip().lower()
    s = re.sub(r"\.0$", "", s)
    s = re.sub(r"\s+", " ", s)
    s = re.sub(r"\b(month|year|day|week|hour)\b", r"\1s", s)
    s = s.replace("mo ", "months ").rstrip(".")
    return s


def norm_c(v):
    if missing(v):
        return ""
    s = str(v).strip()
    s = re.sub(r"^.*?(NM_\d+(\.\d+)?)\s*[:(]?\s*", "", s) if "NM_" in s else s
    s = re.sub(r"\s+", "", s).replace("->", ">").replace("−", "-")
    s = re.sub(r"^chr[0-9xy]{1,2}:", "", s, flags=re.I)  # 좌표 표기의 chr 접두어 유무는 같은 값으로 본다
    s = s.rstrip(")")
    if re.match(r"^\d", s):
        s = "c." + s
    return s.lower()


def norm_p(v):
    if missing(v):
        return ""
    s = re.sub(r"\s+", "", str(v)).replace("(", "").replace(")", "")
    if not s.lower().startswith("p."):
        s = "p." + s
    body = s[2:]
    m = re.match(r"^([A-Z])(\d+)([A-Z*X])$", body)
    if m:
        a, n, b = m.groups()
        b = "Ter" if b in "*X" else AA.get(b, b)
        body = f"{AA.get(a, a)}{n}{b}"
    body = body.replace("*", "Ter").replace("X", "Ter") if "fs" in body or body.endswith(("*", "X")) else body
    return ("p." + body).lower()


def norm_field(key, v):
    if key == "variant_c":
        return norm_c(v)
    if key == "variant_p":
        return norm_p(v)
    if key == "clinvar_id":
        return "" if missing(v) else re.sub(r"\.0$", "", str(v).strip())
    return norm_text(v)


def first_pmid(v):
    m = re.findall(r"\d{6,9}", str(v))
    return m[0] if m else ""


def load(gene):
    crf = pd.read_excel(CRF, sheet_name=gene, header=1, dtype=str)
    crf = crf[crf["Patient ID"].notna()].copy()
    crf["pmid_n"] = crf["Source in PMID"].map(first_pmid)
    crf["korean_note"] = crf["Source in PMID"].fillna("").str.contains("korea", case=False)
    crf["crf_row"] = range(3, len(crf) + 3)  # 엑셀 행 번호(제목행+헤더행 다음부터)
    path = sorted(glob.glob(str(ROOT / f"outputs/current/canonical_patients/{gene}_canonical_patients/*.csv")))[-1]
    can = pd.read_csv(path, dtype=str)
    can["pmid_n"] = can["pmid"].map(first_pmid)
    return crf, can


def variant_key(row, c_col, p_col):
    c, p = norm_c(row[c_col]), norm_p(row[p_col])
    return c or p


def match(gene, crf, can):
    """(crf_idx, can_idx) 쌍, crf-only idx, canonical-only idx를 반환."""
    if gene in CP_GENES:
        group_crf = crf.groupby("Patient ID").groups
        group_can = can.groupby("canonical_patient_id").groups
    else:
        group_crf = crf.groupby("pmid_n").groups
        group_can = can.groupby("pmid_n").groups
    pairs, crf_only, can_only = [], [], []
    for key in set(group_crf) | set(group_can):
        a = list(group_crf.get(key, []))
        b = list(group_can.get(key, []))
        # 1차: 변이(cDNA 또는 protein) 일치, 2차: 남은 행은 같은 그룹 안에서 순서대로 짝짓기(CP 그룹만)
        rest_b = b.copy()
        rest_a = []
        for i in a:
            ci, pi = norm_c(crf.at[i, "Variant cDNA change"]), norm_p(crf.at[i, "Variant protein change"])
            hit = None
            for j in rest_b:
                cj, pj = norm_c(can.at[j, "variant_c"]), norm_p(can.at[j, "variant_p"])
                if (ci and ci == cj) or (pi and pi == pj):
                    hit = j
                    break
            if hit is not None:
                pairs.append((i, hit, "variant_match"))
                rest_b.remove(hit)
            else:
                rest_a.append(i)
        if gene in CP_GENES:
            for i, j in zip(rest_a, rest_b):
                pairs.append((i, j, "same_cp_variant_differs"))
            crf_only += rest_a[len(rest_b):]
            can_only += rest_b[len(rest_a):]
        else:
            crf_only += rest_a
            can_only += rest_b
    return pairs, crf_only, can_only


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    summary, field_rows, diff_rows, only_rows = [], [], [], []
    mut_type = []
    for gene in GENES:
        crf, can = load(gene)
        pairs, crf_only, can_only = match(gene, crf, can)
        mut_type.append(crf["Mutation Type"].fillna("(blank)").str.strip().str.lower()
                        .value_counts().rename_axis("mutation_type").reset_index(name="rows").assign(gene=gene))
        stats = defaultdict(lambda: defaultdict(int))
        for i, j, how in pairs:
            for key, (ccol, kcol) in FIELDS.items():
                if ccol not in crf.columns:
                    continue
                a, b = norm_field(key, crf.at[i, ccol]), norm_field(key, can.at[j, kcol])
                if a == b:
                    status = "same"
                elif a and not b:
                    status = "filled_by_crf"
                elif b and not a:
                    status = "blanked_by_crf"
                else:
                    status = "changed"
                stats[key][status] += 1
                if status != "same":
                    diff_rows.append({
                        "gene": gene, "field": key, "status": status, "match_type": how,
                        "crf_patient_id": crf.at[i, "Patient ID"], "crf_excel_row": crf.at[i, "crf_row"],
                        "canonical_patient_id": can.at[j, "canonical_patient_id"],
                        "pmid": can.at[j, "pmid_n"], "mention_id": can.at[j, "mention_id"],
                        "crf_value": crf.at[i, ccol], "canonical_value": can.at[j, kcol],
                    })
        for key, s in stats.items():
            field_rows.append({"gene": gene, "field": key, **{k: s.get(k, 0) for k in
                               ["same", "changed", "filled_by_crf", "blanked_by_crf"]}})
        for i in crf_only:
            only_rows.append({"gene": gene, "side": "crf_only", "id": crf.at[i, "Patient ID"],
                              "crf_excel_row": crf.at[i, "crf_row"], "pmid": crf.at[i, "pmid_n"],
                              "source_in_pmid": crf.at[i, "Source in PMID"],
                              "variant_c": crf.at[i, "Variant cDNA change"],
                              "variant_p": crf.at[i, "Variant protein change"],
                              "mutation_type": crf.at[i, "Mutation Type"],
                              "id_in_current_canonical": crf.at[i, "Patient ID"] in set(can.canonical_patient_id)})
        for j in can_only:
            only_rows.append({"gene": gene, "side": "canonical_only", "id": can.at[j, "canonical_patient_id"],
                              "pmid": can.at[j, "pmid_n"], "mention_id": can.at[j, "mention_id"],
                              "variant_c": can.at[j, "variant_c"], "variant_p": can.at[j, "variant_p"],
                              "patient_id_in_paper": can.at[j, "patient_id"],
                              "pmid_in_crf": can.at[j, "pmid_n"] in set(crf.pmid_n)})
        summary.append({
            "gene": gene,
            "crf_rows": len(crf), "canonical_rows": len(can),
            "crf_unique_patients": crf["Patient ID"].nunique(),
            "canonical_unique_patients": can.canonical_patient_id.nunique(),
            "crf_pmids": crf.pmid_n.replace("", pd.NA).nunique(), "canonical_pmids": can.pmid_n.nunique(),
            "matched_rows": len(pairs),
            "matched_same_variant": sum(1 for p in pairs if p[2] == "variant_match"),
            "matched_same_cp_variant_differs": sum(1 for p in pairs if p[2] == "same_cp_variant_differs"),
            "crf_only_rows": len(crf_only), "canonical_only_rows": len(can_only),
            "crf_korean_cohort_rows": int(crf.korean_note.sum()),
        })
    pd.DataFrame(summary).to_csv(OUT / "summary_by_gene.csv", index=False)
    pd.DataFrame(field_rows).to_csv(OUT / "field_agreement_by_gene.csv", index=False)
    pd.DataFrame(diff_rows).to_csv(OUT / "matched_row_field_differences.csv", index=False)
    pd.DataFrame(only_rows).to_csv(OUT / "unmatched_rows.csv", index=False)
    pd.concat(mut_type).to_csv(OUT / "crf_mutation_type_counts.csv", index=False)
    print(pd.DataFrame(summary).to_string(index=False))
    print(pd.DataFrame(field_rows).to_string(index=False))


if __name__ == "__main__":
    main()
