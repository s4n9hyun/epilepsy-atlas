#!/usr/bin/env python3
"""
Create conservative canonical patient records from patient_mentions CSV.

This script does not delete the original mention rows. It creates separate
review files:
  - patient_mentions.csv: all mentions with canonical status columns added
  - canonical_patients.csv: patient-variant rows from original_patient rows
    (same canonical_patient_id can appear on multiple variant rows)
  - patient_dedup_links.csv: mapping from canonical patients to source mentions
  - excluded_mentions.csv: cited/aggregate/unclear mentions kept out of canonical DB
  - dedup_summary.json: summary and policy metadata
"""

import argparse
import json
import re
from pathlib import Path

import pandas as pd

from timestamped_outputs import make_timestamp, timestamp_path


MISSING_VALUES = {"", "not specified", "n/a", "na", "none", "unknown", "nan"}
GENES = {
    "SCN1A", "SCN2A", "CDKL5", "EEF1A2", "KCNQ2", "SYNGAP1",
    "PCDH19", "SCN8A", "GABRA1", "STXBP1", "SMC1A", "FOXG1", "KCNT1",
}

NEWLY_REPORTED_PATIENT_PHRASES = [
    # 개선 항목 8 (2026-09-24): 임상 보고는 이전 논문이어도 변이를 처음 보고한 논문이면 original
    "first genotype report",
    "reported for the first time",
    "first reported",
    "newly reported",
    "newly recruited",
    "newly diagnosed",
    "newly characterized",
    "newly characterised",
    "new case",
    "new cases",
    "new patient",
    "new patients",
    "new individual",
    "new individuals",
    "new subject",
    "new subjects",
    "novel patient",
    "novel patients",
]

PREVIOUSLY_REPORTED_PATIENT_PHRASES = [
    "previously described",
    "previously reported",
    "previously published",
    "re-analyzed",
    "re-analysed",
    "reanalyzed",
    "reanalysed",
    "reported by",
    "adapted from",
    "copied from",
    "from prior literature",
]

NEGATED_PREVIOUSLY_REPORTED_PHRASES = [
    "not flagged as previously",
    "not listed as previously",
    "not marked as previously",
    "not previously",
    "not reported previously",
    "not previously reported",
    "not previously described",
    'previously published as "-"',
    "previously published as '-'",
]

GENERIC_DIAGNOSIS_PHRASES = [
    "syngap1-associated syndrome",
    "syngap1 related intellectual disability",
    "syngap1-related intellectual disability",
    "syngap1-id",
    "syngap1 encephalopathy",
    "syngap1 developmental and epileptic encephalopathy",
    "syngap1-associated disorder",
    "syngap1-related disorder",
]

GENERIC_VARIANT_ONLY_PHRASES = [
    "microdeletion",
    "microduplication",
    "copy number variant",
    "cnv",
]


def clean_value(value):
    if pd.isna(value):
        return ""
    value = str(value).strip()
    if value.lower() in MISSING_VALUES:
        return ""
    return value


def output_value(value):
    """CSV에 다시 쓸 값은 Not specified를 빈칸으로 바꾸지 않고 보존한다."""
    if pd.isna(value):
        return ""
    return str(value).strip()


def norm_text(value):
    value = clean_value(value).lower()
    value = re.sub(r"\s+", " ", value)
    return value.strip()


def norm_variant(value):
    value = norm_text(value)
    if not value:
        return ""
    value = re.sub(r"nm_\d+(?:\.\d+)?:", "", value)
    value = value.replace("(", "").replace(")", "")
    value = re.sub(r"\s+", "", value)
    return value


def norm_sex(value):
    value = norm_text(value)
    if value in {"m", "male", "boy"}:
        return "male"
    if value in {"f", "female", "girl"}:
        return "female"
    return ""


def norm_role(value):
    value = norm_text(value)
    if value in {
        "original_patient",
        "cited_patient",
        "aggregate_only",
        "unclear",
        "variant_only",
        "manual_exclude",
        "reanalyzed_previous_patient",
        "reanalysed_previous_patient",
        "previously_reported_patient",
    }:
        return value
    if "original" in value or "own patient" in value or "own cohort" in value:
        return "original_patient"
    if "cited" in value or "literature" in value or "review" in value:
        return "cited_patient"
    if "aggregate" in value or "summary" in value:
        return "aggregate_only"
    if "variant_only" in value or "variant only" in value:
        return "variant_only"
    if "manual_exclude" in value or "manual exclude" in value:
        return "manual_exclude"
    if "reanalys" in value or "re-analy" in value or "previously reported patient" in value:
        return "reanalyzed_previous_patient"
    return "unclear"


def norm_onset(value):
    value = norm_text(value)
    if not value:
        return ""

    match = re.search(r"(\d+(?:\.\d+)?)\s*(day|days|d)\b", value)
    if match:
        return f"d:{float(match.group(1)):g}"

    match = re.search(r"(\d+(?:\.\d+)?)\s*(month|months|mo|mos|m)\b", value)
    if match:
        return f"m:{float(match.group(1)):g}"

    match = re.search(r"(\d+(?:\.\d+)?)\s*(year|years|yr|yrs|y)\b", value)
    if match:
        return f"m:{float(match.group(1)) * 12:g}"

    return value


def confidence_score(value):
    value = norm_text(value)
    if value == "high":
        return 3
    if value == "medium":
        return 2
    if value == "low":
        return 1
    return 0


BENIGN_PAPER_CLASSES = {"benign", "likely benign", "polymorphism", "non-pathogenic", "nonpathogenic"}


def variant_classification_flag(group):
    """VUS 환자 표시: 논문 분류 또는 ClinVar 분류가 VUS/상충이면 표시한다 (포함 여부는 바꾸지 않음)."""
    flags = []
    paper = {clean_value(v).lower() for v in group.get("variant_classification_in_paper", pd.Series(dtype=str))}
    if any("vus" in v or "uncertain" in v for v in paper):
        flags.append("VUS (paper)")
    clinvar = " ".join(clean_value(v).lower() for v in group.get("clinvar_classification", pd.Series(dtype=str)))
    if "uncertain" in clinvar:
        flags.append("VUS (ClinVar)")
    if "conflicting" in clinvar:
        flags.append("Conflicting (ClinVar)")
    return "; ".join(flags)


def first_existing(row, columns):
    for col in columns:
        if col in row:
            value = clean_value(row[col])
            if value:
                return value
    return ""


def infer_gene(input_csv, explicit_gene=None):
    if explicit_gene:
        return explicit_gene.upper()
    name = input_csv.name.upper()
    for gene in GENES:
        if gene in name:
            return gene
    return ""


def make_variant_key(row):
    c_key = norm_variant(first_existing(row, ["variant_c", "variant_c_gnomad", "variant_c_original"]))
    p_key = norm_variant(first_existing(row, ["variant_p"]))
    if c_key and p_key:
        return f"c={c_key}|p={p_key}"
    if c_key:
        return f"c={c_key}"
    if p_key:
        return f"p={p_key}"
    return ""


def make_dedup_key(row):
    gene = norm_text(row.get("gene", ""))
    variant_key = row.get("_variant_key", "")
    sex = row.get("_sex_norm", "")
    onset = row.get("_onset_norm", "")

    if not gene or not variant_key:
        return ""

    # Conservative automatic merge: same gene, same variant, same sex, same onset.
    # If sex or onset is missing, keep as singleton to avoid merging unrelated patients.
    if sex and onset:
        return f"gene={gene}|{variant_key}|sex={sex}|onset={onset}"
    return ""


def choose_primary(group):
    scored = group.copy()
    scored["_score"] = (
        scored["confidence"].map(confidence_score)
        + scored.apply(
            lambda row: 2 if first_existing(row, ["variant_c", "variant_c_gnomad", "variant_c_original"]) else 0,
            axis=1,
        )
        + scored["variant_p"].map(lambda x: 1 if clean_value(x) else 0)
        + scored["phenotype_review_category"].map(lambda x: 1 if clean_value(x) else 0)
        + scored["seizure_types"].map(lambda x: 1 if clean_value(x) else 0)
        + scored["other_phenotypes"].map(lambda x: 1 if clean_value(x) else 0)
    )
    return scored.sort_values(["_score", "pmid", "patient_id"], ascending=[False, True, True]).iloc[0]


def compact_values(values):
    cleaned = sorted({clean_value(value) for value in values if clean_value(value)})
    return "; ".join(cleaned)


def variant_group_key(row):
    key = clean_value(row.get("_variant_key", ""))
    if key:
        return key
    c_value = norm_variant(first_existing(row, ["variant_c", "variant_c_gnomad", "variant_c_original"]))
    p_value = norm_variant(first_existing(row, ["variant_p"]))
    patient_id = norm_text(row.get("patient_id", ""))
    return f"missing_variant|patient={patient_id}|c={c_value}|p={p_value}|mention={clean_value(row.get('mention_id', ''))}"


def same_variant_components(cluster):
    """병합된 한 환자 안에서 표기만 다른 같은 변이를 한 그룹으로 묶는다.

    논문마다 같은 변이를 다르게 쓰는 경우(cDNA 누락, "1502delG" vs "c.1502delG",
    "p.Asn1128fs" vs "p.Asn1128LysfsTer18")가 있어, cDNA가 같거나 단백질 핵심
    표기(frameshift는 위치까지)가 같으면 같은 변이로 본다.
    """
    def c_core(row):
        value = norm_variant(first_existing(row, ["variant_c", "variant_c_gnomad", "variant_c_original"]))
        if re.match(r"^c?\d", value):
            value = "c." + value.lstrip("c")
        value = re.sub(r"(del|dup)[acgt]+$", r"\1", value)  # c.602+1delG == c.602+1del
        if value.startswith("c.") and re.search(r"\d", value):
            return value
        # 서술형 exon CNV: "Exon 1–20 deletion" == "multi exon deletion (1-20)"
        exon = re.search(r"(\d+)\s*[-–~]\s*(\d+)", value) if "exon" in value else None
        kind = "del" if "del" in value else "dup" if "dup" in value else ""
        return f"exon{exon.group(1)}-{exon.group(2)}{kind}" if exon and kind else ""

    def p_core(row):
        value = norm_variant(first_existing(row, ["variant_p"]))
        match = re.match(r"^(p\.[a-z]{3}\d+)[a-z]*fs", value)
        if match:
            return match.group(1) + "fs"
        return value if value.startswith("p.") and re.search(r"\d", value) else ""

    labels = cluster.apply(variant_group_key, axis=1)
    parent = {idx: idx for idx in cluster.index}

    def find(idx):
        while parent[idx] != idx:
            parent[idx] = parent[parent[idx]]
            idx = parent[idx]
        return idx

    seen = {}
    for idx, row in cluster.iterrows():
        for key in (("c", c_core(row)), ("p", p_core(row))):
            if not key[1]:
                continue
            if key in seen:
                parent[find(idx)] = find(seen[key])
            else:
                seen[key] = idx
    roots = {}
    for idx in cluster.index:
        roots.setdefault(find(idx), []).append(labels[idx])
    return cluster.index.map(lambda idx: min(roots[find(idx)]))


def parse_mention_ids(value):
    """Parse manual merge mention_id lists from comma/semicolon/space separated text."""
    return [
        clean_value(part)
        for part in re.split(r"[;,\s]+", clean_value(value))
        if clean_value(part)
    ]


WORKSPACE_ROOT = Path(__file__).resolve().parents[2]
CONFIG_DIR = WORKSPACE_ROOT / "config"


def config_path_from_output(output_dir, filename):
    """config는 출력 위치와 관계없이 항상 작업공간 루트의 config/에서 읽는다."""
    return CONFIG_DIR / filename


def load_manual_merge_rules(output_dir, gene):
    candidate_paths = [config_path_from_output(output_dir, "manual_merge_rules.csv")]
    rules_path = next((path for path in candidate_paths if path.exists()), candidate_paths[0])
    if not rules_path.exists():
        return [], rules_path

    rules_df = pd.read_csv(rules_path, dtype=str).fillna("")
    required_cols = {"gene", "manual_cluster_id", "mention_ids"}
    missing_cols = required_cols - set(rules_df.columns)
    if missing_cols:
        raise ValueError(
            f"manual merge rules file is missing required columns: {sorted(missing_cols)}"
        )

    active_statuses = {"", "active", "apply", "applied"}
    loaded = []
    for _, row in rules_df.iterrows():
        if clean_value(row.get("gene", "")).upper() != gene:
            continue
        status = norm_text(row.get("status", "active"))
        if status not in active_statuses:
            continue
        mention_ids = parse_mention_ids(row.get("mention_ids", ""))
        if len(mention_ids) < 2:
            continue
        loaded.append({
            "manual_cluster_id": clean_value(row.get("manual_cluster_id", "")),
            "mention_ids": mention_ids,
            # 병합된 환자의 canonical 행에 보여줄 mention (비우면 choose_primary 점수로 선택)
            "preferred_mention_id": clean_value(row.get("preferred_mention_id", "")),
            "reason": clean_value(row.get("reason", "")),
            "evidence": clean_value(row.get("evidence", "")),
        })
    return loaded, rules_path


def load_manual_variant_exclusion_rules(output_dir, gene):
    """Load mention-level rules for non-causal or otherwise excluded variants."""
    candidate_paths = [config_path_from_output(output_dir, "manual_variant_exclusions.csv")]
    rules_path = next((path for path in candidate_paths if path.exists()), candidate_paths[0])
    if not rules_path.exists():
        return [], rules_path

    rules_df = pd.read_csv(rules_path, dtype=str).fillna("")
    required_cols = {"gene", "rule_id", "reason"}
    missing_cols = required_cols - set(rules_df.columns)
    if missing_cols:
        raise ValueError(
            f"manual variant exclusion file is missing required columns: {sorted(missing_cols)}"
        )

    active_statuses = {"", "active", "apply", "applied"}
    loaded = []
    for _, row in rules_df.iterrows():
        if clean_value(row.get("gene", "")).upper() != gene:
            continue
        status = norm_text(row.get("status", "active"))
        if status not in active_statuses:
            continue
        loaded.append({
            "rule_id": clean_value(row.get("rule_id", "")),
            "mention_id": clean_value(row.get("mention_id", "")),
            "pmid": clean_value(row.get("pmid", "")),
            "patient_id": clean_value(row.get("patient_id", "")),
            "variant_c": clean_value(row.get("variant_c", "")),
            "variant_p": clean_value(row.get("variant_p", "")),
            "reason": clean_value(row.get("reason", "")) or "manual_variant_exclusion",
            "evidence": clean_value(row.get("evidence", "")),
        })
    return loaded, rules_path


def manual_variant_exclusion_for_row(row, rules):
    mention_id = clean_value(row.get("mention_id", ""))
    pmid = clean_value(row.get("pmid", ""))
    patient_id = norm_text(row.get("patient_id", ""))
    variant_c = norm_variant(first_existing(row, ["variant_c", "variant_c_gnomad", "variant_c_original"]))
    variant_p = norm_variant(row.get("variant_p", ""))

    for rule in rules:
        if rule.get("mention_id"):
            if mention_id == rule["mention_id"]:
                return rule
            continue

        if rule.get("pmid") and pmid != rule["pmid"]:
            continue
        if rule.get("patient_id") and patient_id != norm_text(rule["patient_id"]):
            continue
        if rule.get("variant_c") and variant_c != norm_variant(rule["variant_c"]):
            continue
        if rule.get("variant_p") and variant_p != norm_variant(rule["variant_p"]):
            continue
        if any(rule.get(key) for key in ["pmid", "patient_id", "variant_c", "variant_p"]):
            return rule
    return None


def canonical_id_number(value):
    match = re.fullmatch(r"CP(\d+)", clean_value(value))
    if not match:
        return 0
    return int(match.group(1))


def format_canonical_id(number):
    return f"CP{number:07d}"


def load_reserved_canonical_ids(gene):
    """교수님 CRF 통합본 등 외부에서 이미 쓴 CP 번호. 새 환자에게 다시 주지 않는다."""
    path = CONFIG_DIR / "reserved_canonical_ids" / f"{gene}.csv"
    if not path.exists():
        return set()
    df = pd.read_csv(path, dtype=str).fillna("")
    return {clean_value(v) for v in df.get("canonical_patient_id", []) if canonical_id_number(clean_value(v))}


def load_canonical_id_locks(output_dir, gene):
    """Load optional mention_id -> canonical_patient_id locks for stable CP IDs."""
    # 실행 폴더에 lock 파일을 따로 두면 그것을 우선 사용하고, 없으면 config의 현재 lock을 쓴다.
    candidate_paths = [
        output_dir / f"{gene}_canonical_id_locks.csv",
        CONFIG_DIR / "canonical_id_locks" / f"{gene}.csv",
    ]
    lock_path = next((path for path in candidate_paths if path.exists()), candidate_paths[0])
    if not lock_path.exists():
        return {}, lock_path

    locks_df = pd.read_csv(lock_path, dtype=str).fillna("")
    required_cols = {"mention_id", "canonical_patient_id"}
    missing_cols = required_cols - set(locks_df.columns)
    if missing_cols:
        raise ValueError(
            f"canonical ID lock file is missing required columns: {sorted(missing_cols)}"
        )

    locks = {}
    for _, row in locks_df.iterrows():
        row_gene = clean_value(row.get("gene", "")).upper()
        if row_gene and row_gene != gene:
            continue
        mention_id = clean_value(row.get("mention_id", ""))
        canonical_id = clean_value(row.get("canonical_patient_id", ""))
        if not mention_id or not canonical_id:
            continue
        if not canonical_id_number(canonical_id):
            raise ValueError(f"invalid locked canonical_patient_id: {canonical_id}")
        existing = locks.get(mention_id)
        if existing and existing != canonical_id:
            raise ValueError(
                f"conflicting canonical ID locks for mention_id {mention_id}: "
                f"{existing} vs {canonical_id}"
            )
        locks[mention_id] = canonical_id
    return locks, lock_path


def has_patient_id_conflict_within_paper(group):
    """같은 논문 안에서 서로 다른 patient_id가 있으면 자동 병합하지 않는다."""
    for _, paper_group in group.groupby("pmid", dropna=False):
        patient_ids = {clean_value(value) for value in paper_group["patient_id"] if clean_value(value)}
        if len(patient_ids) > 1:
            return True
    return False


def has_newly_reported_patient_evidence(row):
    text = norm_text(
        " ".join(
            [
                str(row.get("originality_evidence", "")),
                str(row.get("location_in_paper", "")),
                str(row.get("cited_source", "")),
            ]
        )
    )
    return any(phrase in text for phrase in NEWLY_REPORTED_PATIENT_PHRASES)


def has_previously_reported_patient_evidence(row):
    """Experimental guard: detect rows that look like prior/reanalysed patients."""
    if has_newly_reported_patient_evidence(row):
        return False

    evidence_text = norm_text(
        " ".join(
            [
                str(row.get("originality_evidence", "")),
                str(row.get("location_in_paper", "")),
            ]
        )
    )

    if any(phrase in evidence_text for phrase in NEGATED_PREVIOUSLY_REPORTED_PHRASES):
        return False

    # Do not confuse a previously reported variant with a previously reported patient.
    evidence_text = re.sub(
        r"[^.]*\bpreviously\s+(?:described|reported|published)\s+variants?\b[^.]*",
        "",
        evidence_text,
    )
    evidence_text = re.sub(
        r"[^.]*\bvariant\b[^.]*\bpreviously\b[^.]*",
        "",
        evidence_text,
    )

    patient_terms = r"(patient|patients|individual|individuals|participant|participants|cohort|case|cases)"
    prior_terms = r"(previously described|previously reported|previously published|reported by|re-analyzed|re-analysed|reanalyzed|reanalysed)"
    if re.search(patient_terms + r".{0,100}" + prior_terms, evidence_text):
        return True
    if re.search(prior_terms + r".{0,100}" + patient_terms, evidence_text):
        return True

    return False


def other_phenotypes_is_generic_diagnosis(value):
    text = norm_text(value)
    if not text:
        return False
    if any(phrase in text for phrase in GENERIC_DIAGNOSIS_PHRASES):
        remaining = text
        for phrase in GENERIC_DIAGNOSIS_PHRASES:
            remaining = remaining.replace(phrase, "")
        for phrase in GENERIC_VARIANT_ONLY_PHRASES:
            remaining = remaining.replace(phrase, "")
        remaining = re.sub(r"[;,\s-]+", "", remaining)
        return not remaining
    return False


def has_pmid_level_prior_cohort_signal(group):
    """Experimental guard for mixed cohorts with prior + newly reported patients."""
    has_new_statement = group.apply(has_newly_reported_patient_evidence, axis=1).any()
    if not has_new_statement:
        return False

    text = norm_text(
        " ".join(
            str(value)
            for value in group.get("originality_evidence", pd.Series(dtype=str)).fillna("")
        )
    )
    count_prior_terms = (
        r"(previously described|previously reported|previously published|"
        r"re-analyzed|re-analysed|reanalyzed|reanalysed)"
    )
    count_prior_pattern = (
        r"\b\d+\s+(?:of\s+\d+\s+)?"
        r"(?:patients?|individuals?|participants?|cases?)?\s*"
        r"(?:had been|were|was|are|been|had)?\s*"
        + count_prior_terms
    )
    return bool(re.search(count_prior_pattern, text))


def has_patient_level_clinical_detail(row, strict_source_filter=False):
    """환자별 임상정보 없이 variant list만 있는 row는 canonical 후보에서 제외한다."""
    direct_clinical_cols = [
        "seizure_presence",
        "age_of_onset",
        "last_followup_age",
        "seizure_types",
        "phenotype_review_category",
        "response_to_asm",
        "asm_used",
        "responsive_asm",
        "targeted_treatment_response",
        "epilepsy_type",
    ]
    if not strict_source_filter:
        direct_clinical_cols = ["sex"] + direct_clinical_cols

    for col in direct_clinical_cols:
        if clean_value(row.get(col, "")):
            return True

    other = norm_text(row.get("other_phenotypes", ""))
    if not other:
        return False

    if strict_source_filter and other_phenotypes_is_generic_diagnosis(other):
        return False

    non_patient_phrases = [
        "patient-specific phenotypes not provided",
        "no patient-specific",
        "cohort-level phenotype",
        "patient-derived ipsc",
        "patient-derived neuronal",
        "patient-derived cells",
        "ipsc-neurons",
        "haploinsufficiency",
    ]
    if any(phrase in other for phrase in non_patient_phrases):
        return False

    return True


def add_sort_columns(df):
    sorted_df = df.copy()
    role_order = {
        "original_patient": 1,
        "cited_patient": 2,
        "aggregate_only": 3,
        "variant_only": 4,
        "reanalyzed_previous_patient": 5,
        "reanalysed_previous_patient": 5,
        "previously_reported_patient": 5,
        "manual_exclude": 6,
        "unclear": 7,
    }
    variant_c_series = sorted_df.apply(
        lambda row: first_existing(row, ["variant_c", "variant_c_gnomad", "variant_c_original"]),
        axis=1,
    )
    sorted_df["_sort_variant_c_missing"] = variant_c_series.map(lambda value: 1 if not norm_variant(value) else 0)
    variant_p_series = (
        sorted_df["variant_p"]
        if "variant_p" in sorted_df.columns
        else pd.Series([""] * len(sorted_df), index=sorted_df.index)
    )
    sorted_df["_sort_variant_p_missing"] = variant_p_series.map(
        lambda value: 1 if not norm_variant(value) else 0
    )
    sorted_df["_sort_variant_c"] = variant_c_series.map(norm_variant)
    sorted_df["_sort_variant_p"] = variant_p_series.map(norm_variant)
    role_series = (
        sorted_df["mention_source_type"]
        if "mention_source_type" in sorted_df.columns
        else pd.Series([""] * len(sorted_df), index=sorted_df.index)
    )
    sorted_df["_sort_role"] = role_series.map(
        lambda value: role_order.get(norm_role(value), 9)
    )
    return sorted_df


def sort_for_manual_review(df, extra_keys=None):
    if df.empty:
        return df
    sorted_df = add_sort_columns(df)
    sort_keys = ["_sort_variant_c_missing", "_sort_variant_c", "_sort_variant_p_missing", "_sort_variant_p"]
    if "_sort_role" in sorted_df.columns:
        sort_keys.append("_sort_role")
    for key in extra_keys or []:
        if key in sorted_df.columns:
            sort_keys.append(key)
    sorted_df = sorted_df.sort_values(sort_keys, kind="mergesort")
    return sorted_df.drop(columns=[c for c in sorted_df.columns if c.startswith("_sort_")])


def build_outputs(
    input_csv,
    output_dir,
    gene=None,
    experimental_strict_source_filter=False,
    output_timestamp=None,
):
    df = pd.read_csv(input_csv, dtype=str).fillna("")
    df = df.copy()
    source_columns = list(df.columns)

    if "mention_id" not in df.columns:
        df["mention_id"] = [str(i) for i in range(1, len(df) + 1)]
        source_columns.append("mention_id")

    required_cols = [
        "variant_c", "variant_c_gnomad", "variant_p", "variant_type",
        "patient_id", "pmid", "sex", "korean_patient_source", "korean_patient_evidence",
        "age_of_onset", "last_followup_age", "seizure_types",
        "phenotype_review_category", "response_to_asm", "asm_used",
        "responsive_asm", "targeted_treatment_response", "epilepsy_type",
        "other_phenotypes", "mention_source_type", "originality_evidence",
        "cited_source", "location_in_paper", "confidence", "pmcid",
        "variant_c_original", "gnomad_conversion_status", "nm_transcript",
        "gnomad_transcript",
    ]
    for col in required_cols:
        if col not in df.columns:
            df[col] = ""

    inferred_gene = infer_gene(input_csv, gene)
    if "gene" not in df.columns or not df["gene"].map(clean_value).any():
        df["gene"] = inferred_gene

    df["_role_norm"] = df["mention_source_type"].map(norm_role)
    df["_sex_norm"] = df["sex"].map(norm_sex)
    df["_onset_norm"] = df["age_of_onset"].map(norm_onset)
    df["_variant_key"] = df.apply(make_variant_key, axis=1)
    df["_dedup_key"] = df.apply(make_dedup_key, axis=1)

    df["canonical_patient_id"] = ""
    df["canonical_status"] = "excluded"
    df["exclusion_reason"] = ""
    df["manual_exclusion_rule_id"] = ""
    df["manual_exclusion_evidence"] = ""

    exclusion_reasons = {
        "cited_patient": "cited_from_prior_literature",
        "aggregate_only": "aggregate_only_no_individual_patient",
        "variant_only": "variant_only_no_patient_level_clinical_information",
        "reanalyzed_previous_patient": "previously_reported_or_reanalyzed_patient",
        "reanalysed_previous_patient": "previously_reported_or_reanalyzed_patient",
        "previously_reported_patient": "previously_reported_or_reanalyzed_patient",
        "manual_exclude": "manual_excluded_after_source_review",
        "unclear": "unclear_originality",
    }
    for role, reason in exclusion_reasons.items():
        mask = df["_role_norm"] == role
        df.loc[mask, "exclusion_reason"] = reason

    prior_patient_mask = pd.Series(False, index=df.index)
    if experimental_strict_source_filter:
        prior_patient_mask = (
            (df["_role_norm"] == "original_patient")
            & df.apply(has_previously_reported_patient_evidence, axis=1)
        )
        for _, group in df[df["_role_norm"] == "original_patient"].groupby("pmid", dropna=False):
            if not has_pmid_level_prior_cohort_signal(group):
                continue
            newly_reported = group.apply(has_newly_reported_patient_evidence, axis=1)
            prior_patient_mask.loc[group.index] = prior_patient_mask.loc[group.index] | ~newly_reported
        df.loc[prior_patient_mask, "exclusion_reason"] = "previously_reported_or_reanalyzed_patient"

    df["_has_patient_clinical_detail"] = df.apply(
        lambda row: has_patient_level_clinical_detail(
            row,
            strict_source_filter=experimental_strict_source_filter,
        ),
        axis=1,
    )
    variant_only_mask = (
        (df["_role_norm"] == "original_patient")
        & ~prior_patient_mask
        & ~df["_has_patient_clinical_detail"]
    )
    df.loc[variant_only_mask, "exclusion_reason"] = "variant_only_no_patient_level_clinical_information"

    # cDNA·protein이 모두 없는 행은 조용히 버리지 않고 excluded에 사유를 남긴다 (원문 재확인 대상).
    no_variant_mask = df.apply(
        lambda row: not first_existing(row, ["variant_c", "variant_c_gnomad", "variant_c_original", "variant_p"]),
        axis=1,
    )
    df.loc[
        no_variant_mask & (df["exclusion_reason"] == "") & (df["_role_norm"] == "original_patient"),
        "exclusion_reason",
    ] = "variant_not_recorded"

    # 논문이 benign/다형성으로 분류한 변이 행은 원인 변이에서 제외한다 (2026-09-23 변이 분류 기준).
    benign_mask = df.get("variant_classification_in_paper", pd.Series("", index=df.index)).map(
        lambda value: clean_value(value).lower() in BENIGN_PAPER_CLASSES
    )
    df.loc[
        benign_mask & (df["exclusion_reason"] == "") & (df["_role_norm"] == "original_patient"),
        "exclusion_reason",
    ] = "benign_or_polymorphism_in_paper"

    manual_variant_rules, manual_variant_rules_path = load_manual_variant_exclusion_rules(
        output_dir,
        inferred_gene,
    )
    manual_variant_rule_results = []
    manual_variant_exclusion_mask = pd.Series(False, index=df.index)
    if manual_variant_rules:
        for idx, row in df.iterrows():
            rule = manual_variant_exclusion_for_row(row, manual_variant_rules)
            if not rule:
                continue
            manual_variant_exclusion_mask.loc[idx] = True
            df.at[idx, "exclusion_reason"] = rule["reason"]
            df.at[idx, "manual_exclusion_rule_id"] = rule["rule_id"]
            df.at[idx, "manual_exclusion_evidence"] = rule["evidence"]
            manual_variant_rule_results.append({
                "rule_id": rule["rule_id"],
                "mention_id": clean_value(row.get("mention_id", "")),
                "pmid": clean_value(row.get("pmid", "")),
                "patient_id": clean_value(row.get("patient_id", "")),
                "variant_c": clean_value(first_existing(row, ["variant_c", "variant_c_gnomad", "variant_c_original"])),
                "variant_p": clean_value(row.get("variant_p", "")),
                "reason": rule["reason"],
                "status": "applied",
            })

    original = df[
        (df["_role_norm"] == "original_patient")
        & ~prior_patient_mask
        & ~variant_only_mask
        & ~manual_variant_exclusion_mask
        & ~no_variant_mask
        & ~benign_mask
    ].copy()

    canonical_records = []
    canonical_patient_records = []
    link_records = []
    next_id = 1
    key_to_canonical = {}
    manual_rules, manual_rules_path = load_manual_merge_rules(output_dir, inferred_gene)
    manual_rule_results = []
    canonical_id_locks, canonical_id_locks_path = load_canonical_id_locks(output_dir, inferred_gene)
    used_canonical_ids = set()
    canonical_id_source_counts = {"locked": 0, "new": 0}
    reserved_ids = load_reserved_canonical_ids(inferred_gene)
    known_ids = set(canonical_id_locks.values()) | reserved_ids
    if known_ids:
        # 새 번호는 lock·예약 번호(통합본 CP 포함)보다 항상 크게: 다른 환자가 쓰던 번호를 다시 주지 않음
        next_id = max(canonical_id_number(value) for value in known_ids) + 1

    def allocate_canonical_id(cluster):
        nonlocal next_id
        locked_ids = sorted({
            canonical_id_locks[clean_value(mention_id)]
            for mention_id in cluster["mention_id"].tolist()
            if clean_value(mention_id) in canonical_id_locks
        })
        if len(locked_ids) > 1:
            raise ValueError(
                "canonical ID lock conflict: one new cluster contains mentions "
                f"locked to multiple IDs: {locked_ids}"
            )
        if locked_ids:
            canonical_id = locked_ids[0]
            used_canonical_ids.add(canonical_id)
            return canonical_id, "locked"

        while True:
            canonical_id = format_canonical_id(next_id)
            next_id += 1
            if canonical_id not in used_canonical_ids and canonical_id not in known_ids:
                used_canonical_ids.add(canonical_id)
                return canonical_id, "new"

    preferred_display_ids = {rule["preferred_mention_id"] for rule in manual_rules if rule.get("preferred_mention_id")}

    def pick_primary(group):
        preferred = group[group["mention_id"].map(lambda value: clean_value(value) in preferred_display_ids)]
        return preferred.iloc[0] if len(preferred) else choose_primary(group)

    def add_cluster(cluster, match_confidence, match_reason):
        canonical_id, canonical_id_source = allocate_canonical_id(cluster)
        canonical_id_source_counts[canonical_id_source] += 1

        primary = pick_primary(cluster)
        mention_ids = cluster["mention_id"].tolist()
        pmids = compact_values(cluster["pmid"])
        pmcids = compact_values(cluster["pmcid"])
        patient_ids = compact_values(cluster["patient_id"])
        variant_groups = [
            variant_group
            for _, variant_group in cluster.groupby(same_variant_components(cluster), sort=True)
        ]

        df.loc[df["mention_id"].isin(mention_ids), "canonical_patient_id"] = canonical_id
        df.loc[df["mention_id"].isin(mention_ids), "canonical_status"] = "included"
        df.loc[df["mention_id"].isin(mention_ids), "exclusion_reason"] = ""

        dedup_key = clean_value(primary.get("_dedup_key", ""))
        for _, row in cluster.iterrows():
            row_dedup_key = clean_value(row.get("_dedup_key", ""))
            if not row_dedup_key:
                continue
            canonical_ids = key_to_canonical.setdefault(row_dedup_key, [])
            if canonical_id not in canonical_ids:
                canonical_ids.append(canonical_id)

        canonical_patient_records.append({
            "canonical_patient_id": canonical_id,
            "gene": clean_value(primary.get("gene", "")),
            "supporting_pmids": pmids,
            "supporting_pmcids": pmcids,
            "source_patient_ids": patient_ids,
            "mention_count": len(cluster),
            "variant_row_count": len(variant_groups),
            "dedup_key": dedup_key,
            "match_confidence": match_confidence,
            "match_reason": match_reason,
        })

        for variant_group in variant_groups:
            variant_primary = pick_primary(variant_group)
            variant_pmids = compact_values(variant_group["pmid"])
            variant_pmcids = compact_values(variant_group["pmcid"])
            variant_patient_ids = compact_values(variant_group["patient_id"])
            variant_mention_ids = compact_values(variant_group["mention_id"])
            variant_source_values = {
                col: output_value(variant_primary.get(col, ""))
                for col in source_columns
            }
            record = {
                "canonical_patient_id": canonical_id,
                "gene": clean_value(variant_primary.get("gene", "")),
                **variant_source_values,
                "supporting_pmids": variant_pmids,
                "supporting_pmcids": variant_pmcids,
                "source_patient_ids": variant_patient_ids,
                "supporting_mention_ids": variant_mention_ids,
                "mention_count": len(variant_group),
                "canonical_patient_mention_count": len(cluster),
                "canonical_patient_variant_count": len(variant_groups),
                "dedup_key": clean_value(variant_primary.get("_dedup_key", "")),
                "canonical_patient_dedup_keys": compact_values(cluster["_dedup_key"]),
                "match_confidence": match_confidence,
                "match_reason": match_reason,
            }
            if {"variant_classification_in_paper", "clinvar_classification"} & set(source_columns):
                record["variant_classification_flag"] = variant_classification_flag(variant_group)
            canonical_records.append(record)

        for _, row in cluster.iterrows():
            link_records.append({
                "canonical_patient_id": canonical_id,
                "mention_id": clean_value(row["mention_id"]),
                "pmid": clean_value(row.get("pmid", "")),
                "pmcid": clean_value(row.get("pmcid", "")),
                "patient_id": clean_value(row.get("patient_id", "")),
                "mention_source_type": clean_value(row.get("mention_source_type", "")),
                "match_confidence": match_confidence,
                "match_reason": match_reason,
            })

    clustered_mention_ids = set()
    original_by_mention_id = {
        clean_value(row["mention_id"]): idx
        for idx, row in original.iterrows()
    }
    for rule in manual_rules:
        requested_ids = [mid for mid in rule["mention_ids"] if mid in original_by_mention_id]
        already_clustered = [mid for mid in requested_ids if mid in clustered_mention_ids]
        available_ids = [mid for mid in requested_ids if mid not in clustered_mention_ids]
        cluster = original[original["mention_id"].map(lambda value: clean_value(value) in available_ids)].copy()
        missing_ids = [mid for mid in rule["mention_ids"] if mid not in original_by_mention_id]
        result = {
            "manual_cluster_id": rule["manual_cluster_id"],
            "requested_mention_ids": rule["mention_ids"],
            "applied_mention_ids": [],
            "missing_or_ineligible_mention_ids": missing_ids,
            "already_clustered_mention_ids": already_clustered,
            "status": "skipped",
        }
        if len(cluster) >= 2:
            cluster_id = rule["manual_cluster_id"] or "unnamed"
            add_cluster(cluster, "manual", f"manual_merge_rule:{cluster_id}")
            clustered_mention_ids.update(clean_value(value) for value in cluster["mention_id"])
            result["applied_mention_ids"] = [clean_value(value) for value in cluster["mention_id"]]
            result["status"] = "applied"
        elif len(cluster) == 1:
            result["applied_mention_ids"] = [clean_value(cluster.iloc[0]["mention_id"])]
            result["status"] = "skipped_only_one_eligible_mention"
        manual_rule_results.append(result)

    unclustered_original = original[
        ~original["mention_id"].map(lambda value: clean_value(value) in clustered_mention_ids)
    ].copy()
    for _, group in unclustered_original.groupby(["pmid", "patient_id"], sort=True, dropna=False):
        pmid = clean_value(group.iloc[0].get("pmid", ""))
        patient_id = clean_value(group.iloc[0].get("patient_id", ""))
        if not pmid or not patient_id or len(group) <= 1:
            continue
        add_cluster(group, "high", "same_paper_same_patient_id")
        clustered_mention_ids.update(clean_value(value) for value in group["mention_id"])

    remaining_original = original[
        ~original["mention_id"].map(lambda value: clean_value(value) in clustered_mention_ids)
    ].copy()
    original_keyed = remaining_original[remaining_original["_dedup_key"] != ""].copy()
    original_no_key = remaining_original[remaining_original["_dedup_key"] == ""].copy()

    for _, group in original_keyed.groupby("_dedup_key", sort=True):
        if len(group) > 1:
            for _, row in group.iterrows():
                same_paper_conflict = has_patient_id_conflict_within_paper(group)
                add_cluster(
                    pd.DataFrame([row]),
                    "singleton",
                    "conflicting_patient_ids_in_same_paper"
                    if same_paper_conflict
                    else "same_gene_variant_sex_onset_not_auto_merged",
                )
        else:
            add_cluster(group, "singleton", "single_original_patient_with_strong_key")

    for _, row in original_no_key.iterrows():
        add_cluster(pd.DataFrame([row]), "singleton", "single_original_patient_without_merge_key")

    # Cited rows are not included in canonical_patients, but link them as supporting
    # evidence when they match an existing canonical dedup key.
    cited = df[df["_role_norm"] == "cited_patient"].copy()
    for _, row in cited.iterrows():
        canonical_ids = key_to_canonical.get(clean_value(row.get("_dedup_key", "")), [])
        if len(canonical_ids) != 1:
            continue
        canonical_id = canonical_ids[0]
        link_records.append({
            "canonical_patient_id": canonical_id,
            "mention_id": clean_value(row["mention_id"]),
            "pmid": clean_value(row.get("pmid", "")),
            "pmcid": clean_value(row.get("pmcid", "")),
            "patient_id": clean_value(row.get("patient_id", "")),
            "mention_source_type": clean_value(row.get("mention_source_type", "")),
            "match_confidence": "supporting",
            "match_reason": "cited_patient_same_dedup_key",
        })

    output_dir.mkdir(parents=True, exist_ok=True)

    internal_cols = [
        "_role_norm", "_sex_norm", "_onset_norm", "_variant_key", "_dedup_key",
        "_has_patient_clinical_detail",
    ]
    extra_mention_cols = [
        col for col in [
            "gene",
            "canonical_patient_id",
            "canonical_status",
            "exclusion_reason",
            "manual_exclusion_rule_id",
            "manual_exclusion_evidence",
        ]
        if col in df.columns and col not in source_columns
    ]
    mentions = df[source_columns + extra_mention_cols].copy()
    canonical = pd.DataFrame(canonical_records)
    links = pd.DataFrame(link_records)
    excluded = mentions[mentions["canonical_status"] == "excluded"].copy()

    mentions = sort_for_manual_review(mentions, ["pmid", "patient_id", "mention_id"])
    canonical = sort_for_manual_review(canonical, ["supporting_pmids", "patient_id", "canonical_patient_id"])
    excluded = sort_for_manual_review(excluded, ["pmid", "patient_id", "mention_id"])
    if not links.empty:
        link_sort_cols = [col for col in ["canonical_patient_id", "pmid", "patient_id", "mention_id"] if col in links.columns]
        links = links.sort_values(link_sort_cols, kind="mergesort")

    file_prefix = inferred_gene or input_csv.stem
    output_files = {
        "patient_mentions": timestamp_path(
            output_dir / f"{file_prefix}_patient_mentions.csv",
            output_timestamp,
        ),
        "canonical_patients": timestamp_path(
            output_dir / f"{file_prefix}_canonical_patients.csv",
            output_timestamp,
        ),
        "patient_dedup_links": timestamp_path(
            output_dir / f"{file_prefix}_patient_dedup_links.csv",
            output_timestamp,
        ),
        "excluded_mentions": timestamp_path(
            output_dir / f"{file_prefix}_excluded_mentions.csv",
            output_timestamp,
        ),
    }
    mentions.to_csv(output_files["patient_mentions"], index=False, encoding="utf-8-sig")
    canonical.to_csv(output_files["canonical_patients"], index=False, encoding="utf-8-sig")
    links.to_csv(output_files["patient_dedup_links"], index=False, encoding="utf-8-sig")
    excluded.to_csv(output_files["excluded_mentions"], index=False, encoding="utf-8-sig")

    summary = {
        "input_csv": str(input_csv),
        "gene": inferred_gene,
        "total_mentions": int(len(df)),
        "canonical_patients": int(len(canonical_patient_records)),
        "canonical_patient_variant_rows": int(len(canonical)),
        "included_original_mentions": int((mentions["canonical_status"] == "included").sum()),
        "excluded_mentions": int(len(excluded)),
        "linked_cited_mentions": int((links.get("match_reason", pd.Series(dtype=str)) == "cited_patient_same_dedup_key").sum()),
        "excluded_by_reason": excluded["exclusion_reason"].value_counts().to_dict(),
        "merged_clusters": int(sum(1 for item in canonical_patient_records if item["mention_count"] > 1)),
        "manual_merge_rule_file": str(manual_rules_path) if manual_rules_path.exists() else "",
        "manual_merge_rules_loaded": int(len(manual_rules)),
        "manual_merge_rules_applied": int(sum(1 for item in manual_rule_results if item["status"] == "applied")),
        "manual_merge_rule_results": manual_rule_results,
        "manual_variant_exclusion_file": (
            str(manual_variant_rules_path) if manual_variant_rules_path.exists() else ""
        ),
        "manual_variant_exclusion_rules_loaded": int(len(manual_variant_rules)),
        "manual_variant_exclusions_applied": int(len(manual_variant_rule_results)),
        "manual_variant_exclusion_results": manual_variant_rule_results,
        "canonical_id_lock_file": str(canonical_id_locks_path) if canonical_id_locks_path.exists() else "",
        "canonical_id_locks_loaded": int(len(canonical_id_locks)),
        "canonical_id_locked_clusters": int(canonical_id_source_counts.get("locked", 0)),
        "canonical_id_new_clusters": int(canonical_id_source_counts.get("new", 0)),
        "csv_timestamp": output_timestamp or "",
        "output_files": {key: str(path) for key, path in output_files.items()},
        "policy": {
            "canonical_roles": ["original_patient"],
            "excluded_roles": ["cited_patient", "aggregate_only", "unclear"],
            "auto_merge_rule": "no cross-paper original-original auto merge by variant/sex/onset",
            "within_paper_patient_merge": "same PMID + same non-empty patient_id",
            "canonical_output_grain": "one row per canonical patient variant; same canonical_patient_id can appear on multiple rows when one patient has multiple variants",
            "auto_merge_guard": "same variant/sex/onset is treated as a review signal, not an automatic original-patient merge",
            "manual_merge_rule": "explicit source-supported same-patient mappings from config/manual_merge_rules.csv",
            "cited_link_rule": "link cited_patient rows only when their dedup key maps to exactly one canonical patient",
            "missing sex or onset": "kept as singleton",
            "variant_only_exclusion": "original_patient rows without patient-level clinical detail are excluded",
            "manual_variant_exclusion": "source-reviewed non-causal or secondary variants are excluded from canonical causal output while preserved in patient_mentions",
            "experimental_strict_source_filter": (
                "enabled" if experimental_strict_source_filter else "disabled"
            ),
        },
    }
    with open(output_dir / "dedup_summary.json", "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2, ensure_ascii=False)

    return summary


def main():
    parser = argparse.ArgumentParser(description="Create canonical patients from patient_mentions CSV")
    parser.add_argument("input_csv", type=Path, help="patient_mentions CSV")
    parser.add_argument("--output-dir", type=Path, default=None, help="Output directory")
    parser.add_argument("--gene", default=None, help="Gene name if the CSV filename does not include it")
    parser.add_argument(
        "--timestamp",
        default=None,
        help="CSV filename timestamp in YYYYMMDD_HHMMSS format. Defaults to current time.",
    )
    parser.add_argument(
        "--no-timestamp",
        action="store_true",
        help="Write legacy CSV filenames without timestamp.",
    )
    parser.add_argument(
        "--experimental-strict-source-filter",
        action="store_true",
        help=(
            "Experimental: exclude original_patient rows that look previously "
            "reported/reanalysed, and apply stricter variant-only filtering."
        ),
    )
    args = parser.parse_args()

    output_dir = args.output_dir
    if output_dir is None:
        output_dir = args.input_csv.parent / f"{args.input_csv.stem}_canonical_patients"

    summary = build_outputs(
        args.input_csv,
        output_dir,
        gene=args.gene,
        experimental_strict_source_filter=args.experimental_strict_source_filter,
        output_timestamp=None if args.no_timestamp else make_timestamp(args.timestamp),
    )
    print("Deduplication complete")
    print(f"  total_mentions: {summary['total_mentions']}")
    print(f"  canonical_patients: {summary['canonical_patients']}")
    print(f"  canonical_patient_variant_rows: {summary['canonical_patient_variant_rows']}")
    print(f"  included_original_mentions: {summary['included_original_mentions']}")
    print(f"  excluded_mentions: {summary['excluded_mentions']}")
    print(f"  linked_cited_mentions: {summary['linked_cited_mentions']}")
    print(f"  merged_clusters: {summary['merged_clusters']}")
    print(f"  output_dir: {output_dir}")


if __name__ == "__main__":
    main()
