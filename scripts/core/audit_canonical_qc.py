#!/usr/bin/env python3
"""Build QC reports from current canonical patient CSV outputs."""

import argparse
import csv
import json
from collections import Counter, defaultdict
from itertools import combinations
from pathlib import Path

from timestamped_outputs import make_timestamp, resolve_csv_path, timestamp_path


GENES = [
    "CDKL5", "EEF1A2", "KCNQ2", "SCN1A", "SCN2A", "SYNGAP1",
    "PCDH19", "SCN8A", "GABRA1", "STXBP1", "SMC1A", "FOXG1", "KCNT1",
]
MISSING = {"", "not specified", "n/a", "na", "none", "unknown", "nan"}


def clean(value):
    if value is None:
        return ""
    value = str(value).strip()
    return "" if value.lower() in MISSING else value


def yes(value):
    return clean(value).lower() in {"yes", "y", "1", "true"}


def norm_text(value):
    value = clean(value).lower()
    return " ".join(value.split())


def norm_variant(value):
    value = norm_text(value)
    if not value:
        return ""
    value = value.replace("(", "").replace(")", "")
    return "".join(value.split())


def causal_variant_key(row):
    """Use protein consequence first to avoid flagging cDNA-missing duplicate rows."""
    variant_p = norm_variant(row.get("variant_p"))
    if variant_p and variant_p not in {"p.?", "p.="}:
        return f"p={variant_p}"
    variant_c = norm_variant(row.get("variant_c"))
    if variant_c:
        return f"c={variant_c}"
    return f"mention={clean(row.get('mention_id'))}"


def read_csv(path):
    with path.open(newline="", encoding="utf-8-sig") as f:
        return list(csv.DictReader(f))


def write_csv(path, rows, fieldnames, output_timestamp=None):
    path = timestamp_path(path, output_timestamp)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8-sig") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    return path


def compact(values, limit=8):
    kept = []
    seen = set()
    for value in values:
        value = clean(value)
        if not value or value in seen:
            continue
        seen.add(value)
        kept.append(value)
    suffix = ""
    if len(kept) > limit:
        suffix = f"; ...(+{len(kept) - limit})"
        kept = kept[:limit]
    return "; ".join(kept) + suffix


def default_auxiliary_root(canonical_root):
    """유전자별 보조 산출물(patient_mentions, excluded, dedup_links)은 outputs/auxiliary/{GENE}/에 있다."""
    candidate = (canonical_root.parent.parent / "auxiliary").resolve()
    return candidate if candidate.exists() else canonical_root

def find_gene_csv(canonical_root, auxiliary_root, gene, suffix):
    stem = f"{gene}_{suffix}"
    gene_dir = canonical_root / f"{gene}_canonical_patients"
    path = resolve_csv_path(gene_dir, stem)
    if path:
        return path
    if auxiliary_root and auxiliary_root != canonical_root:
        aux_dir = auxiliary_root / gene
        return resolve_csv_path(aux_dir, stem)
    return None


def load_inheritance_models(path):
    defaults = {
        gene: {
            "inheritance_model": "autosomal_dominant",
            "allowed_causal_variant_count": "1",
            "review_rule": "flag_if_more_than_allowed",
            "notes": "",
        }
        for gene in GENES
    }
    if not path.exists():
        return defaults
    rows = read_csv(path)
    for row in rows:
        gene = clean(row.get("gene")).upper()
        if not gene:
            continue
        defaults[gene] = {
            "inheritance_model": clean(row.get("inheritance_model")),
            "allowed_causal_variant_count": clean(row.get("allowed_causal_variant_count")) or "1",
            "review_rule": clean(row.get("review_rule")) or "flag_if_more_than_allowed",
            "notes": clean(row.get("notes")),
        }
    return defaults


def load_canonical(canonical_root, auxiliary_root=None):
    rows = []
    summaries = {}
    auxiliary_root = auxiliary_root or default_auxiliary_root(canonical_root)
    for gene in GENES:
        gene_dir = canonical_root / f"{gene}_canonical_patients"
        aux_dir = auxiliary_root / gene
        canonical_csv = find_gene_csv(canonical_root, auxiliary_root, gene, "canonical_patients")
        mentions_csv = find_gene_csv(canonical_root, auxiliary_root, gene, "patient_mentions")
        excluded_csv = find_gene_csv(canonical_root, auxiliary_root, gene, "excluded_mentions")
        links_csv = find_gene_csv(canonical_root, auxiliary_root, gene, "patient_dedup_links")
        summary_json = gene_dir / "dedup_summary.json"
        if not summary_json.exists() and aux_dir != gene_dir:
            summary_json = aux_dir / "dedup_summary.json"
        if not canonical_csv:
            continue
        gene_rows = read_csv(canonical_csv)
        for row in gene_rows:
            row["_gene_dir"] = str(gene_dir)
            rows.append(row)
        summary = {}
        if summary_json.exists():
            summary = json.loads(summary_json.read_text(encoding="utf-8"))
        summary["_mention_rows"] = len(read_csv(mentions_csv)) if mentions_csv else 0
        summary["_excluded_rows"] = len(read_csv(excluded_csv)) if excluded_csv else 0
        summary["_link_rows"] = len(read_csv(links_csv)) if links_csv else 0
        summaries[gene] = summary
    return rows, summaries


def priority_for_group(rows):
    pmids = {clean(r.get("pmid")) for r in rows if clean(r.get("pmid"))}
    flags = []
    text = " | ".join(
        " ".join(
            clean(r.get(col))
            for col in [
                "title",
                "patient_id",
                "originality_evidence",
                "cited_source",
                "location_in_paper",
                "other_phenotypes",
            ]
        ).lower()
        for r in rows
    )
    if len(rows) >= 3:
        flags.append(">=3_rows_same_key")
    if len(pmids) >= 3:
        flags.append(">=3_pmids_same_key")
    re_report_terms = [
        term
        for term in ["previously reported", "previously published", "follow-up", "same patient", "revisited"]
        if term in text
    ]
    if re_report_terms:
        flags.append("explicit_re_report_terms:" + "|".join(re_report_terms))
    context_terms = [
        term
        for term in ["family", "father", "mother", "proband", "index case", "pedigree", "case report", "review"]
        if term in text
    ]
    if context_terms:
        flags.append("context_terms:" + "|".join(context_terms))
    if len(rows) >= 3 or re_report_terms:
        return "high", "; ".join(flags)
    return "medium", "; ".join(flags)


def build_split_audit(canonical_rows, out_dir, output_timestamp=None):
    by_key = defaultdict(list)
    for row in canonical_rows:
        key = clean(row.get("dedup_key"))
        if not key:
            continue
        if clean(row.get("match_reason")) == "manual_merge_rule":
            continue
        by_key[(clean(row.get("gene")), key)].append(row)

    groups = []
    row_rows = []
    pair_rows = []
    index = 1
    for (gene, dedup_key), rows in sorted(by_key.items()):
        if len(rows) < 2:
            continue
        if len({clean(r.get("canonical_patient_id")) for r in rows}) < 2:
            continue
        priority, flags = priority_for_group(rows)
        group_id = f"{gene}_G{index:04d}"
        index += 1
        groups.append({
            "split_group_id": group_id,
            "gene": gene,
            "dedup_key": dedup_key,
            "review_priority": priority,
            "review_flags": flags,
            "row_count": str(len(rows)),
            "unique_pmid_count": str(len({clean(r.get("pmid")) for r in rows if clean(r.get("pmid"))})),
            "pmids": compact(r.get("pmid") for r in rows),
            "pmcids": compact(r.get("pmcid") for r in rows),
            "titles": compact(r.get("title") for r in rows),
            "canonical_patient_ids": compact(r.get("canonical_patient_id") for r in rows),
            "mention_ids": compact(r.get("mention_id") for r in rows),
            "patient_ids": compact(r.get("patient_id") for r in rows),
            "variant_c": compact(r.get("variant_c") for r in rows),
            "variant_p": compact(r.get("variant_p") for r in rows),
            "sex": compact(r.get("sex") for r in rows),
            "age_of_onset": compact(r.get("age_of_onset") for r in rows),
            "phenotype_categories": compact(r.get("phenotype_review_category") for r in rows),
            "seizure_types": compact((r.get("seizure_types") for r in rows), limit=4),
            "response_to_asm": compact(r.get("response_to_asm") for r in rows),
            "source_locations": compact((r.get("location_in_paper") for r in rows), limit=4),
            "recommended_action": "manual_review_before_any_merge",
        })
        for row in rows:
            item = {
                "split_group_id": group_id,
                "split_group_priority": priority,
                "split_group_flags": flags,
            }
            item.update(row)
            row_rows.append(item)
        for left, right in combinations(rows, 2):
            pair_rows.append({
                "split_group_id": group_id,
                "gene": gene,
                "dedup_key": dedup_key,
                "pair_flags": "different_pmids" if clean(left.get("pmid")) != clean(right.get("pmid")) else "same_pmid",
                "left_canonical_patient_id": clean(left.get("canonical_patient_id")),
                "left_mention_id": clean(left.get("mention_id")),
                "left_pmid": clean(left.get("pmid")),
                "left_title": clean(left.get("title")),
                "left_patient_id": clean(left.get("patient_id")),
                "left_location_in_paper": clean(left.get("location_in_paper")),
                "left_phenotype_category": clean(left.get("phenotype_review_category")),
                "left_seizure_types": clean(left.get("seizure_types")),
                "left_other_phenotypes": clean(left.get("other_phenotypes")),
                "right_canonical_patient_id": clean(right.get("canonical_patient_id")),
                "right_mention_id": clean(right.get("mention_id")),
                "right_pmid": clean(right.get("pmid")),
                "right_title": clean(right.get("title")),
                "right_patient_id": clean(right.get("patient_id")),
                "right_location_in_paper": clean(right.get("location_in_paper")),
                "right_phenotype_category": clean(right.get("phenotype_review_category")),
                "right_seizure_types": clean(right.get("seizure_types")),
                "right_other_phenotypes": clean(right.get("other_phenotypes")),
            })

    group_fields = [
        "split_group_id", "gene", "dedup_key", "review_priority", "review_flags",
        "row_count", "unique_pmid_count", "pmids", "pmcids", "titles",
        "canonical_patient_ids", "mention_ids", "patient_ids", "variant_c",
        "variant_p", "sex", "age_of_onset", "phenotype_categories",
        "seizure_types", "response_to_asm", "source_locations", "recommended_action",
    ]
    row_fields = sorted({key for row in row_rows for key in row})
    pair_fields = [
        "split_group_id", "gene", "dedup_key", "pair_flags",
        "left_canonical_patient_id", "left_mention_id", "left_pmid", "left_title",
        "left_patient_id", "left_location_in_paper", "left_phenotype_category",
        "left_seizure_types", "left_other_phenotypes", "right_canonical_patient_id",
        "right_mention_id", "right_pmid", "right_title", "right_patient_id",
        "right_location_in_paper", "right_phenotype_category", "right_seizure_types",
        "right_other_phenotypes",
    ]
    written_files = [
        write_csv(
            out_dir / "current_same_gene_variant_sex_onset_groups.csv",
            groups,
            group_fields,
            output_timestamp,
        ).name,
        write_csv(
            out_dir / "current_same_gene_variant_sex_onset_rows.csv",
            row_rows,
            row_fields,
            output_timestamp,
        ).name,
        write_csv(
            out_dir / "current_same_gene_variant_sex_onset_pairs.csv",
            pair_rows,
            pair_fields,
            output_timestamp,
        ).name,
    ]
    return groups, row_rows, pair_rows, written_files


def build_multi_variant_review(canonical_rows, inheritance_models, out_dir, output_timestamp=None):
    by_patient = defaultdict(list)
    for row in canonical_rows:
        gene = clean(row.get("gene"))
        canonical_patient_id = clean(row.get("canonical_patient_id"))
        if not gene or not canonical_patient_id:
            continue
        by_patient[(gene, canonical_patient_id)].append(row)

    review_rows = []
    for (gene, canonical_patient_id), rows in sorted(by_patient.items()):
        model = inheritance_models.get(gene, {})
        try:
            allowed_count = int(model.get("allowed_causal_variant_count") or "1")
        except ValueError:
            allowed_count = 1
        variant_keys = {causal_variant_key(row) for row in rows if causal_variant_key(row)}
        observed_count = len(variant_keys)
        if observed_count <= allowed_count:
            continue
        review_rows.append({
            "gene": gene,
            "canonical_patient_id": canonical_patient_id,
            "inheritance_model": model.get("inheritance_model", ""),
            "allowed_causal_variant_count": str(allowed_count),
            "observed_variant_count": str(observed_count),
            "review_reason": "observed_variant_count_exceeds_inheritance_model",
            "recommended_action": "manual_review_required_before_canonical_causal_output",
            "pmids": compact(row.get("pmid") for row in rows),
            "patient_ids": compact(row.get("patient_id") for row in rows),
            "mention_ids": compact(row.get("mention_id") for row in rows),
            "variants": compact(
                f"{clean(row.get('variant_c')) or 'Not specified'} / {clean(row.get('variant_p')) or 'Not specified'}"
                for row in rows
            ),
            "clinvar_pathogenicity": compact(row.get("clinvar_pathogenicity") for row in rows),
            "clinvar_ids": compact(row.get("ClinVar ID") or row.get("clinvar_id") for row in rows),
            "phenotype_categories": compact(row.get("phenotype_review_category") for row in rows),
            "source_locations": compact((row.get("location_in_paper") for row in rows), limit=4),
            "notes": model.get("notes", ""),
        })

    fields = [
        "gene",
        "canonical_patient_id",
        "inheritance_model",
        "allowed_causal_variant_count",
        "observed_variant_count",
        "review_reason",
        "recommended_action",
        "pmids",
        "patient_ids",
        "mention_ids",
        "variants",
        "clinvar_pathogenicity",
        "clinvar_ids",
        "phenotype_categories",
        "source_locations",
        "notes",
    ]
    written_files = [
        write_csv(
            out_dir / "multi_variant_patient_review.csv",
            review_rows,
            fields,
            output_timestamp,
        ).name
    ]
    return review_rows, written_files


def qc_canonical(canonical_rows, summaries, out_dir, output_timestamp=None, multi_variant_counts=None):
    summary_rows = []
    issue_rows = []
    per_gene = defaultdict(list)
    for row in canonical_rows:
        per_gene[clean(row.get("gene"))].append(row)

    required = [
        "nm_number",
        "variant_c",
        "variant_p",
        "sex",
        "age_of_onset",
        "phenotype_review_category",
        "seizure_types",
    ]
    for gene in GENES:
        rows = per_gene.get(gene, [])
        summary = summaries.get(gene, {})
        unique_canonical_patients = len({clean(row.get("canonical_patient_id")) for row in rows if clean(row.get("canonical_patient_id"))})
        if not unique_canonical_patients and summary.get("canonical_patients") not in {"", None}:
            unique_canonical_patients = int(summary.get("canonical_patients", 0) or 0)
        counts = Counter()
        for row in rows:
            for col in required:
                if not clean(row.get(col)):
                    counts[f"missing_{col}"] += 1
            if clean(row.get("variant_c")) and clean(row.get("variant_p")):
                if clean(row.get("variant_c")).lower() == "not specified" and clean(row.get("variant_p")).lower() == "not specified":
                    counts["both_variant_c_p_not_specified"] += 1
            clinvar_pathogenicity = clean(row.get("clinvar_pathogenicity"))
            clinvar_id = clean(row.get("ClinVar ID") or row.get("clinvar_id"))
            if clinvar_pathogenicity and not clinvar_id:
                counts["clinvar_pathogenicity_without_id"] += 1
        summary_rows.append({
            "gene": gene,
            "patient_mentions": str(summary.get("_mention_rows", "")),
            "canonical_patients": str(unique_canonical_patients),
            "canonical_patient_variant_rows": str(len(rows)),
            "excluded_mentions": str(summary.get("_excluded_rows", "")),
            "dedup_links": str(summary.get("_link_rows", "")),
            "merged_clusters": str(summary.get("merged_clusters", "")),
            "manual_merge_rules_loaded": str(summary.get("manual_merge_rules_loaded", "")),
            "manual_merge_rules_applied": str(summary.get("manual_merge_rules_applied", "")),
            "multi_variant_review_groups": str((multi_variant_counts or {}).get(gene, 0)),
            **{key: str(counts.get(key, 0)) for key in [
                "missing_nm_number",
                "missing_variant_c",
                "missing_variant_p",
                "missing_sex",
                "missing_age_of_onset",
                "missing_phenotype_review_category",
                "missing_seizure_types",
                "both_variant_c_p_not_specified",
                "clinvar_pathogenicity_without_id",
            ]},
        })
        for row in rows:
            issues = []
            for col in required:
                if not clean(row.get(col)):
                    issues.append(f"missing_{col}")
            clinvar_pathogenicity = clean(row.get("clinvar_pathogenicity"))
            clinvar_id = clean(row.get("ClinVar ID") or row.get("clinvar_id"))
            if clinvar_pathogenicity and not clinvar_id:
                issues.append("clinvar_pathogenicity_without_id")
            if issues:
                issue_rows.append({
                    "gene": gene,
                    "canonical_patient_id": clean(row.get("canonical_patient_id")),
                    "mention_id": clean(row.get("mention_id")),
                    "pmid": clean(row.get("pmid")),
                    "patient_id": clean(row.get("patient_id")),
                    "variant_c": clean(row.get("variant_c")),
                    "variant_p": clean(row.get("variant_p")),
                    "issues": "; ".join(issues),
                    "location_in_paper": clean(row.get("location_in_paper")),
                })
    written_files = [
        write_csv(
            out_dir / "canonical_qc_summary.csv",
            summary_rows,
            list(summary_rows[0]),
            output_timestamp,
        ).name
    ]
    issue_fields = [
        "gene", "canonical_patient_id", "mention_id", "pmid", "patient_id",
        "variant_c", "variant_p", "issues", "location_in_paper",
    ]
    written_files.append(
        write_csv(
            out_dir / "canonical_qc_issues.csv",
            issue_rows,
            issue_fields,
            output_timestamp,
        ).name
    )
    return summary_rows, issue_rows, written_files


def write_markdown(out_dir, groups, summary_rows, issue_rows, csv_files, multi_variant_rows=None):
    priority_counts = Counter(row["review_priority"] for row in groups)
    gene_counts = Counter(row["gene"] for row in groups)
    multi_variant_rows = multi_variant_rows or []
    multi_variant_counts = Counter(row["gene"] for row in multi_variant_rows)
    lines = [
        "# Current Canonical QC",
        "",
        "Generated from current `outputs/current/canonical_patients/{GENE}_canonical_patients/*_canonical_patients.csv` outputs.",
        "",
        "## Split Candidates",
        "",
        f"- Candidate groups: `{len(groups)}`",
        f"- high: `{priority_counts.get('high', 0)}`",
        f"- medium: `{priority_counts.get('medium', 0)}`",
        "",
        "| Gene | groups |",
        "|---|---:|",
    ]
    for gene in GENES:
        lines.append(f"| {gene} | {gene_counts.get(gene, 0)} |")
    lines.extend([
        "",
        "## Multi-Variant Patients Requiring Manual Review",
        "",
        f"- Candidate patients: `{len(multi_variant_rows)}`",
        "",
        "| Gene | patients |",
        "|---|---:|",
    ])
    for gene in GENES:
        lines.append(f"| {gene} | {multi_variant_counts.get(gene, 0)} |")
    lines.extend([
        "",
        "## Canonical Counts",
        "",
        "| Gene | mentions | canonical patients | patient-variant rows | excluded | links | manual rules | multi-variant review | QC issue rows |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|",
    ])
    issue_counts = Counter(row["gene"] for row in issue_rows)
    for row in summary_rows:
        lines.append(
            f"| {row['gene']} | {row['patient_mentions']} | {row['canonical_patients']} | "
            f"{row['canonical_patient_variant_rows']} | {row['excluded_mentions']} | {row['dedup_links']} | "
            f"{row['manual_merge_rules_applied']} | {row.get('multi_variant_review_groups', 0)} | "
            f"{issue_counts.get(row['gene'], 0)} |"
        )
    lines.extend([
        "",
        "## Output Files",
        "",
    ])
    for name in csv_files:
        lines.append(f"- `{name}`")
    (out_dir / "current_canonical_qc.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main():
    parser = argparse.ArgumentParser(description="Audit current canonical patient outputs")
    parser.add_argument(
        "--canonical-root",
        type=Path,
        default=Path("outputs/current/canonical_patients"),
        help="Root containing {GENE}_canonical_patients directories",
    )
    parser.add_argument(
        "--auxiliary-root",
        type=Path,
        default=None,
        help=(
            "Optional root containing auxiliary per-gene CSVs. Defaults to "
            "outputs/auxiliary/ when present."
        ),
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("runs/qc"),
        help="Report output directory",
    )
    parser.add_argument(
        "--inheritance-config",
        type=Path,
        default=Path("config/gene_inheritance_models.csv"),
        help="CSV containing gene inheritance model and allowed causal variant count.",
    )
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
    args = parser.parse_args()
    output_timestamp = None if args.no_timestamp else make_timestamp(args.timestamp)
    auxiliary_root = args.auxiliary_root or default_auxiliary_root(args.canonical_root)
    inheritance_models = load_inheritance_models(args.inheritance_config)

    canonical_rows, summaries = load_canonical(args.canonical_root, auxiliary_root)
    groups, _, _, split_files = build_split_audit(canonical_rows, args.output_dir, output_timestamp)
    multi_variant_rows, multi_variant_files = build_multi_variant_review(
        canonical_rows,
        inheritance_models,
        args.output_dir,
        output_timestamp,
    )
    multi_variant_counts = Counter(row["gene"] for row in multi_variant_rows)
    summary_rows, issue_rows, qc_files = qc_canonical(
        canonical_rows,
        summaries,
        args.output_dir,
        output_timestamp,
        multi_variant_counts,
    )
    write_markdown(
        args.output_dir,
        groups,
        summary_rows,
        issue_rows,
        split_files + multi_variant_files + qc_files,
        multi_variant_rows,
    )

    print("Canonical QC complete")
    print(f"  canonical_rows: {len(canonical_rows)}")
    print(f"  split_candidate_groups: {len(groups)}")
    print(f"  multi_variant_review_groups: {len(multi_variant_rows)}")
    print(f"  qc_issue_rows: {len(issue_rows)}")
    print(f"  output_dir: {args.output_dir}")


if __name__ == "__main__":
    main()
