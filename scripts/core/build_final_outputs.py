#!/usr/bin/env python3
"""Build a consolidated final output package from current project outputs."""

import argparse
import csv
import shutil
import sqlite3
from collections import Counter, defaultdict
from datetime import date
from pathlib import Path

from timestamped_outputs import make_timestamp, resolve_csv_path, resolve_latest_path, timestamp_path


GENES = [
    "CDKL5", "EEF1A2", "KCNQ2", "SCN1A", "SCN2A", "SYNGAP1",
    "PCDH19", "SCN8A", "GABRA1", "STXBP1", "SMC1A", "FOXG1", "KCNT1",
]


def read_csv(path):
    with path.open(newline="", encoding="utf-8-sig") as f:
        return list(csv.DictReader(f))


def write_csv(path, rows, fieldnames=None, output_timestamp=None):
    path = timestamp_path(path, output_timestamp)
    path.parent.mkdir(parents=True, exist_ok=True)
    if fieldnames is None:
        fieldnames = list(rows[0]) if rows else []
    with path.open("w", newline="", encoding="utf-8-sig") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    return path


def union_fieldnames(rows, preferred=None):
    fields = []
    seen = set()
    for name in preferred or []:
        if name and name not in seen:
            fields.append(name)
            seen.add(name)
    for row in rows:
        for name in row:
            if name not in seen:
                fields.append(name)
                seen.add(name)
    return fields


def add_source_gene(rows, gene):
    out = []
    for row in rows:
        item = dict(row)
        if not item.get("gene"):
            item["gene"] = gene
        item["source_gene"] = gene
        out.append(item)
    return out


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


def combine_gene_csvs(canonical_root, suffix, auxiliary_root=None):
    rows = []
    preferred = []
    auxiliary_root = auxiliary_root or default_auxiliary_root(canonical_root)
    for gene in GENES:
        path = find_gene_csv(canonical_root, auxiliary_root, gene, suffix)
        if not path:
            continue
        gene_rows = read_csv(path)
        if not preferred and gene_rows:
            preferred = ["source_gene"] + list(gene_rows[0])
        rows.extend(add_source_gene(gene_rows, gene))
    return rows, union_fieldnames(rows, preferred)


def fetch_external_dataset_links(db_dir):
    rows = []
    preferred = [
        "source_gene",
        "downloaded",
        "extracted",
        "source_db",
        "link_id",
        "pmid",
        "pmcid",
        "gene",
        "source_type",
        "repository",
        "url",
        "doi",
        "link_text",
        "evidence_location",
        "expected_content",
        "access_status",
        "manual_download_required",
        "manual_download_reason",
        "extraction_status",
        "external_dataset_extracted",
        "local_file_path",
        "notes",
        "discovered_at",
        "updated_at",
    ]
    for gene in GENES:
        db_path = db_dir / f"{gene}.sqlite"
        if not db_path.exists():
            continue
        con = sqlite3.connect(db_path)
        con.row_factory = sqlite3.Row
        try:
            db_rows = [
                dict(row)
                for row in con.execute(
                    """
                    select link_id, pmid, pmcid, gene, source_type, repository, url,
                           doi, link_text, evidence_location, expected_content,
                           access_status, manual_download_required,
                           manual_download_reason, extraction_status, local_file_path,
                           notes, discovered_at, updated_at
                    from external_dataset_links
                    order by pmid, link_id
                    """
                )
            ]
        except sqlite3.OperationalError:
            db_rows = []
        finally:
            con.close()
        for row in db_rows:
            row["source_gene"] = gene
            row["source_db"] = str(db_path)
            if not row.get("gene"):
                row["gene"] = gene
            row["downloaded"] = (
                "Yes"
                if str(row.get("access_status", "")).lower() == "downloaded"
                or bool(row.get("local_file_path"))
                else "No"
            )
            row["extracted"] = (
                "Yes" if str(row.get("extraction_status", "")).lower() == "processed" else "No"
            )
            row["external_dataset_extracted"] = (
                "Yes" if str(row.get("extraction_status", "")).lower() == "processed" else "No"
            )
            rows.append(row)
    return rows, preferred


def copy_if_exists(src, dst, output_timestamp=None):
    src = resolve_latest_path(src)
    if src:
        dst = timestamp_path(dst, output_timestamp) if dst.suffix.lower() == ".csv" else dst
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src, dst)
        return dst.name
    return False


def copy_supporting_reports(out_dir, coverage_dir, qc_dir, canonical_root, output_timestamp=None):
    copied = []
    workspace_root = canonical_root.parent
    manual_merge_rules = workspace_root / "config" / "manual_merge_rules.csv"
    if not manual_merge_rules.exists():
        manual_merge_rules = canonical_root / "manual_merge_rules.csv"
    mapping = [
        (coverage_dir / "coverage_summary.csv", out_dir / "coverage_summary.csv"),
        (coverage_dir / "coverage_gene_rows.csv", out_dir / "coverage_gene_rows.csv"),
        (coverage_dir / "coverage_unique_pmids.csv", out_dir / "coverage_unique_pmids.csv"),
        (
            coverage_dir / "downloaded_not_extracted_gene_rows.csv",
            out_dir / "downloaded_not_extracted_gene_rows.csv",
        ),
        (coverage_dir / "not_downloaded_unique_pmids.csv", out_dir / "not_downloaded_unique_pmids.csv"),
        (coverage_dir / "coverage_qc.md", out_dir / "coverage_qc.md"),
        (qc_dir / "canonical_qc_summary.csv", out_dir / "canonical_qc_summary.csv"),
        (qc_dir / "canonical_qc_issues.csv", out_dir / "canonical_qc_issues.csv"),
        (qc_dir / "canonical_actionable_qc_summary.csv", out_dir / "canonical_actionable_qc_summary.csv"),
        (qc_dir / "canonical_actionable_qc_issues.csv", out_dir / "canonical_actionable_qc_issues.csv"),
        (
            qc_dir / "current_same_gene_variant_sex_onset_groups.csv",
            out_dir / "remaining_split_candidate_groups.csv",
        ),
        (
            qc_dir / "multi_variant_patient_review.csv",
            out_dir / "multi_variant_patient_review.csv",
        ),
        (
            qc_dir / "remaining_split_candidate_decisions_20260608.csv",
            out_dir / "remaining_split_candidate_decisions_20260608.csv",
        ),
        (
            qc_dir / "medium_review_new_manual_rules_applied.csv",
            out_dir / "medium_review_new_manual_rules_applied.csv",
        ),
        (qc_dir / "medium_split_review_20260608.md", out_dir / "medium_split_review_20260608.md"),
        (qc_dir / "current_canonical_qc.md", out_dir / "current_canonical_qc.md"),
        (manual_merge_rules, out_dir / "manual_merge_rules.csv"),
        (canonical_root / "README.md", out_dir / "canonical_patients_README.md"),
    ]
    for src, dst in mapping:
        copied_name = copy_if_exists(src, dst, output_timestamp)
        if copied_name:
            copied.append(copied_name)
    return copied


def count_by_gene(rows):
    counts = Counter(row.get("gene") or row.get("source_gene") or "" for row in rows)
    return {gene: counts.get(gene, 0) for gene in GENES}


def count_unique_canonical_by_gene(rows):
    values = defaultdict(set)
    for row in rows:
        gene = row.get("gene") or row.get("source_gene") or ""
        canonical_id = row.get("canonical_patient_id") or ""
        if gene and canonical_id:
            values[gene].add(canonical_id)
    return {gene: len(values.get(gene, set())) for gene in GENES}


def build_summary_rows(combined):
    rows = []
    patient_counts = count_unique_canonical_by_gene(combined["canonical_patients"])
    variant_row_counts = count_by_gene(combined["canonical_patients"])
    for gene in GENES:
        rows.append(
            {
                "gene": gene,
                "canonical_patients": str(patient_counts[gene]),
                "canonical_patient_variant_rows": str(variant_row_counts[gene]),
                "patient_mentions": str(count_by_gene(combined["patient_mentions"])[gene]),
                "excluded_mentions": str(count_by_gene(combined["excluded_mentions"])[gene]),
                "dedup_links": str(count_by_gene(combined["patient_dedup_links"])[gene]),
                "external_dataset_links": str(count_by_gene(combined["external_dataset_links"])[gene]),
            }
        )
    rows.append(
        {
            "gene": "TOTAL",
            "canonical_patients": str(sum(patient_counts.values())),
            "canonical_patient_variant_rows": str(len(combined["canonical_patients"])),
            "patient_mentions": str(len(combined["patient_mentions"])),
            "excluded_mentions": str(len(combined["excluded_mentions"])),
            "dedup_links": str(len(combined["patient_dedup_links"])),
            "external_dataset_links": str(len(combined["external_dataset_links"])),
        }
    )
    return rows


def read_summary_lookup(path):
    if not path.exists():
        return {}
    return {row["scope"]: row for row in read_csv(path)}


def markdown_table(rows, fields):
    lines = [
        "| " + " | ".join(fields) + " |",
        "| " + " | ".join("---" for _ in fields) + " |",
    ]
    for row in rows:
        lines.append("| " + " | ".join(str(row.get(field, "")) for field in fields) + " |")
    return lines


def write_summary_md(out_dir, summary_rows, coverage_summary, copied_reports, core_files):
    lines = [
        "# Final Consolidated Outputs",
        "",
        f"Generated on {date.today().isoformat()} from the current per-gene SQLite DBs and canonical patient outputs.",
        "",
        "## Core Counts",
        "",
    ]
    lines.extend(
        markdown_table(
            summary_rows,
            [
                "gene",
                "canonical_patients",
                "canonical_patient_variant_rows",
                "patient_mentions",
                "excluded_mentions",
                "dedup_links",
                "external_dataset_links",
            ],
        )
    )
    if coverage_summary:
        lines.extend(["", "## Download And Extraction Coverage", ""])
        coverage_rows = [
            coverage_summary[key]
            for key in ["gene_row", "unique_pmid"]
            if key in coverage_summary
        ]
        lines.extend(
            markdown_table(
                coverage_rows,
                [
                    "scope",
                    "total",
                    "downloaded_and_extracted",
                    "downloaded_not_extracted",
                    "not_downloaded",
                    "with_pmcid",
                    "without_pmcid",
                ],
            )
        )
    lines.extend(
        [
            "",
            "## Included Files",
            "",
        ]
    )
    for name in core_files:
        lines.append(f"- `{name}`")
    if (out_dir / "sample_validation_20260608.md").exists():
        lines.append("- `sample_validation_20260608.md`")
    for name in sorted(copied_reports):
        lines.append(f"- `{name}`")
    lines.extend(
        [
            "",
            "## Notes",
            "",
            "- `patient_mentions` includes all extracted mention rows with canonical status.",
            "- `canonical_patients` is patient-variant grain: the same canonical_patient_id can appear on multiple rows when one patient has multiple variants.",
            "- Unique patient counts are calculated by unique `canonical_patient_id` in `canonical_patients`.",
            "- `excluded_mentions` is the subset excluded from canonical patient counts.",
            "- `dedup_links` maps extracted mention rows to canonical patient IDs.",
            "- Coverage totals include all current DB article rows, including downloaded papers that have not yet been GPT-extracted.",
        ]
    )
    (out_dir / "final_output_summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main():
    parser = argparse.ArgumentParser(description="Build final consolidated output files")
    parser.add_argument("--canonical-root", type=Path, default=Path("outputs/current/canonical_patients"))
    parser.add_argument(
        "--auxiliary-root",
        type=Path,
        default=None,
        help=(
            "Optional root containing auxiliary per-gene CSVs. Defaults to "
            "outputs/auxiliary/ when present."
        ),
    )
    parser.add_argument("--db-dir", type=Path, default=Path("databases"))
    parser.add_argument("--coverage-dir", type=Path, default=Path("coverage_qc_20260608"))
    parser.add_argument(
        "--qc-dir",
        type=Path,
        default=Path("outputs/current/qc"),
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("runs/final_outputs"),
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

    outputs = {
        "canonical_patients": ("canonical_patients", "all_genes_canonical_patients.csv"),
        "patient_mentions": ("patient_mentions", "all_genes_patient_mentions.csv"),
        "excluded_mentions": ("excluded_mentions", "all_genes_excluded_mentions.csv"),
        "patient_dedup_links": ("patient_dedup_links", "all_genes_patient_dedup_links.csv"),
    }

    combined = {}
    core_files = []
    for key, (suffix, filename) in outputs.items():
        rows, fields = combine_gene_csvs(args.canonical_root, suffix, auxiliary_root)
        combined[key] = rows
        core_files.append(write_csv(args.output_dir / filename, rows, fields, output_timestamp).name)

    external_rows, external_fields = fetch_external_dataset_links(args.db_dir)
    combined["external_dataset_links"] = external_rows
    core_files.append(
        write_csv(
            args.output_dir / "all_genes_external_dataset_links.csv",
            external_rows,
            external_fields,
            output_timestamp,
        ).name
    )

    copied = copy_supporting_reports(
        args.output_dir,
        args.coverage_dir,
        args.qc_dir,
        args.canonical_root,
        output_timestamp,
    )
    summary_rows = build_summary_rows(combined)
    core_files.append(write_csv(args.output_dir / "final_output_summary.csv", summary_rows, output_timestamp=output_timestamp).name)
    coverage_summary_path = resolve_latest_path(args.coverage_dir / "coverage_summary.csv")
    coverage_summary = read_summary_lookup(coverage_summary_path) if coverage_summary_path else {}
    write_summary_md(args.output_dir, summary_rows, coverage_summary, copied, core_files)

    print(f"Wrote final outputs to {args.output_dir}")
    for row in summary_rows:
        print(
            f"{row['gene']}: canonical={row['canonical_patients']} "
            f"canonical_variant_rows={row['canonical_patient_variant_rows']} "
            f"mentions={row['patient_mentions']} excluded={row['excluded_mentions']} "
            f"dedup_links={row['dedup_links']} external_links={row['external_dataset_links']}"
        )


if __name__ == "__main__":
    main()
