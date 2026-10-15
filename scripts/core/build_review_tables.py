#!/usr/bin/env python3
"""
Export patient_mentions from a gene SQLite DB and build manual-review CSVs.

One command performs:
  1. DB -> outputs/working_canonical/{GENE}_canonical_patients/{GENE}_patient_mentions.csv
  2. patient_mentions CSV -> canonical patient candidates
  3. manual-review sorting by variant_c_gnomad / variant_p / source type / PMID
"""

import argparse
import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent

if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from deduplicate_patients import build_outputs
from export_patient_mentions import export_patient_mentions
from timestamped_outputs import make_timestamp, timestamp_path


def build_review_tables(
    db_path,
    output_dir=None,
    patient_csv=None,
    gene=None,
    experimental_strict_source_filter=False,
    output_timestamp=None,
):
    db_path = Path(db_path)
    gene_name = (gene or db_path.stem).upper()
    output_root = Path(output_dir or "outputs/working_canonical")
    output_root.mkdir(parents=True, exist_ok=True)

    canonical_dir = output_root / f"{gene_name}_canonical_patients"
    canonical_dir.mkdir(parents=True, exist_ok=True)
    patient_csv_path = (
        Path(patient_csv)
        if patient_csv
        else timestamp_path(canonical_dir / f"{gene_name}_patient_mentions.csv", output_timestamp)
    )

    exported_csv = export_patient_mentions(db_path, patient_csv_path)
    summary = build_outputs(
        exported_csv,
        canonical_dir,
        gene=gene_name,
        experimental_strict_source_filter=experimental_strict_source_filter,
        output_timestamp=output_timestamp,
    )

    return {
        "gene": gene_name,
        "patient_mentions_csv": str(exported_csv),
        "canonical_dir": str(canonical_dir),
        **summary,
    }


def main():
    parser = argparse.ArgumentParser(
        description="Export DB patient_mentions and build canonical manual-review CSVs"
    )
    parser.add_argument("db_path", type=Path, help="SQLite DB path, e.g. databases/SYNGAP1.sqlite")
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("outputs/working_canonical"),
        help="Directory for exported CSV and canonical review folder",
    )
    parser.add_argument(
        "--patient-csv",
        type=Path,
        default=None,
        help="Optional explicit path for exported patient_mentions CSV",
    )
    parser.add_argument("--gene", default=None, help="Gene name if DB filename is not the gene")
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

    summary = build_review_tables(
        args.db_path,
        output_dir=args.output_dir,
        patient_csv=args.patient_csv,
        gene=args.gene,
        experimental_strict_source_filter=args.experimental_strict_source_filter,
        output_timestamp=None if args.no_timestamp else make_timestamp(args.timestamp),
    )

    print("Review table build complete")
    print(f"  gene: {summary['gene']}")
    print(f"  patient_mentions_csv: {summary['patient_mentions_csv']}")
    print(f"  canonical_dir: {summary['canonical_dir']}")
    print(f"  total_mentions: {summary['total_mentions']}")
    print(f"  canonical_patients: {summary['canonical_patients']}")
    print(f"  canonical_patient_variant_rows: {summary.get('canonical_patient_variant_rows', '')}")
    print(f"  excluded_mentions: {summary['excluded_mentions']}")
    print(f"  linked_cited_mentions: {summary['linked_cited_mentions']}")
    print(f"  merged_clusters: {summary['merged_clusters']}")


if __name__ == "__main__":
    main()
