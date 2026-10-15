#!/usr/bin/env python3
"""patient_mentions 전체를 CSV로 export한다."""

import argparse
import sqlite3
from pathlib import Path

from timestamped_outputs import make_timestamp, timestamp_path


def ordered_patient_columns(db_path):
    conn = sqlite3.connect(db_path)
    all_columns = [row[1] for row in conn.execute("PRAGMA table_info(patient_mentions)")]
    conn.close()

    leading_columns = [
        "nm_number",
        "variant_c_source",
        "variant_c",
        "variant_p",
        "variant_type",
        "clinvar_pathogenicity",
        "ClinVar ID",
        "patient_id",
        "pmid",
        "korean_patient_source",
        "korean_patient_evidence",
    ]
    trailing_columns = [
        "gnomad_conversion_status",
        "nm_transcript",
        "gnomad_transcript",
        "pmcid",
        "mention_id",
    ]
    middle_columns = [
        col for col in all_columns
        if col not in leading_columns
        and col not in trailing_columns
        and col not in {"variant_c_gnomad", "variant_c_original", "clinvar", "clinvar_id"}
    ]
    output_only_columns = {"variant_c", "variant_c_source"}
    return [
        col for col in leading_columns
        if col in output_only_columns
        or col in all_columns
        or (col == "ClinVar ID" and "clinvar_id" in all_columns)
    ] + middle_columns + [
        col for col in trailing_columns
        if col in all_columns
    ]


def export_csv(db_path, output_path, columns):
    conn = sqlite3.connect(db_path)
    nm_number_expr = """
        CASE
            WHEN variant_c_gnomad IS NOT NULL
             AND instr(variant_c_gnomad, ':') > 0
             AND substr(variant_c_gnomad, 1, instr(variant_c_gnomad, ':') - 1) LIKE 'NM_%'
            THEN substr(variant_c_gnomad, 1, instr(variant_c_gnomad, ':') - 1)
            WHEN nm_number IS NOT NULL
             AND trim(nm_number) != ''
             AND lower(trim(nm_number)) NOT IN ('not specified', 'n/a', 'na', 'none', 'unknown')
            THEN nm_number
            WHEN nm_transcript IS NOT NULL
             AND trim(nm_transcript) != ''
             AND lower(trim(nm_transcript)) NOT IN ('not specified', 'n/a', 'na', 'none', 'unknown')
            THEN nm_transcript
            WHEN gnomad_transcript IS NOT NULL
             AND trim(gnomad_transcript) != ''
             AND lower(trim(gnomad_transcript)) NOT IN ('not specified', 'n/a', 'na', 'none', 'unknown')
            THEN gnomad_transcript
            WHEN variant_c_original IS NOT NULL
             AND instr(variant_c_original, ':') > 0
             AND substr(variant_c_original, 1, instr(variant_c_original, ':') - 1) LIKE 'NM_%'
            THEN substr(variant_c_original, 1, instr(variant_c_original, ':') - 1)
            ELSE ''
        END
    """
    variant_c_expr = """
        CASE
            WHEN variant_c_gnomad IS NOT NULL
             AND trim(variant_c_gnomad) != ''
             AND lower(trim(variant_c_gnomad)) NOT IN ('not specified', 'n/a', 'na', 'none', 'unknown')
            THEN
                CASE
                    WHEN instr(variant_c_gnomad, ':') > 0
                    THEN substr(variant_c_gnomad, instr(variant_c_gnomad, ':') + 1)
                    ELSE variant_c_gnomad
                END
            WHEN variant_c_original IS NOT NULL
             AND trim(variant_c_original) != ''
             AND lower(trim(variant_c_original)) NOT IN ('not specified', 'n/a', 'na', 'none', 'unknown')
            THEN
                CASE
                    -- transcript 접두어(NM_/NC_/ENST…:)만 떼고, chrX: 같은 좌표 표기는 그대로 둔다
                    WHEN instr(variant_c_original, ':') > 0
                     AND (upper(substr(variant_c_original, 1, 3)) IN ('NM_', 'NC_')
                          OR upper(substr(variant_c_original, 1, 4)) = 'ENST')
                    THEN substr(variant_c_original, instr(variant_c_original, ':') + 1)
                    ELSE variant_c_original
                END
            ELSE 'Not specified'
        END
    """
    variant_c_source_expr = """
        CASE
            WHEN variant_c_gnomad IS NOT NULL
             AND trim(variant_c_gnomad) != ''
             AND lower(trim(variant_c_gnomad)) NOT IN ('not specified', 'n/a', 'na', 'none', 'unknown')
            THEN 'converted_with_nm_number'
            WHEN variant_c_original IS NOT NULL
             AND trim(variant_c_original) != ''
             AND lower(trim(variant_c_original)) NOT IN ('not specified', 'n/a', 'na', 'none', 'unknown')
            THEN 'raw_text'
            ELSE 'not_specified'
        END
    """
    select_exprs = []
    for col in columns:
        if col == "nm_number":
            select_exprs.append(f"{nm_number_expr} AS nm_number")
        elif col == "variant_c":
            select_exprs.append(f"{variant_c_expr} AS variant_c")
        elif col == "variant_c_source":
            select_exprs.append(f"{variant_c_source_expr} AS variant_c_source")
        elif col == "clinvar_pathogenicity":
            select_exprs.append(
                """
                CASE
                    WHEN clinvar_pathogenicity IN (
                        'Pathogenic',
                        'Likely pathogenic',
                        'Pathogenic/Likely pathogenic'
                    )
                    THEN clinvar_pathogenicity
                    ELSE ''
                END AS clinvar_pathogenicity
                """
            )
        elif col == "ClinVar ID":
            select_exprs.append(
                """
                CASE
                    WHEN clinvar_pathogenicity IN (
                        'Pathogenic',
                        'Likely pathogenic',
                        'Pathogenic/Likely pathogenic'
                    )
                     AND clinvar_id IS NOT NULL
                    THEN clinvar_id
                    ELSE ''
                END AS "ClinVar ID"
                """
            )
        else:
            select_exprs.append(col)
    select_cols = ", ".join(select_exprs)
    query = f"""
        SELECT {select_cols}
        FROM patient_mentions
        ORDER BY variant_c, variant_p, patient_id, pmid
    """

    rows = conn.execute(query)
    with open(output_path, "w", encoding="utf-8", newline="") as f:
        import csv
        writer = csv.writer(f)
        writer.writerow(columns)
        writer.writerows(rows)
    conn.close()


def export_patient_mentions(db_path, output_path=None, output_timestamp=None):
    db_path = Path(db_path)
    gene = db_path.stem
    output = Path(output_path) if output_path else timestamp_path(
        Path(f"{gene}_patient_mentions.csv"),
        output_timestamp,
    )
    columns = ordered_patient_columns(db_path)
    export_csv(db_path, output, columns)
    return output


def main():
    parser = argparse.ArgumentParser(description="Export full patient_mentions CSV")
    parser.add_argument("db_path", help="SQLite DB path")
    parser.add_argument("--output", default=None, help="CSV 출력 경로")
    parser.add_argument(
        "--timestamp",
        default=None,
        help="CSV filename timestamp in YYYYMMDD_HHMMSS format. Defaults to current time.",
    )
    parser.add_argument(
        "--no-timestamp",
        action="store_true",
        help="Write legacy CSV filename without timestamp when --output is not set.",
    )
    args = parser.parse_args()

    output = export_patient_mentions(
        args.db_path,
        args.output,
        output_timestamp=None if args.no_timestamp or args.output else make_timestamp(args.timestamp),
    )
    print(output)


if __name__ == "__main__":
    main()
