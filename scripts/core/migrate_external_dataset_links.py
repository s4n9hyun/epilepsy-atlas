#!/usr/bin/env python3
"""Create external_dataset_links tables in gene SQLite databases."""

import argparse
import sqlite3
from datetime import datetime
from pathlib import Path


SCHEMA_SQL = """
CREATE TABLE IF NOT EXISTS external_dataset_links (
    link_id INTEGER PRIMARY KEY AUTOINCREMENT,
    pmid TEXT NOT NULL,
    pmcid TEXT,
    gene TEXT,
    source_type TEXT NOT NULL DEFAULT 'external_dataset',
    repository TEXT,
    url TEXT NOT NULL,
    doi TEXT,
    link_text TEXT,
    evidence_location TEXT,
    expected_content TEXT,
    access_status TEXT NOT NULL DEFAULT 'unknown',
    manual_download_required INTEGER NOT NULL DEFAULT 0,
    manual_download_reason TEXT,
    extraction_status TEXT NOT NULL DEFAULT 'not_processed',
    local_file_path TEXT,
    notes TEXT,
    discovered_at TEXT NOT NULL,
    updated_at TEXT NOT NULL,
    UNIQUE(pmid, url)
);

CREATE INDEX IF NOT EXISTS idx_external_dataset_links_pmid
    ON external_dataset_links (pmid);
CREATE INDEX IF NOT EXISTS idx_external_dataset_links_pmcid
    ON external_dataset_links (pmcid);
CREATE INDEX IF NOT EXISTS idx_external_dataset_links_manual
    ON external_dataset_links (manual_download_required);
CREATE INDEX IF NOT EXISTS idx_external_dataset_links_extraction
    ON external_dataset_links (extraction_status);
"""


SEED_SYNGAP1_30541864 = {
    "pmid": "30541864",
    "pmcid": "PMC6340340",
    "gene": "SYNGAP1",
    "source_type": "external_dataset",
    "repository": "Dryad",
    "url": "https://datadryad.org/dataset/doi:10.5061/dryad.ck70sj0",
    "doi": "10.5061/dryad.ck70sj0",
    "link_text": "Data: SYNGAP1 encephalopathy: a distinctive generalized developmental and epileptic encephalopathy",
    "evidence_location": "Article data availability / external Dryad dataset link",
    "expected_content": "Patient-level genotype/phenotype table; reported to include more than 60 individuals",
    "access_status": "manual_download_needed",
    "manual_download_required": 1,
    "manual_download_reason": "Dataset is outside the PMC supplementary package and must be downloaded from Dryad separately.",
    "extraction_status": "not_processed",
    "local_file_path": "",
    "notes": "Created after identifying external patient/variant data source not captured by supplementary-file download.",
}


def create_table(db_path: Path) -> None:
    with sqlite3.connect(db_path) as conn:
        conn.executescript(SCHEMA_SQL)


def seed_known_link(db_path: Path) -> None:
    now = datetime.now().isoformat(timespec="seconds")
    row = dict(SEED_SYNGAP1_30541864)
    row["discovered_at"] = now
    row["updated_at"] = now

    with sqlite3.connect(db_path) as conn:
        conn.execute(
            """
            INSERT INTO external_dataset_links (
                pmid, pmcid, gene, source_type, repository, url, doi,
                link_text, evidence_location, expected_content, access_status,
                manual_download_required, manual_download_reason, extraction_status,
                local_file_path, notes, discovered_at, updated_at
            )
            VALUES (
                :pmid, :pmcid, :gene, :source_type, :repository, :url, :doi,
                :link_text, :evidence_location, :expected_content, :access_status,
                :manual_download_required, :manual_download_reason, :extraction_status,
                :local_file_path, :notes, :discovered_at, :updated_at
            )
            ON CONFLICT(pmid, url) DO UPDATE SET
                pmcid = excluded.pmcid,
                gene = excluded.gene,
                source_type = excluded.source_type,
                repository = excluded.repository,
                doi = excluded.doi,
                link_text = excluded.link_text,
                evidence_location = excluded.evidence_location,
                expected_content = excluded.expected_content,
                access_status = excluded.access_status,
                manual_download_required = excluded.manual_download_required,
                manual_download_reason = excluded.manual_download_reason,
                extraction_status = excluded.extraction_status,
                notes = excluded.notes,
                updated_at = excluded.updated_at
            """,
            row,
        )


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Create external_dataset_links tables in gene SQLite databases."
    )
    parser.add_argument(
        "--db-dir",
        default="databases",
        help="Directory containing gene SQLite databases.",
    )
    parser.add_argument(
        "--seed-known",
        action="store_true",
        help="Seed currently known external dataset links.",
    )
    args = parser.parse_args()

    db_dir = Path(args.db_dir)
    db_paths = sorted(db_dir.glob("*.sqlite"))
    if not db_paths:
        raise SystemExit(f"No SQLite databases found in {db_dir}")

    for db_path in db_paths:
        create_table(db_path)
        if args.seed_known and db_path.name == "SYNGAP1.sqlite":
            seed_known_link(db_path)
        print(f"updated {db_path}")


if __name__ == "__main__":
    main()
