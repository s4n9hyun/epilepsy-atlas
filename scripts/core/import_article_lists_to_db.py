#!/usr/bin/env python3
"""
Import gene-specific PubMed article list CSVs into gene SQLite databases.

Usage:
    python scripts/core/import_article_lists_to_db.py data/pubmed_article_lists_epilepsy --db-dir databases
"""

import argparse
import csv
import sqlite3
from datetime import datetime
from pathlib import Path


def infer_gene_from_filename(path):
    name = path.name
    if "_pubmed_" in name:
        return name.split("_pubmed_", 1)[0].upper()
    return path.stem.split("_", 1)[0].upper()


def clean(value):
    if value is None:
        return ""
    return str(value).strip()


def init_schema(conn):
    conn.executescript("""
        CREATE TABLE IF NOT EXISTS pubmed_articles (
            pmid TEXT PRIMARY KEY,
            pmcid TEXT,
            gene TEXT NOT NULL,
            title TEXT,
            full_package_available TEXT,
            downloaded TEXT,
            extracted TEXT NOT NULL DEFAULT 'No',
            extraction_status TEXT,
            extracted_at TEXT,
            update_date TEXT,
            source_csv TEXT,
            imported_at TEXT
        );
    """)

    existing_cols = {
        row[1] for row in conn.execute("PRAGMA table_info(pubmed_articles)").fetchall()
    }
    migrations = {
        "extracted": "ALTER TABLE pubmed_articles ADD COLUMN extracted TEXT NOT NULL DEFAULT 'No'",
        "extraction_status": "ALTER TABLE pubmed_articles ADD COLUMN extraction_status TEXT",
        "extracted_at": "ALTER TABLE pubmed_articles ADD COLUMN extracted_at TEXT",
    }
    for col, sql in migrations.items():
        if col not in existing_cols:
            conn.execute(sql)

    conn.executescript("""
        CREATE INDEX IF NOT EXISTS idx_pubmed_articles_pmcid
            ON pubmed_articles (pmcid);
        CREATE INDEX IF NOT EXISTS idx_pubmed_articles_downloaded
            ON pubmed_articles (downloaded);
        CREATE INDEX IF NOT EXISTS idx_pubmed_articles_extracted
            ON pubmed_articles (extracted);
        CREATE INDEX IF NOT EXISTS idx_pubmed_articles_full_package
            ON pubmed_articles (full_package_available);
    """)

    conn.commit()


def import_csv(csv_path, db_dir):
    gene = infer_gene_from_filename(csv_path)
    db_path = db_dir / f"{gene}.sqlite"
    db_path.parent.mkdir(parents=True, exist_ok=True)

    imported_at = datetime.now().isoformat()
    inserted_or_updated = 0

    conn = sqlite3.connect(db_path)
    try:
        init_schema(conn)
        with open(csv_path, newline="", encoding="utf-8-sig") as f:
            reader = csv.DictReader(f)
            with conn:
                for row in reader:
                    pmid = clean(row.get("pmid"))
                    if not pmid:
                        continue
                    conn.execute("""
                        INSERT INTO pubmed_articles (
                            pmid, pmcid, gene, title, full_package_available,
                            downloaded, update_date, source_csv, imported_at
                        )
                        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
                        ON CONFLICT(pmid) DO UPDATE SET
                            pmcid=excluded.pmcid,
                            gene=excluded.gene,
                            title=excluded.title,
                            full_package_available=excluded.full_package_available,
                            downloaded=excluded.downloaded,
                            update_date=excluded.update_date,
                            source_csv=excluded.source_csv,
                            imported_at=excluded.imported_at
                    """, (
                        pmid,
                        clean(row.get("pmcid")),
                        gene,
                        clean(row.get("title")),
                        clean(row.get("full_package_available")),
                        clean(row.get("downloaded")),
                        clean(row.get("update_date")),
                        str(csv_path),
                        imported_at,
                    ))
                    inserted_or_updated += 1
    finally:
        conn.close()

    return {
        "gene": gene,
        "db_path": db_path,
        "rows": inserted_or_updated,
    }


def main():
    parser = argparse.ArgumentParser(description="Import gene PubMed article CSVs into SQLite DBs")
    parser.add_argument("csv_dir", type=Path, help="Directory containing *_pubmed_epilepsy_articles.csv")
    parser.add_argument("--db-dir", type=Path, default=Path("databases"), help="SQLite DB output directory")
    args = parser.parse_args()

    csv_files = sorted(args.csv_dir.glob("*.csv"))
    if not csv_files:
        raise SystemExit(f"No CSV files found in {args.csv_dir}")

    total = 0
    for csv_path in csv_files:
        result = import_csv(csv_path, args.db_dir)
        total += result["rows"]
        print(f"{result['gene']}: {result['rows']} rows -> {result['db_path']}")

    print(f"Total imported/updated rows: {total}")


if __name__ == "__main__":
    main()
