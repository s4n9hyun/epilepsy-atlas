#!/usr/bin/env python3
"""유전자별 SQLite DB 상태를 논문 기준으로 요약한다."""

import argparse
import sqlite3
from pathlib import Path


def table_exists(conn, table_name):
    row = conn.execute("""
        SELECT 1
        FROM sqlite_master
        WHERE type = 'table' AND name = ?
    """, (table_name,)).fetchone()
    return row is not None


def count_where(conn, table, where_clause):
    return conn.execute(f"SELECT COUNT(*) FROM {table} WHERE {where_clause}").fetchone()[0]


def summarize_db(db_path, include_files=False):
    conn = sqlite3.connect(db_path)

    total = count_where(conn, "pubmed_articles", "1 = 1")
    downloaded = count_where(conn, "pubmed_articles", "downloaded = 'Yes'")
    extracted = count_where(conn, "pubmed_articles", "extracted = 'Yes'")

    patient_mentions = 0
    processed_articles = 0
    if table_exists(conn, "patient_mentions"):
        patient_mentions = count_where(conn, "patient_mentions", "1 = 1")
    if table_exists(conn, "processed_articles"):
        processed_articles = count_where(conn, "processed_articles", "1 = 1")

    summary = {
        "gene": db_path.stem,
        "total": total,
        "downloaded": downloaded,
        "not_downloaded": total - downloaded,
        "extracted": extracted,
        "processed_articles": processed_articles,
        "patient_mentions": patient_mentions,
    }

    if include_files and table_exists(conn, "article_files"):
        summary["files"] = count_where(conn, "article_files", "1 = 1")

    conn.close()
    return summary


def main():
    parser = argparse.ArgumentParser(description="Show gene DB status by article counts")
    parser.add_argument("--db-dir", default="databases", help="유전자별 SQLite DB 디렉토리")
    parser.add_argument("--include-files", action="store_true", help="article_files 파일 row 수도 함께 표시")
    args = parser.parse_args()

    db_paths = sorted(Path(args.db_dir).glob("*.sqlite"))
    if not db_paths:
        print(f"No SQLite DBs found in {args.db_dir}")
        return

    headers = [
        "gene", "total", "downloaded", "not_downloaded",
        "extracted", "processed_articles", "patient_mentions",
    ]
    if args.include_files:
        headers.append("files")

    print("\t".join(headers))
    for db_path in db_paths:
        summary = summarize_db(db_path, include_files=args.include_files)
        print("\t".join(str(summary.get(header, 0)) for header in headers))


if __name__ == "__main__":
    main()
