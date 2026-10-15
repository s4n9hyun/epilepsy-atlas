#!/usr/bin/env python3
"""
article_files/{PMID} 폴더의 추출 후보 파일을 유전자별 DB에 기록한다.

DB에는 PDF와 supplementary 문서/표 파일의 상대 경로만 저장하고,
이미지/XML처럼 현재 추출에 쓰지 않는 파일은 제외한다.
"""

import argparse
import sqlite3
from datetime import datetime
from pathlib import Path

SUPPORTED_EXTENSIONS = {
    '.pdf', '.docx', '.doc', '.pptx', '.ppt',
    '.xlsx', '.xls', '.csv', '.tsv', '.txt', '.md'
}


def is_supplementary_file(filename):
    nl = filename.lower()
    supp_keywords = [
        'supplement', 'supp', 'esm', 'moesm', 'mmc',
        'data_sheet', 'checklist', 'appendix', 'supporting',
        'suppl', 'suplet', 'additional',
    ]
    if any(k in nl for k in supp_keywords):
        return True
    return False


def ensure_schema(conn):
    conn.executescript("""
        CREATE TABLE IF NOT EXISTS article_files (
            file_id INTEGER PRIMARY KEY AUTOINCREMENT,
            article_id TEXT,
            pmid TEXT,
            pmcid TEXT,
            file_name TEXT,
            relative_path TEXT,
            file_type TEXT,
            is_supplementary INTEGER,
            file_size INTEGER,
            indexed_at TEXT
        );

        CREATE INDEX IF NOT EXISTS idx_article_files_article_id
            ON article_files (article_id);
        CREATE INDEX IF NOT EXISTS idx_article_files_pmid
            ON article_files (pmid);
    """)


def index_database(db_path, article_root, workspace):
    conn = sqlite3.connect(db_path)
    conn.row_factory = sqlite3.Row
    ensure_schema(conn)

    rows = conn.execute("""
        SELECT pmid, pmcid
        FROM pubmed_articles
        WHERE downloaded = 'Yes' AND COALESCE(pmid, '') != ''
    """).fetchall()

    indexed_at = datetime.now().isoformat()
    file_count = 0

    with conn:
        conn.execute("DELETE FROM article_files")
        for row in rows:
            pmid = str(row["pmid"])
            pmcid = row["pmcid"] or ""
            article_dir = article_root / pmid
            if not article_dir.exists():
                continue

            for path in sorted(article_dir.rglob("*")):
                suffix = path.suffix.lower()
                if not path.is_file() or suffix not in SUPPORTED_EXTENSIONS:
                    continue

                relative_path = path.relative_to(workspace).as_posix()
                relative_article_path = path.relative_to(article_dir).as_posix()

                conn.execute("""
                    INSERT INTO article_files (
                        article_id, pmid, pmcid, file_name, relative_path,
                        file_type, is_supplementary, file_size, indexed_at
                    )
                    VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
                """, (
                    pmid,
                    pmid,
                    pmcid,
                    path.name,
                    relative_path,
                    suffix.lstrip("."),
                    int(is_supplementary_file(relative_article_path)),
                    path.stat().st_size,
                    indexed_at,
                ))
                file_count += 1

    conn.close()
    return len(rows), file_count


def main():
    parser = argparse.ArgumentParser(description="Index extractable article files into gene SQLite DBs")
    parser.add_argument("--db-dir", default="databases", help="유전자별 SQLite DB 디렉토리")
    parser.add_argument("--article-dir", default="article_files", help="PMID별 논문 파일 디렉토리")
    args = parser.parse_args()

    workspace = Path(__file__).resolve().parents[2]
    db_dir = workspace / args.db_dir
    article_root = workspace / args.article_dir

    for db_path in sorted(db_dir.glob("*.sqlite")):
        article_count, file_count = index_database(db_path, article_root, workspace)
        print(f"{db_path.name}: articles={article_count}, files={file_count}")


if __name__ == "__main__":
    main()
