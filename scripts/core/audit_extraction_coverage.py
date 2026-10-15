#!/usr/bin/env python3
"""Audit download and extraction coverage across per-gene SQLite databases."""

import argparse
import csv
import sqlite3
from collections import Counter, defaultdict
from pathlib import Path


GENES = [
    "CDKL5", "EEF1A2", "KCNQ2", "SCN1A", "SCN2A", "SYNGAP1",
    "PCDH19", "SCN8A", "GABRA1", "STXBP1", "SMC1A", "FOXG1", "KCNT1",
]


def clean(value):
    if value is None:
        return ""
    return str(value).strip()


def yes(value):
    return clean(value).lower() in {"yes", "y", "1", "true"}


def write_csv(path, rows, fieldnames):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8-sig") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def fetch_gene(db_path, gene):
    con = sqlite3.connect(db_path)
    con.row_factory = sqlite3.Row
    articles = [
        dict(row)
        for row in con.execute(
            """
            select pmid, pmcid, gene, title, downloaded, extracted,
                   extraction_status, extracted_at
            from pubmed_articles
            """
        )
    ]
    file_counts = {
        clean(row["pmid"]): dict(row)
        for row in con.execute(
            """
            select pmid,
                   count(*) as file_count,
                   sum(case when is_supplementary=1 then 1 else 0 end) as supplementary_file_count,
                   sum(case when lower(file_name) like '%.pdf' then 1 else 0 end) as pdf_file_count,
                   sum(case when file_name='main.pdf' then 1 else 0 end) as main_pdf_count,
                   sum(case when file_name='main.rendered.pdf' then 1 else 0 end) as rendered_pdf_count
            from article_files
            group by pmid
            """
        )
    }
    # 추출 전 DB(예: 신규 유전자)에는 patient_mentions 테이블이 아직 없다.
    has_mentions = con.execute(
        "select 1 from sqlite_master where type='table' and name='patient_mentions'"
    ).fetchone()
    mention_counts = {} if not has_mentions else {
        clean(row["pmid"]): dict(row)
        for row in con.execute(
            """
            select pmid,
                   count(*) as patient_mention_count,
                   sum(case when mention_source_type='original_patient' then 1 else 0 end) as original_patient_mentions,
                   sum(case when mention_source_type='cited_patient' then 1 else 0 end) as cited_patient_mentions
            from patient_mentions
            group by pmid
            """
        )
    }
    external_counts = {}
    try:
        external_counts = {
            clean(row["pmid"]): dict(row)
            for row in con.execute(
                """
                select pmid, count(*) as external_dataset_link_count
                from external_dataset_links
                group by pmid
                """
            )
        }
    except sqlite3.OperationalError:
        external_counts = {}
    con.close()

    rows = []
    for article in articles:
        pmid = clean(article.get("pmid"))
        files = file_counts.get(pmid, {})
        mentions = mention_counts.get(pmid, {})
        externals = external_counts.get(pmid, {})
        file_count = int(files.get("file_count") or 0)
        mention_count = int(mentions.get("patient_mention_count") or 0)
        downloaded = yes(article.get("downloaded")) or file_count > 0
        extracted = yes(article.get("extracted")) or mention_count > 0
        status = "not_downloaded"
        if downloaded and extracted:
            status = "downloaded_and_extracted"
        elif downloaded and not extracted:
            status = "downloaded_not_extracted"
        elif not downloaded and extracted:
            status = "extracted_without_download_flag"
        rows.append({
            "gene": gene,
            "pmid": pmid,
            "pmcid": clean(article.get("pmcid")),
            "title": clean(article.get("title")),
            "downloaded_flag": clean(article.get("downloaded")),
            "extracted_flag": clean(article.get("extracted")),
            "coverage_status": status,
            "file_count": str(file_count),
            "supplementary_file_count": str(int(files.get("supplementary_file_count") or 0)),
            "pdf_file_count": str(int(files.get("pdf_file_count") or 0)),
            "main_pdf_count": str(int(files.get("main_pdf_count") or 0)),
            "rendered_pdf_count": str(int(files.get("rendered_pdf_count") or 0)),
            "patient_mention_count": str(mention_count),
            "original_patient_mentions": str(int(mentions.get("original_patient_mentions") or 0)),
            "cited_patient_mentions": str(int(mentions.get("cited_patient_mentions") or 0)),
            "external_dataset_link_count": str(int(externals.get("external_dataset_link_count") or 0)),
            "extraction_status": clean(article.get("extraction_status")),
            "extracted_at": clean(article.get("extracted_at")),
        })
    return rows


def unique_pmid_rows(gene_rows):
    grouped = defaultdict(list)
    for row in gene_rows:
        grouped[row["pmid"]].append(row)
    rows = []
    for pmid, items in sorted(grouped.items()):
        genes = sorted({row["gene"] for row in items})
        pmcid = next((row["pmcid"] for row in items if row["pmcid"]), "")
        downloaded = any(row["coverage_status"] in {"downloaded_and_extracted", "downloaded_not_extracted"} for row in items)
        extracted = any(row["coverage_status"] in {"downloaded_and_extracted", "extracted_without_download_flag"} for row in items)
        status = "not_downloaded"
        if downloaded and extracted:
            status = "downloaded_and_extracted"
        elif downloaded and not extracted:
            status = "downloaded_not_extracted"
        elif extracted and not downloaded:
            status = "extracted_without_download_flag"
        rows.append({
            "pmid": pmid,
            "pmcid": pmcid,
            "genes": ";".join(genes),
            "gene_row_count": str(len(items)),
            "coverage_status": status,
            "any_pmcid": "Yes" if pmcid else "No",
            "total_files": str(sum(int(row["file_count"]) for row in items)),
            "total_supplementary_files": str(sum(int(row["supplementary_file_count"]) for row in items)),
            "total_patient_mentions": str(sum(int(row["patient_mention_count"]) for row in items)),
            "total_external_dataset_links": str(sum(int(row["external_dataset_link_count"]) for row in items)),
            "titles": " | ".join(dict.fromkeys(row["title"] for row in items if row["title"])),
        })
    return rows


def build_summary(gene_rows, unique_rows):
    rows = []
    for scope, items in [("gene_row", gene_rows), ("unique_pmid", unique_rows)]:
        total = len(items)
        counts = Counter(row["coverage_status"] for row in items)
        pmcid_count = sum(1 for row in items if row.get("pmcid") or row.get("any_pmcid") == "Yes")
        rows.append({
            "scope": scope,
            "total": str(total),
            "downloaded_and_extracted": str(counts.get("downloaded_and_extracted", 0)),
            "downloaded_not_extracted": str(counts.get("downloaded_not_extracted", 0)),
            "not_downloaded": str(counts.get("not_downloaded", 0)),
            "extracted_without_download_flag": str(counts.get("extracted_without_download_flag", 0)),
            "with_pmcid": str(pmcid_count),
            "without_pmcid": str(total - pmcid_count),
        })
    for gene in GENES:
        items = [row for row in gene_rows if row["gene"] == gene]
        counts = Counter(row["coverage_status"] for row in items)
        rows.append({
            "scope": f"gene:{gene}",
            "total": str(len(items)),
            "downloaded_and_extracted": str(counts.get("downloaded_and_extracted", 0)),
            "downloaded_not_extracted": str(counts.get("downloaded_not_extracted", 0)),
            "not_downloaded": str(counts.get("not_downloaded", 0)),
            "extracted_without_download_flag": str(counts.get("extracted_without_download_flag", 0)),
            "with_pmcid": str(sum(1 for row in items if row["pmcid"])),
            "without_pmcid": str(sum(1 for row in items if not row["pmcid"])),
        })
    return rows


def write_markdown(out_dir, summary_rows):
    lines = [
        "# Extraction Coverage QC",
        "",
        "Coverage is computed from current per-gene SQLite DBs using article flags, indexed article files, and patient mention rows.",
        "",
        "| Scope | total | downloaded+extracted | downloaded not extracted | not downloaded | PMCID | no PMCID |",
        "|---|---:|---:|---:|---:|---:|---:|",
    ]
    for row in summary_rows:
        lines.append(
            f"| {row['scope']} | {row['total']} | {row['downloaded_and_extracted']} | "
            f"{row['downloaded_not_extracted']} | {row['not_downloaded']} | "
            f"{row['with_pmcid']} | {row['without_pmcid']} |"
        )
    lines.extend([
        "",
        "## Files",
        "",
        "- `coverage_summary.csv`",
        "- `coverage_gene_rows.csv`",
        "- `coverage_unique_pmids.csv`",
        "- `downloaded_not_extracted_gene_rows.csv`",
        "- `not_downloaded_unique_pmids.csv`",
    ])
    (out_dir / "coverage_qc.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main():
    parser = argparse.ArgumentParser(description="Audit article download/extraction coverage")
    parser.add_argument("--db-dir", type=Path, default=Path("databases"))
    parser.add_argument("--output-dir", type=Path, default=Path("coverage_qc"))
    args = parser.parse_args()

    gene_rows = []
    for gene in GENES:
        db_path = args.db_dir / f"{gene}.sqlite"
        if db_path.exists():
            gene_rows.extend(fetch_gene(db_path, gene))
    unique_rows = unique_pmid_rows(gene_rows)
    summary_rows = build_summary(gene_rows, unique_rows)

    fields_gene = [
        "gene", "pmid", "pmcid", "title", "downloaded_flag", "extracted_flag",
        "coverage_status", "file_count", "supplementary_file_count", "pdf_file_count",
        "main_pdf_count", "rendered_pdf_count", "patient_mention_count",
        "original_patient_mentions", "cited_patient_mentions",
        "external_dataset_link_count", "extraction_status", "extracted_at",
    ]
    fields_unique = [
        "pmid", "pmcid", "genes", "gene_row_count", "coverage_status", "any_pmcid",
        "total_files", "total_supplementary_files", "total_patient_mentions",
        "total_external_dataset_links", "titles",
    ]
    fields_summary = [
        "scope", "total", "downloaded_and_extracted", "downloaded_not_extracted",
        "not_downloaded", "extracted_without_download_flag", "with_pmcid", "without_pmcid",
    ]
    write_csv(args.output_dir / "coverage_gene_rows.csv", gene_rows, fields_gene)
    write_csv(args.output_dir / "coverage_unique_pmids.csv", unique_rows, fields_unique)
    write_csv(args.output_dir / "coverage_summary.csv", summary_rows, fields_summary)
    write_csv(
        args.output_dir / "downloaded_not_extracted_gene_rows.csv",
        [row for row in gene_rows if row["coverage_status"] == "downloaded_not_extracted"],
        fields_gene,
    )
    write_csv(
        args.output_dir / "not_downloaded_unique_pmids.csv",
        [row for row in unique_rows if row["coverage_status"] == "not_downloaded"],
        fields_unique,
    )
    write_markdown(args.output_dir, summary_rows)

    print("Extraction coverage QC complete")
    print(f"  gene_rows: {len(gene_rows)}")
    print(f"  unique_pmids: {len(unique_rows)}")
    print(f"  output_dir: {args.output_dir}")


if __name__ == "__main__":
    main()
