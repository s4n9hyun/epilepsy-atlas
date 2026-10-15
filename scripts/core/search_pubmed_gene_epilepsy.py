#!/usr/bin/env python3
"""Search PubMed for gene AND epilepsy article lists.

The output schema is compatible with scripts/core/import_article_lists_to_db.py.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import time
from collections import Counter
from datetime import datetime
from pathlib import Path

import requests


EMAIL = "researcher@snu.ac.kr"
BATCH_SIZE = 200
REQUEST_DELAY_SECONDS = 0.34
ESEARCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
ESUMMARY_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
IDCONV_URL = "https://www.ncbi.nlm.nih.gov/pmc/utils/idconv/v1.0/"
PMC_OA_URL = "https://www.ncbi.nlm.nih.gov/pmc/utils/oa/oa.fcgi"
PMC_ESEARCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
TOOL_NAME = "epilepsy_pubmed_gene_search"
OA_FILTER = "open_access[filter] NOT pmc embargo[filter]"
PDF_FILTER = "open_access[filter] AND has_pdf[filter] NOT pmc embargo[filter]"


def request_json(session: requests.Session, url: str, params: dict[str, str]) -> dict:
    last_error: Exception | None = None
    for attempt in range(1, 4):
        try:
            response = session.get(url, params=params, timeout=45)
            response.raise_for_status()
            return response.json()
        except Exception as exc:  # noqa: BLE001
            last_error = exc
            if attempt == 3:
                break
            time.sleep(1.5 * attempt)
    raise RuntimeError(f"Request failed after retries: {last_error}")


def request_text(session: requests.Session, url: str, params: dict[str, str]) -> str:
    last_error: Exception | None = None
    for attempt in range(1, 4):
        try:
            response = session.get(url, params=params, timeout=30)
            response.raise_for_status()
            return response.text
        except Exception as exc:  # noqa: BLE001
            last_error = exc
            if attempt == 3:
                break
            time.sleep(1.5 * attempt)
    raise RuntimeError(f"Request failed after retries: {last_error}")


def build_query(gene: str, mode: str) -> str:
    if mode == "epilepsy":
        return f"{gene} AND epilepsy"
    if mode == "epilepsy_or_seizure":
        return f"{gene} AND (epilepsy OR seizure)"
    raise ValueError(f"Unknown query mode: {mode}")


def pmcid_to_uid(pmcid: str) -> str | None:
    match = re.fullmatch(r"PMC(\d+)", (pmcid or "").strip())
    if not match:
        return None
    return match.group(1)


def query_pmc_subset(
    session: requests.Session,
    pmc_uids: list[str],
    filter_clause: str,
    email: str,
) -> set[str]:
    if not pmc_uids:
        return set()
    uid_terms = " OR ".join(f"{uid}[uid]" for uid in pmc_uids)
    term = f"({uid_terms}) AND {filter_clause}"
    data = request_json(
        session,
        PMC_ESEARCH_URL,
        {
            "db": "pmc",
            "term": term,
            "retmode": "json",
            "retmax": str(len(pmc_uids)),
            "email": email,
            "tool": TOOL_NAME,
        },
    )
    return set(data.get("esearchresult", {}).get("idlist", []))


def check_pmc_oa_batches(
    session: requests.Session,
    pmcids: list[str],
    cache: dict[str, dict[str, str]],
    email: str,
    chunk_size: int = 100,
) -> None:
    total = len(pmcids)
    for start in range(0, total, chunk_size):
        batch = pmcids[start : start + chunk_size]
        checked_at = datetime.now().isoformat(timespec="seconds")
        uid_by_pmcid = {pmcid: pmcid_to_uid(pmcid) for pmcid in batch}
        valid_uids = [uid for uid in uid_by_pmcid.values() if uid]
        try:
            oa_uids = query_pmc_subset(session, valid_uids, OA_FILTER, email)
            pdf_uids = query_pmc_subset(session, valid_uids, PDF_FILTER, email)
        except Exception as exc:  # noqa: BLE001
            detail = f"pmc_esearch_error={type(exc).__name__}: {exc}"
            for pmcid, uid in uid_by_pmcid.items():
                cache[pmcid] = {
                    "pmcid": pmcid,
                    "pmc_uid": uid or "",
                    "oa_status": "Error",
                    "full_package_available": "Error",
                    "has_pdf_in_pmc_oa": "Error",
                    "oa_checked_at": checked_at,
                    "oa_access_basis": detail,
                }
        else:
            for pmcid, uid in uid_by_pmcid.items():
                if not uid:
                    cache[pmcid] = {
                        "pmcid": pmcid,
                        "pmc_uid": "",
                        "oa_status": "Error",
                        "full_package_available": "Error",
                        "has_pdf_in_pmc_oa": "Error",
                        "oa_checked_at": checked_at,
                        "oa_access_basis": "invalid_pmcid_format",
                    }
                    continue
                is_oa = uid in oa_uids
                has_pdf = uid in pdf_uids
                cache[pmcid] = {
                    "pmcid": pmcid,
                    "pmc_uid": uid,
                    "oa_status": "OA" if is_oa else "Not OA",
                    "full_package_available": "Yes" if is_oa else "No",
                    "has_pdf_in_pmc_oa": "Yes" if has_pdf else "No",
                    "oa_checked_at": checked_at,
                    "oa_access_basis": (
                        "PMC open_access[filter], not pmc embargo[filter]"
                        if is_oa
                        else "PMCID present but not in PMC open_access[filter]"
                    ),
                }
        print(f"    PMC OA checked {min(start + chunk_size, total)}/{total}", flush=True)


def status_fields_for_pmcid(
    pmcid: str,
    cache: dict[str, dict[str, str]],
    checked_at: str,
) -> dict[str, str]:
    if not pmcid:
        return {
            "full_package_available": "No PMCID",
            "oa_status": "No PMCID",
            "has_pdf_in_pmc_oa": "No PMCID",
            "oa_checked_at": checked_at,
            "oa_access_basis": "No PMCID in PubMed search result",
        }
    cached = cache.get(pmcid)
    if not cached:
        return {
            "full_package_available": "Not checked",
            "oa_status": "Not checked",
            "has_pdf_in_pmc_oa": "Not checked",
            "oa_checked_at": checked_at,
            "oa_access_basis": "PMCID was not checked",
        }
    return {
        "full_package_available": cached.get("full_package_available", "Error"),
        "oa_status": cached.get("oa_status", "Error"),
        "has_pdf_in_pmc_oa": cached.get("has_pdf_in_pmc_oa", "Error"),
        "oa_checked_at": cached.get("oa_checked_at", checked_at),
        "oa_access_basis": cached.get("oa_access_basis", "PMC OA check error"),
    }


def esearch(session: requests.Session, query: str, email: str) -> tuple[int, list[str]]:
    params = {
        "db": "pubmed",
        "term": query,
        "retmax": "10000",
        "retmode": "json",
        "email": email,
    }
    data = request_json(session, ESEARCH_URL, params)
    result = data["esearchresult"]
    count = int(result["count"])
    pmids = result["idlist"]
    if count > len(pmids):
        raise RuntimeError(f"Search returned {count} hits but fetched {len(pmids)} PMIDs")
    time.sleep(REQUEST_DELAY_SECONDS)
    return count, pmids


def esummary(session: requests.Session, pmids: list[str], email: str) -> dict[str, dict[str, str]]:
    rows: dict[str, dict[str, str]] = {}
    for start in range(0, len(pmids), BATCH_SIZE):
        batch = pmids[start : start + BATCH_SIZE]
        params = {
            "db": "pubmed",
            "id": ",".join(batch),
            "retmode": "json",
            "email": email,
        }
        data = request_json(session, ESUMMARY_URL, params)
        result = data.get("result", {})
        for pmid in batch:
            item = result.get(pmid, {})
            pmcid = ""
            doi = ""
            for article_id in item.get("articleids", []):
                idtype = article_id.get("idtype")
                if idtype == "pmc":
                    pmcid = article_id.get("value", "")
                elif idtype == "doi":
                    doi = article_id.get("value", "")
            rows[pmid] = {
                "pmid": pmid,
                "pmcid": pmcid,
                "doi": doi,
                "title": item.get("title", ""),
                "pubdate": item.get("pubdate", ""),
                "source": item.get("source", ""),
            }
        print(f"    esummary {start + 1}-{start + len(batch)}/{len(pmids)}", flush=True)
        time.sleep(REQUEST_DELAY_SECONDS)
    return rows


def fill_missing_pmcids(
    session: requests.Session,
    info: dict[str, dict[str, str]],
    email: str,
) -> int:
    missing = [pmid for pmid, row in info.items() if not row.get("pmcid")]
    found = 0
    for start in range(0, len(missing), BATCH_SIZE):
        batch = missing[start : start + BATCH_SIZE]
        params = {
            "ids": ",".join(batch),
            "format": "json",
            "email": email,
        }
        try:
            data = request_json(session, IDCONV_URL, params)
        except RuntimeError as exc:
            print(f"    ID converter failed for batch {start + 1}: {exc}", flush=True)
            time.sleep(REQUEST_DELAY_SECONDS)
            continue
        for rec in data.get("records", []):
            pmid = rec.get("pmid", "")
            pmcid = rec.get("pmcid", "")
            if pmid in info and pmcid and not info[pmid].get("pmcid"):
                info[pmid]["pmcid"] = pmcid
                found += 1
        time.sleep(REQUEST_DELAY_SECONDS)
    return found


def check_oa_package(session: requests.Session, pmcid: str, email: str) -> str:
    if not pmcid:
        return "No PMCID"
    normalized = pmcid if pmcid.startswith("PMC") else f"PMC{pmcid}"
    try:
        response = session.get(PMC_OA_URL, params={"id": normalized, "email": email}, timeout=20)
    except requests.Timeout:
        return "Timeout"
    except Exception:  # noqa: BLE001
        return "Error"
    if response.status_code == 404:
        return "Not checked"
    if response.status_code != 200:
        return "Error"
    text = response.text
    if ".tar.gz" in text and "oa_package" in text:
        return "Yes"
    return "No"


def write_csv(path: Path, fieldnames: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8-sig") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def load_genes(args: argparse.Namespace) -> list[str]:
    genes: list[str] = []
    if args.gene_file:
        genes.extend(
            line.strip()
            for line in args.gene_file.read_text(encoding="utf-8").splitlines()
            if line.strip() and not line.strip().startswith("#")
        )
    genes.extend(args.genes or [])
    genes = [gene.strip().upper() for gene in genes if gene.strip()]
    duplicates = [gene for gene, count in Counter(genes).items() if count > 1]
    if duplicates:
        raise SystemExit(f"Duplicate genes: {', '.join(sorted(duplicates))}")
    if not genes:
        raise SystemExit("No genes provided.")
    return genes


def main() -> None:
    parser = argparse.ArgumentParser(description="Search PubMed for gene AND epilepsy article lists")
    parser.add_argument("--genes", nargs="*", help="Gene symbols to search")
    parser.add_argument("--gene-file", type=Path, help="Optional newline-separated gene list")
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument(
        "--query-mode",
        choices=["epilepsy", "epilepsy_or_seizure"],
        default="epilepsy",
    )
    parser.add_argument("--check-oa", action="store_true", help="Check PMC OA full-package availability")
    parser.add_argument("--email", default=EMAIL)
    args = parser.parse_args()

    genes = load_genes(args)
    search_date = datetime.now().isoformat(timespec="seconds")
    session = requests.Session()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    article_fields = [
        "gene",
        "pmid",
        "pmcid",
        "doi",
        "title",
        "pubdate",
        "source",
        "full_package_available",
        "oa_status",
        "has_pdf_in_pmc_oa",
        "oa_checked_at",
        "oa_access_basis",
        "downloaded",
        "update_date",
        "query",
        "search_date",
    ]
    summary_fields = [
        "gene",
        "total_articles",
        "with_pmcid",
        "full_package_yes",
        "full_package_no",
        "full_package_not_checked",
        "no_pmcid",
        "fp_rate_of_total",
        "fp_rate_of_pmcid",
        "oa_yes",
        "oa_no",
        "oa_no_pmcid",
        "oa_error",
        "oa_has_pdf_yes",
        "query",
    ]

    all_rows: list[dict[str, str]] = []
    summary_rows: list[dict[str, str]] = []
    oa_cache: dict[str, dict[str, str]] = {}

    print(f"PubMed search started: {search_date}", flush=True)
    print(f"Genes: {', '.join(genes)}", flush=True)
    print(f"Query mode: {args.query_mode}", flush=True)
    print(f"OA full-package check: {'on' if args.check_oa else 'off'}", flush=True)

    for idx, gene in enumerate(genes, 1):
        query = build_query(gene, args.query_mode)
        print(f"\n[{idx}/{len(genes)}] {query}", flush=True)
        count, pmids = esearch(session, query, args.email)
        print(f"  esearch: {count} PMIDs", flush=True)
        info = esummary(session, pmids, args.email)
        added = fill_missing_pmcids(session, info, args.email)
        if added:
            print(f"  ID converter added PMCIDs: {added}", flush=True)

        pmcids = sorted({row["pmcid"] for row in info.values() if row.get("pmcid")})
        if args.check_oa:
            to_check = [pmcid for pmcid in pmcids if pmcid not in oa_cache]
            print(f"  PMCIDs: {len(pmcids)} ({len(to_check)} new OA checks)", flush=True)
            check_pmc_oa_batches(session, to_check, oa_cache, args.email)
        else:
            print(f"  PMCIDs: {len(pmcids)} (OA checks skipped)", flush=True)

        rows: list[dict[str, str]] = []
        run_checked_at = datetime.now().isoformat(timespec="seconds")
        for pmid in pmids:
            item = info[pmid]
            pmcid = item.get("pmcid", "")
            status_fields = status_fields_for_pmcid(pmcid, oa_cache, run_checked_at)
            rows.append(
                {
                    "gene": gene,
                    "pmid": pmid,
                    "pmcid": pmcid,
                    "doi": item.get("doi", ""),
                    "title": item.get("title", ""),
                    "pubdate": item.get("pubdate", ""),
                    "source": item.get("source", ""),
                    **status_fields,
                    "downloaded": "No",
                    "update_date": search_date[:10],
                    "query": query,
                    "search_date": search_date,
                }
            )

        rows.sort(key=lambda row: int(row["pmid"]), reverse=True)
        write_csv(args.output_dir / f"{gene}_pubmed_epilepsy_articles.csv", article_fields, rows)
        all_rows.extend(rows)

        total = len(rows)
        with_pmcid = sum(1 for row in rows if row["pmcid"])
        fp_yes = sum(1 for row in rows if row["full_package_available"] == "Yes")
        fp_no = sum(1 for row in rows if row["full_package_available"] == "No")
        fp_not_checked = sum(
            1
            for row in rows
            if row["full_package_available"] not in {"Yes", "No", "No PMCID"}
        )
        no_pmcid = sum(1 for row in rows if row["full_package_available"] == "No PMCID")
        oa_error = sum(1 for row in rows if row.get("oa_status") == "Error")
        has_pdf_yes = sum(1 for row in rows if row.get("has_pdf_in_pmc_oa") == "Yes")
        summary_rows.append(
            {
                "gene": gene,
                "total_articles": str(total),
                "with_pmcid": str(with_pmcid),
                "full_package_yes": str(fp_yes),
                "full_package_no": str(fp_no),
                "full_package_not_checked": str(fp_not_checked),
                "no_pmcid": str(no_pmcid),
                "fp_rate_of_total": f"{fp_yes / total * 100:.1f}%" if total else "0.0%",
                "fp_rate_of_pmcid": f"{fp_yes / with_pmcid * 100:.1f}%" if with_pmcid else "0.0%",
                "oa_yes": str(fp_yes),
                "oa_no": str(fp_no),
                "oa_no_pmcid": str(no_pmcid),
                "oa_error": str(oa_error),
                "oa_has_pdf_yes": str(has_pdf_yes),
                "query": query,
            }
        )
        print(
            "  summary: "
            f"total={total}, pmcid={with_pmcid}, oa_yes={fp_yes}, "
            f"oa_no={fp_no}, oa_not_checked={fp_not_checked}, no_pmcid={no_pmcid}",
            flush=True,
        )

    all_rows.sort(key=lambda row: (row["gene"], -int(row["pmid"])))
    unique: dict[str, dict[str, str]] = {}
    for row in all_rows:
        pmid = row["pmid"]
        existing = unique.setdefault(
            pmid,
            {
                "pmid": pmid,
                "pmcid": row["pmcid"],
                "doi": row["doi"],
                "title": row["title"],
                "pubdate": row["pubdate"],
                "source": row["source"],
                "full_package_available": row["full_package_available"],
                "oa_status": row.get("oa_status", ""),
                "has_pdf_in_pmc_oa": row.get("has_pdf_in_pmc_oa", ""),
                "oa_checked_at": row.get("oa_checked_at", ""),
                "oa_access_basis": row.get("oa_access_basis", ""),
                "downloaded": "No",
                "genes": "",
                "gene_count": "0",
                "search_date": row["search_date"],
            },
        )
        genes_for_article = set(existing["genes"].split(";")) if existing["genes"] else set()
        genes_for_article.add(row["gene"])
        existing["genes"] = ";".join(sorted(genes_for_article))
        existing["gene_count"] = str(len(genes_for_article))
        if not existing["pmcid"] and row["pmcid"]:
            existing["pmcid"] = row["pmcid"]
        if existing["full_package_available"] != "Yes" and row["full_package_available"] == "Yes":
            existing["full_package_available"] = "Yes"
            existing["oa_status"] = row.get("oa_status", "")
            existing["has_pdf_in_pmc_oa"] = row.get("has_pdf_in_pmc_oa", "")
            existing["oa_checked_at"] = row.get("oa_checked_at", "")
            existing["oa_access_basis"] = row.get("oa_access_basis", "")

    unique_fields = [
        "pmid",
        "pmcid",
        "doi",
        "title",
        "pubdate",
        "source",
        "full_package_available",
        "oa_status",
        "has_pdf_in_pmc_oa",
        "oa_checked_at",
        "oa_access_basis",
        "downloaded",
        "genes",
        "gene_count",
        "search_date",
    ]
    unique_rows = sorted(unique.values(), key=lambda row: int(row["pmid"]), reverse=True)

    write_csv(args.output_dir / "seven_gene_article_rows.csv", article_fields, all_rows)
    write_csv(args.output_dir / "seven_gene_unique_articles.csv", unique_fields, unique_rows)
    write_csv(args.output_dir / "seven_gene_summary.csv", summary_fields, summary_rows)

    def count_by(rows: list[dict[str, str]], field: str) -> dict[str, int]:
        values = sorted({row.get(field, "") for row in rows})
        return {value: sum(1 for row in rows if row.get(field, "") == value) for value in values}

    (args.output_dir / "search_metadata.json").write_text(
        json.dumps(
            {
                "search_date": search_date,
                "genes": genes,
                "query_mode": args.query_mode,
                "gene_article_rows": len(all_rows),
                "unique_pmids": len(unique_rows),
                "oa_full_package_check": args.check_oa,
                "oa_check_method": "NCBI E-utilities PMC esearch with open_access[filter] and NOT pmc embargo[filter]",
                "oa_status_counts_unique": count_by(unique_rows, "oa_status"),
                "oa_has_pdf_counts_unique": count_by(unique_rows, "has_pdf_in_pmc_oa"),
                "per_gene_oa_yes": {
                    row["gene"]: int(row["full_package_yes"])
                    for row in summary_rows
                },
                "output_dir": str(args.output_dir),
            },
            indent=2,
            ensure_ascii=False,
        )
        + "\n",
        encoding="utf-8",
    )

    print("\nFINAL SUMMARY", flush=True)
    for row in summary_rows:
        print(
            f"  {row['gene']}: total={row['total_articles']}, "
            f"pmcid={row['with_pmcid']}, oa_yes={row['full_package_yes']}, "
            f"no_pmcid={row['no_pmcid']}",
            flush=True,
        )
    print(f"  gene_article_rows: {len(all_rows)}", flush=True)
    print(f"  unique_pmids: {len(unique_rows)}", flush=True)
    print(f"  output_dir: {args.output_dir}", flush=True)


if __name__ == "__main__":
    main()
