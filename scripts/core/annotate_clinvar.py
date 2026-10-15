#!/usr/bin/env python3
"""Add compact ClinVar annotations to patient_mentions rows.

The script uses variant_c_gnomad as an exact HGVS query against NCBI ClinVar.
It keeps the legacy `clinvar` column for traceability, but also writes two
review-facing columns:
  - clinvar_pathogenicity: Pathogenic, Likely pathogenic, or
    Pathogenic/Likely pathogenic. Other classifications are left blank.
  - clinvar_id: ClinVar Variation ID only (= the number in the VCV accession,
    e.g. VCV000192252 -> 192252). The esummary `measure_id` is the AlleleID and
    is kept only in the legacy `clinvar` string as "AlleleID=".
"""

import argparse
import json
import re
import sqlite3
import time
from pathlib import Path

try:
    import requests
except ImportError:
    requests = None


NCBI_ESEARCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
NCBI_ESUMMARY = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
MISSING_VARIANTS = {"", "not specified", "n/a", "na", "none", "unknown"}
PATHOGENICITY_MAP = {
    "pathogenic": "Pathogenic",
    "likely pathogenic": "Likely pathogenic",
    "pathogenic/likely pathogenic": "Pathogenic/Likely pathogenic",
    "likely pathogenic/pathogenic": "Pathogenic/Likely pathogenic",
}
PATHOGENICITY_ORDER = {
    "Pathogenic": 1,
    "Likely pathogenic": 2,
    "Pathogenic/Likely pathogenic": 3,
}


def ensure_clinvar_columns(conn):
    columns = {row[1] for row in conn.execute("PRAGMA table_info(patient_mentions)")}
    for column in ["clinvar", "clinvar_pathogenicity", "clinvar_id"]:
        if column not in columns:
            conn.execute(f"ALTER TABLE patient_mentions ADD COLUMN {column} TEXT")


def usable_variant(value):
    if value is None:
        return ""
    value = str(value).strip()
    if value.lower() in MISSING_VARIANTS:
        return ""
    if not value.startswith("NM_") or ":c." not in value:
        return ""
    return value


def load_cache(cache_path):
    if not cache_path or not cache_path.exists():
        return {}
    return json.loads(cache_path.read_text(encoding="utf-8"))


def save_cache(cache_path, cache):
    if not cache_path:
        return
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    cache_path.write_text(json.dumps(cache, indent=2, ensure_ascii=False), encoding="utf-8")


def compact_unique(values, order=None):
    kept = []
    seen = set()
    for value in values:
        value = str(value or "").strip()
        if not value or value in seen:
            continue
        seen.add(value)
        kept.append(value)
    if order:
        kept.sort(key=lambda item: order.get(item, 999))
    return "; ".join(kept)


def parse_clinvar_string(value):
    value = str(value or "").strip()
    entry = {
        "status": "not_queried",
        "clinvar": "not_queried",
        "clinvar_pathogenicity": "",
        "clinvar_id": "",
        "classification_raw": "",
        "uids": [],
    }
    if not value:
        return entry

    entry["clinvar"] = value
    if value.startswith("listed:"):
        entry["status"] = "listed"
        # VCV accession 번호가 Variation ID다 (예전 문자열의 "VariationID="는 실제로 AlleleID였음)
        entry["uids"] = [
            str(int(match.group(1)))
            for match in re.finditer(r"VCV0*([0-9]+)(?:\.\d+)?", value)
        ]
        entry["clinvar_id"] = compact_unique(entry["uids"])
    elif value == "not_listed":
        entry["status"] = "not_listed"
    elif value == "not_queried":
        entry["status"] = "not_queried"
    elif value.startswith("query_error:"):
        entry["status"] = "query_error"
    else:
        entry["status"] = "unknown"
    return entry


def normalize_cache_entry(value):
    if isinstance(value, dict):
        entry = parse_clinvar_string(value.get("clinvar", ""))
        entry.update({k: value.get(k, entry.get(k, "")) for k in entry})
        if entry.get("status") == "listed":
            entry["clinvar_id"] = parse_clinvar_string(entry.get("clinvar", "")).get("clinvar_id", "")
        if isinstance(entry.get("uids"), str):
            entry["uids"] = [part for part in re.split(r"[,;\s]+", entry["uids"]) if part]
        entry["clinvar_pathogenicity"] = normalize_pathogenicity(
            entry.get("clinvar_pathogenicity", "")
        ) or str(entry.get("clinvar_pathogenicity", "")).strip()
        return entry
    return parse_clinvar_string(value)


def normalize_pathogenicity(description):
    text = re.sub(r"\s+", " ", str(description or "").strip())
    return PATHOGENICITY_MAP.get(text.lower(), "")


def summarize_summary_items(items, fallback=None):
    fallback = normalize_cache_entry(fallback or {})
    clinvar_parts = []
    variation_ids = []
    raw_classifications = []
    compact_pathogenicities = []

    for item in items:
        if item.get("error"):
            continue
        accession = item.get("accession_version") or item.get("accession") or ""
        allele_ids = [
            str(variation.get("measure_id"))
            for variation in item.get("variation_set", []) or []
            if variation.get("measure_id")
        ]
        if accession:
            match = re.match(r"VCV0*([0-9]+)", accession)
            if match:
                variation_ids.append(str(int(match.group(1))))
            suffix = f"; AlleleID={compact_unique(allele_ids)}" if allele_ids else ""
            clinvar_parts.append(f"{accession}{suffix}")

        classification = (
            item.get("germline_classification", {}) or {}
        ).get("description", "")
        if classification:
            raw_classifications.append(classification)
        pathogenicity = normalize_pathogenicity(classification)
        if pathogenicity:
            compact_pathogenicities.append(pathogenicity)

    clinvar = fallback.get("clinvar", "not_queried")
    if clinvar_parts:
        clinvar = "listed: " + " | ".join(clinvar_parts)

    return {
        "status": "listed" if clinvar.startswith("listed:") else fallback.get("status", "not_queried"),
        "clinvar": clinvar,
        "clinvar_pathogenicity": compact_unique(
            compact_pathogenicities,
            order=PATHOGENICITY_ORDER,
        ),
        "clinvar_id": compact_unique(variation_ids) or fallback.get("clinvar_id", ""),
        "classification_raw": compact_unique(raw_classifications),
        "uids": fallback.get("uids", []),
    }


def request_json_with_retries(session, url, params, attempts=4):
    last_error = None
    for attempt in range(1, attempts + 1):
        try:
            response = session.get(url, params=params, timeout=30)
            if response.status_code == 429:
                wait = max(2.0, attempt * 2.0)
                time.sleep(wait)
                continue
            response.raise_for_status()
            data = response.json()
            if "ERROR" in data.get("esearchresult", {}):
                raise RuntimeError(data["esearchresult"]["ERROR"].strip())
            return data
        except Exception as exc:
            last_error = exc
            if attempt < attempts:
                time.sleep(attempt * 1.5)
    raise last_error


def fetch_summaries_by_uids(session, uids, email="", tool="patient_phenotype_pipeline", batch_size=50):
    summaries = {}
    clean_uids = sorted({str(uid).strip() for uid in uids if str(uid).strip()}, key=lambda x: int(x))
    for start in range(0, len(clean_uids), batch_size):
        batch = clean_uids[start:start + batch_size]
        params = {
            "db": "clinvar",
            "id": ",".join(batch),
            "retmode": "json",
            "tool": tool,
        }
        if email:
            params["email"] = email
        data = request_json_with_retries(session, NCBI_ESUMMARY, params)
        result = data.get("result", {})
        for uid in result.get("uids", []):
            summaries[str(uid)] = result.get(str(uid), {})
        time.sleep(0.34)
    return summaries


def exact_items(items, variant):
    """ClinVar 텍스트 검색은 다른 변이(다른 위치·transcript 번호) 기록도 돌려준다.
    기록 제목의 transcript(버전 무시)와 c. 표기가 조회한 변이와 정확히 같은 단일 변이 기록만 남긴다."""
    m = re.match(r"(NM_\d+)(?:\.\d+)?:(c\..+)$", str(variant).strip())
    if not m:
        return items
    nm, cdna = m.group(1), m.group(2)
    same_nm, other_nm = [], []
    for item in items:
        title = str(item.get("title", ""))
        t = re.match(r"(NM_\d+)(?:\.\d+)?\([^)]*\):(c\.\S+)", title)
        if not t or t.group(2) != cdna or len(item.get("variation_set", []) or []) > 1:
            continue
        (same_nm if t.group(1) == nm else other_nm).append(item)
    # 같은 transcript 정확 일치를 우선, 없으면 제목만 다른 transcript로 적힌 같은 c. 기록 (예: STXBP1 NM_001032221)
    return same_nm or other_nm


def query_clinvar(variant, session, email="", tool="patient_phenotype_pipeline"):
    search_terms = [
        variant,
        f'"{variant}"[variant name]',
        f'"{variant}"',
    ]
    ids = []
    for term in search_terms:
        params = {
            "db": "clinvar",
            "term": term,
            "retmode": "json",
            "retmax": "10",
            "tool": tool,
        }
        if email:
            params["email"] = email
        search_json = request_json_with_retries(session, NCBI_ESEARCH, params)
        ids = search_json.get("esearchresult", {}).get("idlist", [])
        if ids:
            break

    if not ids:
        return {
            "status": "not_listed",
            "clinvar": "not_listed",
            "clinvar_pathogenicity": "",
            "clinvar_id": "",
            "classification_raw": "",
            "uids": [],
        }

    summaries = fetch_summaries_by_uids(session, ids, email=email, tool=tool)
    items = exact_items([summaries[uid] for uid in ids if uid in summaries], variant)
    if not items:
        return {
            "status": "not_listed",
            "clinvar": "not_listed (search hits were other variants)",
            "clinvar_pathogenicity": "",
            "clinvar_id": "",
            "classification_raw": "",
            "uids": [str(uid) for uid in ids],
        }
    entry = summarize_summary_items(items, fallback={"clinvar": f"listed: UID={','.join(ids)}"})
    entry["uids"] = [str(item.get("uid")) for item in items if item.get("uid")]
    return entry


def needs_detailed_refresh(entry):
    return (
        entry.get("status") == "listed"
        and not str(entry.get("classification_raw", "")).strip()
        and bool(entry.get("uids"))
    )


def annotate_db(db_path, cache, cache_path=None, delay_seconds=0.34, email=""):
    if requests is None:
        raise RuntimeError("requests package is required")

    conn = sqlite3.connect(db_path)
    try:
        ensure_clinvar_columns(conn)
        variants = [
            usable_variant(row[0])
            for row in conn.execute("""
                SELECT DISTINCT variant_c_gnomad
                FROM patient_mentions
                WHERE variant_c_gnomad IS NOT NULL
            """)
        ]
        variants = sorted({variant for variant in variants if variant})

        session = requests.Session()
        entries = {variant: normalize_cache_entry(cache.get(variant, "")) for variant in variants}

        uid_to_variants = {}
        for variant, entry in entries.items():
            if needs_detailed_refresh(entry):
                for uid in entry.get("uids", []):
                    uid_to_variants.setdefault(str(uid), []).append(variant)

        summary_queried = 0
        if uid_to_variants:
            summaries = fetch_summaries_by_uids(session, uid_to_variants, email=email)
            summary_queried = len(summaries)
            for variant, entry in list(entries.items()):
                if not needs_detailed_refresh(entry):
                    continue
                items = [
                    summaries[uid]
                    for uid in entry.get("uids", [])
                    if uid in summaries
                ]
                if items:
                    updated = summarize_summary_items(items, fallback=entry)
                    updated["uids"] = entry.get("uids", [])
                    entries[variant] = updated
                    cache[variant] = updated
            save_cache(cache_path, cache)

        searched = 0
        for variant, entry in list(entries.items()):
            should_search = (
                entry.get("status") in {"unknown", "query_error", "not_queried"}
                or not entry.get("clinvar")
            )
            if not should_search:
                continue
            try:
                entries[variant] = query_clinvar(variant, session, email=email)
            except Exception as exc:
                fallback = normalize_cache_entry(entry)
                fallback["status"] = "query_error"
                fallback["clinvar"] = fallback.get("clinvar") or f"query_error: {type(exc).__name__}"
                if not fallback["clinvar"].startswith("listed:"):
                    fallback["clinvar"] = f"query_error: {type(exc).__name__}"
                entries[variant] = fallback
            cache[variant] = entries[variant]
            searched += 1
            save_cache(cache_path, cache)
            time.sleep(delay_seconds)

        for variant, entry in entries.items():
            cache[variant] = entry
            conn.execute(
                """
                UPDATE patient_mentions
                SET clinvar = ?,
                    clinvar_pathogenicity = ?,
                    clinvar_id = ?
                WHERE variant_c_gnomad = ?
                """,
                (
                    entry.get("clinvar", "not_queried"),
                    entry.get("clinvar_pathogenicity", ""),
                    entry.get("clinvar_id", ""),
                    variant,
                ),
            )
        conn.commit()

        legacy_clinvar_values = [
            row[0]
            for row in conn.execute("""
                SELECT DISTINCT clinvar
                FROM patient_mentions
                WHERE clinvar LIKE 'listed:%'
                  AND (
                      clinvar_id IS NULL
                      OR trim(clinvar_id) = ''
                      OR clinvar_pathogenicity IS NULL
                      OR trim(clinvar_pathogenicity) = ''
                  )
            """)
            if row[0]
        ]
        legacy_entries = {
            clinvar: normalize_cache_entry(clinvar)
            for clinvar in legacy_clinvar_values
        }
        legacy_uids = []
        for entry in legacy_entries.values():
            legacy_uids.extend(entry.get("uids", []))
        legacy_summaries = fetch_summaries_by_uids(session, legacy_uids, email=email) if legacy_uids else {}
        for clinvar, entry in legacy_entries.items():
            items = [
                legacy_summaries[uid]
                for uid in entry.get("uids", [])
                if uid in legacy_summaries
            ]
            updated = summarize_summary_items(items, fallback=entry) if items else entry
            updated["uids"] = entry.get("uids", [])
            conn.execute(
                """
                UPDATE patient_mentions
                SET clinvar_pathogenicity = ?,
                    clinvar_id = ?
                WHERE clinvar = ?
                """,
                (
                    updated.get("clinvar_pathogenicity", ""),
                    updated.get("clinvar_id", ""),
                    clinvar,
                ),
            )
        conn.commit()

        conn.execute("""
            UPDATE patient_mentions
            SET clinvar = 'not_queried',
                clinvar_pathogenicity = '',
                clinvar_id = ''
            WHERE variant_c_gnomad IS NULL
               OR trim(variant_c_gnomad) = ''
               OR lower(trim(variant_c_gnomad)) IN ('not specified', 'n/a', 'na', 'none', 'unknown')
        """)
        conn.execute("""
            UPDATE patient_mentions
            SET clinvar = 'not_queried'
            WHERE clinvar IS NULL
               OR trim(clinvar) = ''
        """)
        conn.execute("""
            UPDATE patient_mentions
            SET clinvar_pathogenicity = ''
            WHERE clinvar_pathogenicity IS NULL
               OR trim(clinvar_pathogenicity) NOT IN (
                   'Pathogenic',
                   'Likely pathogenic',
                   'Pathogenic/Likely pathogenic'
               )
        """)
        conn.execute("""
            UPDATE patient_mentions
            SET clinvar_id = ''
            WHERE clinvar_id IS NULL
        """)
        conn.commit()

        listed = conn.execute(
            "SELECT count(*) FROM patient_mentions WHERE clinvar LIKE 'listed:%'"
        ).fetchone()[0]
        not_listed = conn.execute(
            "SELECT count(*) FROM patient_mentions WHERE clinvar = 'not_listed'"
        ).fetchone()[0]
        errors = conn.execute(
            "SELECT count(*) FROM patient_mentions WHERE clinvar LIKE 'query_error:%'"
        ).fetchone()[0]
        not_queried = conn.execute(
            "SELECT count(*) FROM patient_mentions WHERE clinvar = 'not_queried'"
        ).fetchone()[0]
        pathogenicity_rows = conn.execute("""
            SELECT count(*)
            FROM patient_mentions
            WHERE clinvar_pathogenicity IS NOT NULL
              AND trim(clinvar_pathogenicity) != ''
        """).fetchone()[0]
        clinvar_id_rows = conn.execute("""
            SELECT count(*)
            FROM patient_mentions
            WHERE clinvar_id IS NOT NULL
              AND trim(clinvar_id) != ''
        """).fetchone()[0]
        return {
            "db": str(db_path),
            "variants": len(variants),
            "summary_queried": summary_queried,
            "searched": searched,
            "listed_rows": listed,
            "not_listed_rows": not_listed,
            "not_queried_rows": not_queried,
            "query_error_rows": errors,
            "pathogenicity_rows": pathogenicity_rows,
            "clinvar_id_rows": clinvar_id_rows,
        }
    finally:
        conn.close()


def main():
    default_cache = (
        Path("databases") / "clinvar_cache.json"
        if Path("databases").exists()
        else Path("clinvar_cache.json")
    )
    parser = argparse.ArgumentParser(description="Annotate patient_mentions with compact ClinVar columns")
    parser.add_argument("db_paths", nargs="+", type=Path, help="SQLite DB path(s)")
    parser.add_argument(
        "--cache",
        type=Path,
        default=default_cache,
        help="JSON cache path to avoid repeated ClinVar queries",
    )
    parser.add_argument("--email", default="", help="Optional email for NCBI E-utilities")
    parser.add_argument("--delay", type=float, default=0.34, help="Delay between new NCBI searches")
    args = parser.parse_args()

    cache = load_cache(args.cache)
    for db_path in args.db_paths:
        summary = annotate_db(
            db_path,
            cache,
            cache_path=args.cache,
            delay_seconds=args.delay,
            email=args.email,
        )
        save_cache(args.cache, cache)
        print(
            f"{Path(db_path).stem}: variants={summary['variants']} "
            f"summary_queried={summary['summary_queried']} "
            f"searched={summary['searched']} listed_rows={summary['listed_rows']} "
            f"clinvar_id_rows={summary['clinvar_id_rows']} "
            f"pathogenicity_rows={summary['pathogenicity_rows']} "
            f"not_listed_rows={summary['not_listed_rows']} "
            f"not_queried_rows={summary['not_queried_rows']} "
            f"query_error_rows={summary['query_error_rows']}"
        )


if __name__ == "__main__":
    main()
