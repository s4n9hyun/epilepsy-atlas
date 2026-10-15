#!/usr/bin/env python3
"""
환자별 변이 + Phenotype 정보 추출 (Responses API)

논문 폴더의 main PDF와 supplementary 파일을 OpenAI Responses API에 함께 업로드해
유전자별 patient mention, variant, phenotype 정보를 JSON으로 추출한다.

추출 후에는 cDNA 표기를 gnomAD 기준 transcript로 보정하고 ClinVar ID/pathogenicity를 조회한 뒤,
결과를 유전자별 SQLite DB의 patient_mentions / processed_articles 테이블에 저장한다.
이 스크립트는 중복 제거를 하지 않고, 논문별 patient mention을 보존하는 데 집중한다.
"""

import atexit
import json
import re
import argparse
import logging
import time
import sqlite3
import csv
from pathlib import Path
import os
from datetime import datetime
import urllib.parse
import tempfile
import shutil
import zipfile
import html
import xml.etree.ElementTree as ET
from html.parser import HTMLParser

try:
    import requests
except ImportError:
    requests = None

try:
    from openai import OpenAI
except ImportError:
    OpenAI = None

try:
    from annotate_clinvar import (
        fetch_summaries_by_uids,
        load_cache,
        needs_detailed_refresh,
        normalize_cache_entry,
        query_clinvar,
        save_cache,
        summarize_summary_items,
        usable_variant,
    )
except ImportError:
    try:
        from scripts.core.annotate_clinvar import (
            fetch_summaries_by_uids,
            load_cache,
            needs_detailed_refresh,
            normalize_cache_entry,
            query_clinvar,
            save_cache,
            summarize_summary_items,
            usable_variant,
        )
    except ImportError:
        fetch_summaries_by_uids = load_cache = needs_detailed_refresh = normalize_cache_entry = None
        query_clinvar = save_cache = summarize_summary_items = usable_variant = None

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)


class InvalidApiKeyError(RuntimeError):
    """OpenAI API key가 잘못되어 더 진행할 수 없을 때 사용."""


class FatalApiQuotaError(RuntimeError):
    """OpenAI API quota가 없어 배치를 계속 진행할 수 없을 때 사용."""


class OutputTruncatedError(RuntimeError):
    """응답이 출력 한도에 걸려 잘렸을 때 사용 (잘린 JSON을 결과로 쓰지 않는다)."""


def is_quota_error(error):
    """quota 부족처럼 재시도/다음 논문 진행이 의미 없는 OpenAI 오류인지 확인한다."""
    message = str(error)
    return (
        "insufficient_quota" in message
        or "exceeded your current quota" in message
    )


def load_env_file(env_path):
    """간단한 KEY=VALUE 형식의 .env를 환경변수로 읽는다."""
    if not env_path.exists():
        return
    for line in env_path.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if not line or line.startswith("#") or "=" not in line:
            continue
        key, value = line.split("=", 1)
        os.environ.setdefault(key.strip(), value.strip().strip('"').strip("'"))


WORKSPACE_ROOT = Path(__file__).resolve().parents[2]

load_env_file(WORKSPACE_ROOT / ".env")
OPENAI_API_KEY = os.getenv("OPENAI_API_KEY")
TARGET_GENE = "EEF1A2"
RESPONSE_API_TIMEOUT_SECONDS = int(os.getenv("RESPONSE_API_TIMEOUT_SECONDS", "180"))
# 추출 모델: 2026-09-28부터 gpt-6-sol (이전 6개 유전자는 gpt-5.4). --model 또는 OPENAI_MODEL로 바꿀 수 있다.
EXTRACTION_MODEL = os.getenv("OPENAI_MODEL", "gpt-6-sol")
# 여러 유전자가 섞인 큰 코호트 표는 대상 유전자 행만 남겨 입력한다 (입력 한도·출력 크기 대비).
SPREADSHEET_FILTER_MIN_BYTES = 150 * 1024
SPREADSHEET_FILTER_MIN_ROWS = 300
SPREADSHEET_HEADER_ROWS = 5
LARGE_TEXT_FILTER_MIN_BYTES = 1024 * 1024
MAX_OUTPUT_CHUNKS = 20
CLINVAR_CACHE_PATH = WORKSPACE_ROOT / "databases" / "clinvar_cache.json"
GENE_REVIEW_GUIDANCE_PATH = WORKSPACE_ROOT / "config" / "gene_phenotype_review_guidance.csv"
NCBI_EMAIL = os.getenv("NCBI_EMAIL", "")

# gnomAD/ClinVar와 맞추기 위해 현재 MANE Select RefSeq transcript를 기준으로 쓴다.
# Mutalyzer 변환 결과를 어느 RefSeq transcript 기준으로 저장할지 정한다.
GNOMAD_TRANSCRIPTS = {
    "SCN1A":   "NM_001165963.4",
    "SCN2A":   "NM_001040142.2",
    "CDKL5":   "NM_001323289.2",
    "EEF1A2":  "NM_001958.5",
    "KCNQ2":   "NM_172107.4",
    "SYNGAP1": "NM_006772.3",
    "PCDH19":  "NM_001184880.2",  # MANE Select (2026-09-28 정정; 이전 NM_001105243.2)
    "SCN8A":   "NM_001330260.2",  # MANE Select (이전 NM_014191.4)
    "GABRA1":  "NM_001127644.2",  # MANE Select (이전 NM_000806.5, CDS 동일)
    "STXBP1":  "NM_001032221.6",  # MANE Select (이전 NM_003165.6)
    "SMC1A":   "NM_006306.4",
    "FOXG1":   "NM_005249.5",
    "KCNT1":   "NM_020822.3",
}

AA_1_TO_3 = {
    "A": "Ala",
    "R": "Arg",
    "N": "Asn",
    "D": "Asp",
    "C": "Cys",
    "Q": "Gln",
    "E": "Glu",
    "G": "Gly",
    "H": "His",
    "I": "Ile",
    "L": "Leu",
    "K": "Lys",
    "M": "Met",
    "F": "Phe",
    "P": "Pro",
    "S": "Ser",
    "T": "Thr",
    "W": "Trp",
    "Y": "Tyr",
    "V": "Val",
    "*": "Ter",
    "X": "Ter",
}

AA_3 = set(AA_1_TO_3.values()) - {"Ter"}

PDF_EXTENSIONS = {'.pdf'}
NXML_EXTENSIONS = {'.nxml', '.xml'}
HTML_EXTENSIONS = {'.html', '.htm'}

SUPPLEMENTARY_EXTENSIONS = {
    '.docx', '.doc', '.pptx', '.ppt',
    '.xlsx', '.xls', '.csv', '.tsv',
    '.txt', '.md'
}

SUPPORTED_EXTENSIONS = PDF_EXTENSIONS | SUPPLEMENTARY_EXTENSIONS | NXML_EXTENSIONS | HTML_EXTENSIONS
# ZIP 안 문서나 확장자 없이 받은 supplementary(downloadSupplement 등)를 magic bytes로 판별한다.
FILE_SIGNATURES = {b"%PDF": "pdf", b"PK\x03\x04": "zip", b"\xd0\xcf\x11\xe0": "ole"}
MAX_ZIP_MEMBERS = 40
# 저널이 ZIP에 함께 넣는 심사 기록은 환자 정보가 아니므로 제외한다.
PEER_REVIEW_FILE_PATTERN = re.compile(
    r"(author'?s? response|decision letter|original submission|revision \d|reviewer|peer.?review)", re.I
)
MAX_UPLOAD_FILE_BYTES = 32 * 1024 * 1024
MAX_NXML_TEXT_CHARS = int(os.getenv("MAX_NXML_TEXT_CHARS", "500000"))

KNOWN_GENE_SYMBOLS = {
    "SCN1A",
    "SCN2A",
    "CDKL5",
    "EEF1A2",
    "KCNQ2",
    "SYNGAP1",
    "PCDH19",
    "SCN8A",
    "GABRA1",
    "STXBP1",
    "SMC1A",
    "FOXG1",
    "KCNT1",
    "HNRNPU",
}

CONTEXT_RETRY_ERROR_SNIPPETS = (
    "context_length_exceeded",
    "Request too large",
    "tokens per min",
    "invalid_file",
    "badly formatted or corrupted",
    "timed out",
    "timeout",
    "ReadTimeout",
)


def invalid_downloaded_html_reason(file_path):
    """PMC direct download이 실제 파일 대신 HTML challenge를 저장한 경우를 감지한다."""
    path = Path(file_path)
    if path.suffix.lower() in {'.html', '.htm'}:
        return ""
    try:
        head = path.read_bytes()[:65536]
    except OSError:
        return ""
    text = head.decode('utf-8', errors='ignore')
    lowered = text.lower()
    if 'pow_challenge' in text or 'preparing to download' in lowered:
        return "PMC proof-of-work challenge HTML"
    if (
        'recaptcha/challengepage' in lowered
        or 'checking your browser before accessing pmc.ncbi.nlm.nih.gov' in lowered
        or 'g-recaptcha' in lowered
    ):
        return "PMC/Google reCAPTCHA browser-check HTML"
    return ""


def is_invalid_downloaded_html_file(file_path):
    return bool(invalid_downloaded_html_reason(file_path))


def clean_text(value):
    if value is None:
        return ""
    return str(value).strip()


def is_missing_text(value):
    return clean_text(value).lower() in {"", "not specified", "n/a", "na", "none", "unknown", "nan"}


def extract_nm_number(*values):
    """variant/transcript fields에서 RefSeq NM accession만 분리한다."""
    for value in values:
        text = clean_text(value)
        if is_missing_text(text):
            continue
        match = re.search(r"\bNM_\d+(?:\.\d+)?\b", text)
        if match:
            return match.group(0)
    return ""


def normalize_bool(value):
    if isinstance(value, bool):
        return int(value)
    text = clean_text(value).lower()
    return int(text in {"1", "true", "yes", "y", "manual_download_required"})


def canonical_dataset_url(url, doi=""):
    url = clean_text(url)
    doi = clean_text(doi)
    if doi.lower().startswith("doi:"):
        doi = doi[4:].strip()
    if doi.startswith("10.5061/dryad."):
        return f"https://datadryad.org/dataset/doi:{doi}"
    if url.lower().startswith("https://doi.org/10.5061/dryad."):
        dryad_doi = url.split("doi.org/", 1)[1]
        return f"https://datadryad.org/dataset/doi:{dryad_doi}"
    return url


def should_keep_external_dataset_link(link, target_gene=""):
    """외부 patient/variant dataset landing/download 링크가 아닌 false positive를 제외한다."""
    repository = clean_text(link.get("repository")).lower()
    source_type = clean_text(link.get("source_type")).lower()
    expected_content = clean_text(link.get("expected_content")).lower()
    evidence_location = clean_text(link.get("evidence_location")).lower()
    link_text = clean_text(link.get("link_text")).lower()
    url = clean_text(link.get("url")).lower()
    combined_text = " ".join([
        repository,
        source_type,
        expected_content,
        evidence_location,
        link_text,
        url,
    ])

    explicit_dataset_url = (
        "datadryad.org/dataset/doi:" in url
        or "doi.org/10.5061/dryad." in url
        or "figshare" in url
        or "doi.org/10.6084/m9.figshare" in url
        or "mendeley.com" in url
        or "doi.org/10.17632/" in url
        or "zenodo.org/record" in url
        or "zenodo.org/records" in url
        or "doi.org/10.5281/zenodo" in url
        or "ega-archive.org/studies/egas" in url
        or "study_id=phs" in url
        or "term=phs" in url
        or "bioproject/prjna" in url
        or "sra/prjna" in url
        or "github.com/helbig-lab/scn2a" in url
    )
    if not explicit_dataset_url:
        return False

    target_gene = clean_text(target_gene).upper()
    other_genes = KNOWN_GENE_SYMBOLS - {target_gene}
    if target_gene and any(gene.lower() in combined_text for gene in other_genes):
        return False

    if "code availability" in evidence_location or "code availability" in combined_text:
        return False
    if "#supplementary-material" in url:
        return False
    if "supplementary material" in evidence_location and "frontiersin.org/articles/" in url:
        return False
    if "recruitment" in combined_text or "contact page" in combined_text:
        return False
    if "varsome.com" in url:
        return False
    if "clinvar" in url or "genbank" in url or "arrayexpress" in url:
        return False
    if "geo/query" in url or "ncbi.nlm.nih.gov/geo" in url:
        return False
    if "lovd.nl" in url or "mdpi.com/article" in url:
        return False
    if "db.cngb.org/cnsa" in url or "base.sfari.org" in url:
        return False
    if "genomicsengland.co.uk" in url:
        return False
    if "nuccore" in url or "getentry.ddbj" in url:
        return False
    if "clirinx.com" in url or "mlink.cc" in url or "/dss" in url:
        return False
    if url.rstrip("/") in {"http://www.ebi.ac.uk/ega", "https://www.ebi.ac.uk/ega"}:
        return False
    if "doi.org/10." in url and repository in {"unknown", "", "publisher"}:
        return False

    non_patient_dataset_terms = (
        "rna-seq",
        "rnaseq",
        "single-cell",
        "single-nuclei",
        "transcriptomic",
        "expression data",
        "microarray",
        "hi-c",
        "contact maps",
        "atac-seq",
        "mouse",
        "mice",
        "animal",
        "pre-clinical",
        "cell line",
        "protein",
        "electrophysiology",
        "eeg-lfp",
    )
    if any(term in combined_text for term in non_patient_dataset_terms):
        return False

    if repository in {"github", "zenodo", "unknown", "", "publisher"}:
        content_has_patient_variant = (
            "patient" in expected_content
            and ("variant" in expected_content or "genotype" in expected_content)
        )
        if not content_has_patient_variant:
            return False

    if expected_content in {"unclear", ""}:
        return False

    if "doi.org/10." in url and ("citation" in evidence_location or "footer" in evidence_location):
        return False

    dataset_terms = (
        "patient",
        "variant",
        "genotype",
        "phenotype",
        "cohort",
        "sequencing",
        "supplement",
        "supplementary",
        "dataset",
        "data set",
        "dryad",
        "figshare",
        "zenodo",
        "dbgap",
        "ega",
    )
    return any(term in combined_text for term in dataset_terms)

def is_supplementary_file(filename):
    """파일명이 supplementary 자료인지 판별"""
    nl = filename.lower()
    supp_keywords = [
        'supplement', 'supp', 'esm', 'moesm', 'mmc',
        'data_sheet', 'checklist', 'appendix', 'supporting',
        'suppl', 'suplet', 'additional',
    ]
    if any(k in nl for k in supp_keywords):
        return True
    if re.search(r'[-_](s|sup|supp|esm|moesm)\d+', nl):
        return True
    return False


def file_signature(path):
    """앞 4바이트로 pdf / zip / ole(doc·xls·ppt) 여부를 판별한다."""
    try:
        with open(path, "rb") as f:
            head = f.read(4)
    except OSError:
        return None
    return next((kind for sig, kind in FILE_SIGNATURES.items() if head.startswith(sig)), None)


def ooxml_extension(names):
    """ZIP 목록이 Office Open XML 문서면 확장자(.docx/.xlsx/.pptx)를 반환한다."""
    if "[Content_Types].xml" not in names:
        return None
    for prefix, ext in (("word/", ".docx"), ("xl/", ".xlsx"), ("ppt/", ".pptx")):
        if any(n.startswith(prefix) for n in names):
            return ext
    return None


def ole_document_extension(data):
    """OLE 파일 안 stream 이름(UTF-16)으로 .xls/.doc/.ppt를 구분한다."""
    for stream, ext in (("Workbook", ".xls"), ("Book", ".xls"), ("WordDocument", ".doc"),
                        ("PowerPoint Document", ".ppt")):
        if stream.encode("utf-16-le") in data:
            return ext
    return None


def is_main_fulltext_html_file(filename):
    """PMC 본문 HTML처럼 논문 본문으로 쓸 수 있는 HTML만 허용한다."""
    nl = filename.lower()
    if nl.startswith("pubmed_"):
        return False
    return any(token in nl for token in ("pmc", "fulltext", "main", "article"))


def is_main_fulltext_text_file(filename):
    """ScienceDirect 등에서 정리 추출한 본문 TXT 파일을 main source로 허용한다."""
    nl = filename.lower()
    return nl in {"main.fulltext.txt", "main_article.txt", "article.fulltext.txt"}


def is_rendered_fulltext_pdf_file(filename):
    """로컬 HTML에서 렌더링한 보관용 PDF인지 판별한다."""
    nl = filename.lower()
    return nl in {"main.rendered.pdf", "main.fulltext.rendered.pdf"} or nl.endswith(".rendered.pdf")


class HTMLTextExtractor(HTMLParser):
    """간단한 본문 HTML 텍스트 추출기."""

    BLOCK_TAGS = {
        "address", "article", "aside", "blockquote", "br", "caption", "div",
        "figcaption", "figure", "h1", "h2", "h3", "h4", "h5", "h6", "li",
        "p", "section", "td", "th", "tr",
    }
    SKIP_TAGS = {"script", "style", "svg", "noscript"}

    def __init__(self):
        super().__init__(convert_charrefs=True)
        self.parts = []
        self._skip_depth = 0

    def handle_starttag(self, tag, attrs):
        tag = tag.lower()
        if tag in self.SKIP_TAGS:
            self._skip_depth += 1
            return
        if tag in self.BLOCK_TAGS:
            self.parts.append("\n")

    def handle_endtag(self, tag):
        tag = tag.lower()
        if tag in self.SKIP_TAGS and self._skip_depth:
            self._skip_depth -= 1
            return
        if tag in self.BLOCK_TAGS:
            self.parts.append("\n")

    def handle_data(self, data):
        if self._skip_depth:
            return
        if data and data.strip():
            self.parts.append(data.strip())
            self.parts.append(" ")

    def text(self):
        return " ".join("".join(self.parts).split())


def flatten_patients(result):
    """추출 결과 JSON을 patient mention row 목록으로 변환"""
    rows = []
    pmcid = result.get("pmcid", "")
    pmid = result.get("pmid", "")
    for patient in result.get("patients", []):
        phenotype = patient.get("phenotype", {}) or {}
        row = {
            "nm_number": extract_nm_number(
                patient.get("variant_c_gnomad", ""),
                patient.get("nm_transcript", ""),
                patient.get("gnomad_transcript", ""),
                patient.get("variant_c_original", ""),
            ),
            "variant_c_gnomad": patient.get("variant_c_gnomad", ""),
            "variant_p": normalize_variant_p(patient.get("variant_p", "")),
            "variant_type": patient.get("variant_type", ""),
            "clinvar": patient.get("clinvar", ""),
            "clinvar_pathogenicity": patient.get("clinvar_pathogenicity", ""),
            "clinvar_id": patient.get("clinvar_id", ""),
            "patient_id": patient.get("patient_id", ""),
            "pmid": pmid,
            "sex": patient.get("sex", ""),
            "korean_patient_source": patient.get("korean_patient_source", ""),
            "korean_patient_evidence": patient.get("korean_patient_evidence", ""),
            "seizure_presence": phenotype.get("seizure_presence", ""),
            "age_of_onset": phenotype.get("age_of_onset", ""),
            "last_followup_age": phenotype.get("last_followup_age", ""),
            "seizure_types": phenotype.get("seizure_types", ""),
            "phenotype_review_category": phenotype.get("phenotype_review_category", ""),
            "response_to_asm": phenotype.get("response_to_asm", ""),
            "asm_used": phenotype.get("asm_used", ""),
            "responsive_asm": phenotype.get("responsive_asm", ""),
            "targeted_treatment_response": phenotype.get("targeted_treatment_response", ""),
            "epilepsy_type": phenotype.get("epilepsy_type", ""),
            "other_phenotypes": phenotype.get("other_phenotypes", ""),
            "mention_source_type": patient.get("mention_source_type", ""),
            "originality_evidence": patient.get("originality_evidence", ""),
            "cited_source": patient.get("cited_source", ""),
            "location_in_paper": patient.get("location_in_paper", ""),
            "confidence": patient.get("confidence", ""),
            "used_supplementary": result.get("used_supplementary", ""),
            "files_processed": json.dumps(result.get("files_processed", []), ensure_ascii=False),
            "variant_c_original": patient.get("variant_c_original", ""),
            "gnomad_conversion_status": patient.get("gnomad_conversion_status", ""),
            "nm_transcript": patient.get("nm_transcript", ""),
            "gnomad_transcript": patient.get("gnomad_transcript", ""),
            "pmcid": pmcid,
            "raw_patient_json": json.dumps(patient, ensure_ascii=False),
            "extraction_timestamp": result.get("timestamp", ""),
            "variant_classification_in_paper": patient.get("variant_classification_in_paper", ""),
            "variant_evidence_source": patient.get("variant_evidence_source", ""),
            "variant_phase": patient.get("variant_phase", ""),
            "clinvar_classification": patient.get("clinvar_classification", ""),
        }
        if "mention_id" in patient:
            row["_mention_id"] = patient.get("mention_id")
        rows.append(row)
    return rows


def normalize_variant_p(variant_p):
    """protein variant를 HGVS 3-letter amino acid 표기로 가볍게 정규화한다."""
    if variant_p is None:
        return ""
    value = str(variant_p).strip()
    if not value:
        return ""
    if value in {"Not specified", "N/A", "Not Applicable", "p."}:
        return "Not specified"
    if value == "p.?":
        return "Not specified"
    if ";" in value:
        return "Not specified"

    # protein variant가 아닌 설명문은 유지하지 않고 명시되지 않은 값으로 처리한다.
    if not value.startswith("p") and not re.match(r"^[A-Z][a-z]{2}\d+", value) and not re.match(r"^[A-Z]\d+", value):
        return "Not specified"

    if value.startswith("p") and not value.startswith("p."):
        value = "p." + value[1:]
    elif not value.startswith("p."):
        value = "p." + value

    body = value[2:]
    if body.startswith("(") and body.endswith(")"):
        body = body[1:-1]
    if " (" in body:
        body = body.split(" (", 1)[0]
    body = body.replace(" ", "").replace("-FS", "fs").replace("-fs", "fs")
    body = re.sub(r"FS", "fs", body)

    def aa_token(token):
        if token in AA_3 or token == "Ter":
            return token
        if token in AA_1_TO_3:
            return AA_1_TO_3[token]
        return token

    aa_pattern = r"([A-Z][a-z]{2}|[A-Z]|\*|X)"

    # p.V306_W308del -> p.Val306_Trp308del
    m = re.fullmatch(rf"{aa_pattern}(\d+)_{aa_pattern}(\d+)(del.*)", body)
    if m:
        body = f"{aa_token(m.group(1))}{m.group(2)}_{aa_token(m.group(3))}{m.group(4)}{m.group(5)}"
        return f"p.{body}"

    # p.K114Sfs*20 / p.Val1195Alaf*27 / p.T488fs
    m = re.fullmatch(rf"{aa_pattern}(\d+)(?:{aa_pattern})?(fs|f)(?:\*|X|Ter)?(\d*)", body)
    if m:
        old_aa = aa_token(m.group(1))
        pos = m.group(2)
        new_aa = aa_token(m.group(3)) if m.group(3) else ""
        suffix = m.group(5)
        body = f"{old_aa}{pos}{new_aa}fsTer{suffix}" if suffix else f"{old_aa}{pos}{new_aa}fs"
        return f"p.{body}"

    # p.H95PfsXS처럼 stop 이후 글자가 깨진 경우: 위치/frameshift까지만 보존
    m = re.fullmatch(rf"{aa_pattern}(\d+)(?:{aa_pattern})?fsX[A-Z]+", body)
    if m:
        old_aa = aa_token(m.group(1))
        pos = m.group(2)
        new_aa = aa_token(m.group(3)) if m.group(3) else ""
        body = f"{old_aa}{pos}{new_aa}fs"
        return f"p.{body}"

    # p.R687X / p.Arg687* / p.E122K / p.Glu122Lys
    m = re.fullmatch(rf"{aa_pattern}(\d+){aa_pattern}", body)
    if m:
        body = f"{aa_token(m.group(1))}{m.group(2)}{aa_token(m.group(3))}"
        return f"p.{body}"

    # p.G391fs처럼 새 amino acid가 없는 frameshift
    m = re.fullmatch(rf"{aa_pattern}(\d+)fs", body)
    if m:
        body = f"{aa_token(m.group(1))}{m.group(2)}fs"
        return f"p.{body}"

    return f"p.{body}"


def normalize_to_gnomad(variant_c, gene_name, nm_source=None):
    """
    Mutalyzer API를 이용해 gnomAD 기준 transcript로 변이 표기 변환

    Args:
        variant_c:  논문에서 추출한 cDNA 표기 (e.g. "c.2528A>G" 또는 "NM_001165963.4:c.2528A>G")
        gene_name:  유전자명 (e.g. "SCN1A")
        nm_source:  논문에 명시된 NM 번호 (GPT가 별도 추출한 경우)

    Returns:
        (variant_c_gnomad, gnomad_conversion_status)
        status:
            'already_gnomad'      - gnomAD 기준 transcript 표기를 Mutalyzer로 확인
            'converted'           - Mutalyzer로 변환 성공
            'assumed_gnomad'      - NM 정보 없어 gnomAD 기준으로 가정하고 Mutalyzer로 확인
            'not_specified'       - 변이 정보 없음
            'unknown_gene'        - GNOMAD_TRANSCRIPTS에 없는 유전자
            'no_gnomad_found'     - Mutalyzer 응답에서 gnomAD 기준 transcript 미발견
            'api_error'           - API 오류
    """
    gnomad_transcript = GNOMAD_TRANSCRIPTS.get(gene_name.upper())
    if not gnomad_transcript:
        return variant_c, 'unknown_gene'

    if not variant_c or str(variant_c).strip() in ('Not specified', '', 'N/A'):
        return variant_c, 'not_specified'

    variant_c = str(variant_c).strip()

    # variant_c 자체에 NM prefix가 있는지 확인
    nm_match = re.match(r'(NM_\d+(?:\.\d+)?):(.+)', variant_c)
    if nm_match:
        nm_source = nm_match.group(1)
        c_part = nm_match.group(2)
    else:
        c_part = variant_c

    if not looks_like_cdna(c_part):
        return 'Not specified', 'invalid_hgvs'

    if requests is None:
        logging.warning("requests 패키지가 없어 Mutalyzer 변환을 건너뜁니다.")
        return f"{gnomad_transcript}:{c_part}", 'api_error'

    # NM 정보가 없거나 이미 기준 transcript인 경우도 Mutalyzer로 HGVS 문법을 확인한다.
    # c.1532-2_1532del처럼 c.로 시작하지만 유효하지 않은 표기를 걸러내기 위함이다.
    if not nm_source or nm_source == gnomad_transcript:
        query = f"{gnomad_transcript}:{c_part}"
        status = 'assumed_gnomad' if not nm_source else 'already_gnomad'
        try:
            encoded = urllib.parse.quote(query, safe='')
            url = f"https://mutalyzer.nl/api/normalize/{encoded}"
            r = requests.get(url, timeout=15)
            if r.status_code != 200:
                return query, 'api_error'
            return query, status
        except requests.Timeout:
            return query, 'api_error'
        except Exception as e:
            logging.warning(f"Mutalyzer 오류 ({query}): {e}")
            return query, 'api_error'

    # 다른 transcript로 쓰인 표기는 Mutalyzer API로 gnomAD 기준 transcript로 변환한다.
    query = f"{nm_source}:{c_part}"
    try:
        encoded = urllib.parse.quote(query, safe='')
        url = f"https://mutalyzer.nl/api/normalize/{encoded}"
        r = requests.get(url, timeout=15)

        if r.status_code != 200:
            return f"{gnomad_transcript}:{c_part}", 'api_error'

        data = r.json()
        gnomad_base = gnomad_transcript.split('.')[0]  # e.g. "NM_001165963"

        # equivalent_descriptions 에서 gnomAD 기준 NM 찾기
        for desc_item in data.get('equivalent_descriptions', {}).get('c', []):
            desc = desc_item.get('description', '')
            if desc.startswith(gnomad_transcript + ':'):
                return desc, 'converted'
            # 버전이 다르더라도 같은 base NM이면 반환
            if re.match(rf'{re.escape(gnomad_base)}\.\d+:', desc):
                return desc, 'converted'

        return f"{gnomad_transcript}:{c_part}", 'no_gnomad_found'

    except requests.Timeout:
        return f"{gnomad_transcript}:{c_part}", 'api_error'
    except Exception as e:
        logging.warning(f"Mutalyzer 오류 ({query}): {e}")
        return f"{gnomad_transcript}:{c_part}", 'api_error'


def looks_like_cdna(value):
    """기본적인 HGVS cDNA 모양인지 빠르게 판별한다."""
    if not value or str(value).strip() in {"", "Not specified", "N/A"}:
        return False
    return bool(re.match(r"^c\.[0-9*?+-]", str(value).strip()))


class PatientGeneDatabase:
    """유전자별 SQLite DB에 patient mentions를 append"""

    # variant/환자 정보가 먼저 보이고, 논문/파일 메타데이터는 뒤에 오도록 정렬한다.
    MENTION_COLUMNS = [
        "nm_number", "variant_c_gnomad", "variant_p",
        "variant_type",
        "clinvar", "clinvar_pathogenicity", "clinvar_id", "patient_id", "pmid",
        "sex", "korean_patient_source", "korean_patient_evidence",
        "seizure_presence", "age_of_onset",
        "last_followup_age", "seizure_types", "phenotype_review_category",
        "response_to_asm", "asm_used", "responsive_asm",
        "targeted_treatment_response", "epilepsy_type", "other_phenotypes",
        "mention_source_type", "originality_evidence",
        "cited_source", "location_in_paper",
        "confidence", "used_supplementary",
        "files_processed", "variant_c_original", "gnomad_conversion_status",
        "nm_transcript", "gnomad_transcript", "pmcid",
        "raw_patient_json", "extraction_timestamp",
        "variant_classification_in_paper", "variant_evidence_source", "variant_phase",
        "clinvar_classification",
    ]

    EXTERNAL_DATASET_COLUMNS = [
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
        "local_file_path",
        "notes",
        "discovered_at",
        "updated_at",
    ]

    def __init__(self, db_path):
        self.db_path = Path(db_path)
        self.db_path.parent.mkdir(parents=True, exist_ok=True)
        # 같은 DB에 여러 프로세스가 병렬로 쓸 수 있어 잠금 대기를 길게 둔다
        self.conn = sqlite3.connect(self.db_path, timeout=120)
        self.conn.row_factory = sqlite3.Row
        self._init_schema()

    def close(self):
        self.conn.close()

    @staticmethod
    def article_id(pmcid="", pmid=""):
        if pmid:
            return pmid
        if pmcid:
            return pmcid
        return ""

    def has_successful_article(self, pmcid="", pmid=""):
        article_id = self.article_id(pmcid, pmid)
        if not article_id:
            return False
        row = self.conn.execute("""
            SELECT status
            FROM processed_articles
            WHERE article_id = ?
        """, (article_id,)).fetchone()
        return row is not None and row["status"] == "success"

    def pmid_for_pmcid(self, pmcid):
        """pubmed_articles에서 PMCID에 대응하는 PMID를 찾는다."""
        table_exists = self.conn.execute("""
            SELECT 1
            FROM sqlite_master
            WHERE type = 'table' AND name = 'pubmed_articles'
        """).fetchone()
        if not table_exists:
            return ""
        row = self.conn.execute("""
            SELECT pmid
            FROM pubmed_articles
            WHERE pmcid = ?
        """, (pmcid,)).fetchone()
        return row["pmid"] if row and row["pmid"] else ""

    def _init_schema(self):
        # patient_mentions는 원본 추출 결과를 보존하는 append-only 성격의 테이블이다.
        # 주요 분석 컬럼을 먼저 두고, 내부 ID/논문 보조 메타데이터는 오른쪽으로 보낸다.
        self.conn.executescript("""
            CREATE TABLE IF NOT EXISTS patient_mentions (
                nm_number TEXT,
                variant_c_gnomad TEXT,
                variant_p TEXT,
                variant_type TEXT,
                clinvar TEXT,
                clinvar_pathogenicity TEXT,
                clinvar_id TEXT,
                patient_id TEXT,
                pmid TEXT,
                sex TEXT,
                korean_patient_source TEXT,
                korean_patient_evidence TEXT,
                seizure_presence TEXT,
                age_of_onset TEXT,
                last_followup_age TEXT,
                seizure_types TEXT,
                phenotype_review_category TEXT,
                response_to_asm TEXT,
                asm_used TEXT,
                responsive_asm TEXT,
                targeted_treatment_response TEXT,
                epilepsy_type TEXT,
                other_phenotypes TEXT,
                mention_source_type TEXT,
                originality_evidence TEXT,
                cited_source TEXT,
                location_in_paper TEXT,
                confidence TEXT,
                used_supplementary TEXT,
                files_processed TEXT,
                raw_patient_json TEXT,
                extraction_timestamp TEXT,
                variant_c_original TEXT,
                gnomad_conversion_status TEXT,
                nm_transcript TEXT,
                gnomad_transcript TEXT,
                pmcid TEXT,
                mention_id INTEGER PRIMARY KEY AUTOINCREMENT
            );

            CREATE TABLE IF NOT EXISTS processed_articles (
                article_id TEXT PRIMARY KEY,
                pmcid TEXT,
                pmid TEXT,
                gene TEXT,
                status TEXT,
                patient_count INTEGER,
                used_supplementary TEXT,
                files_processed TEXT,
                error TEXT,
                processed_at TEXT
            );

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

            CREATE INDEX IF NOT EXISTS idx_mentions_pmcid
                ON patient_mentions (pmcid);
            CREATE INDEX IF NOT EXISTS idx_mentions_source_type
                ON patient_mentions (mention_source_type);
            CREATE INDEX IF NOT EXISTS idx_article_files_article_id
                ON article_files (article_id);
            CREATE INDEX IF NOT EXISTS idx_article_files_pmid
                ON article_files (pmid);
            CREATE INDEX IF NOT EXISTS idx_external_dataset_links_pmid
                ON external_dataset_links (pmid);
            CREATE INDEX IF NOT EXISTS idx_external_dataset_links_pmcid
                ON external_dataset_links (pmcid);
            CREATE INDEX IF NOT EXISTS idx_external_dataset_links_manual
                ON external_dataset_links (manual_download_required);
            CREATE INDEX IF NOT EXISTS idx_external_dataset_links_extraction
                ON external_dataset_links (extraction_status);
        """)
        self._ensure_pubmed_article_columns()
        self._ensure_pipeline_columns()
        self.conn.commit()

    def _ensure_pipeline_columns(self):
        def ensure_column(table, column, column_sql):
            table_exists = self.conn.execute("""
                SELECT 1
                FROM sqlite_master
                WHERE type = 'table' AND name = ?
            """, (table,)).fetchone()
            if not table_exists:
                return
            existing_cols = {
                row[1] for row in self.conn.execute(f"PRAGMA table_info({table})").fetchall()
            }
            if column not in existing_cols:
                self.conn.execute(f"ALTER TABLE {table} ADD COLUMN {column_sql}")

        ensure_column("patient_mentions", "variant_c_original", "variant_c_original TEXT")
        ensure_column("patient_mentions", "nm_number", "nm_number TEXT")
        ensure_column("patient_mentions", "clinvar", "clinvar TEXT")
        ensure_column("patient_mentions", "clinvar_pathogenicity", "clinvar_pathogenicity TEXT")
        ensure_column("patient_mentions", "clinvar_id", "clinvar_id TEXT")
        ensure_column("patient_mentions", "variant_type", "variant_type TEXT")
        ensure_column("patient_mentions", "pmid", "pmid TEXT")
        ensure_column("patient_mentions", "korean_patient_source", "korean_patient_source TEXT")
        ensure_column("patient_mentions", "korean_patient_evidence", "korean_patient_evidence TEXT")
        ensure_column("patient_mentions", "mention_source_type", "mention_source_type TEXT")
        ensure_column("patient_mentions", "gnomad_transcript", "gnomad_transcript TEXT")
        ensure_column("patient_mentions", "gnomad_conversion_status", "gnomad_conversion_status TEXT")
        ensure_column("patient_mentions", "phenotype_review_category", "phenotype_review_category TEXT")
        ensure_column("patient_mentions", "targeted_treatment_response", "targeted_treatment_response TEXT")
        ensure_column("patient_mentions", "cited_source", "cited_source TEXT")
        ensure_column("patient_mentions", "variant_classification_in_paper", "variant_classification_in_paper TEXT")
        ensure_column("patient_mentions", "variant_evidence_source", "variant_evidence_source TEXT")
        ensure_column("patient_mentions", "variant_phase", "variant_phase TEXT")
        ensure_column("patient_mentions", "clinvar_classification", "clinvar_classification TEXT")
        ensure_column("processed_articles", "pmid", "pmid TEXT")
        ensure_column("processed_articles", "supplement_note", "supplement_note TEXT")

        self.conn.execute("""
            CREATE INDEX IF NOT EXISTS idx_mentions_pmid
            ON patient_mentions (pmid)
        """)

    def refresh_article_files(self, article_dir, pmcid="", pmid=""):
        """논문 폴더 안의 실제 파일 목록을 article_files 테이블에 저장한다."""
        article_path = Path(article_dir)
        article_id = pmid or pmcid
        if not article_id or not article_path.exists():
            return

        indexed_at = datetime.now().isoformat()
        files = sorted(
            path for path in article_path.rglob("*")
            if (
                path.is_file()
                and path.suffix.lower() in SUPPORTED_EXTENSIONS
                and not is_invalid_downloaded_html_file(path)
            )
        )

        with self.conn:
            self.conn.execute("""
                DELETE FROM article_files
                WHERE article_id = ?
            """, (article_id,))
            for path in files:
                suffix = path.suffix.lower()
                relative_path = path.relative_to(WORKSPACE_ROOT).as_posix()
                relative_article_path = path.relative_to(article_path).as_posix()
                self.conn.execute("""
                    INSERT INTO article_files (
                        article_id, pmid, pmcid, file_name, relative_path,
                        file_type, is_supplementary, file_size, indexed_at
                    )
                    VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
                """, (
                    article_id,
                    pmid,
                    pmcid,
                    path.name,
                    relative_path,
                    suffix.lstrip("."),
                    int(is_supplementary_file(relative_article_path)),
                    path.stat().st_size,
                    indexed_at,
                ))

    def _ensure_pubmed_article_columns(self):
        # article list를 먼저 import한 DB라면 추출 상태 컬럼만 가볍게 보강한다.
        table_exists = self.conn.execute("""
            SELECT 1
            FROM sqlite_master
            WHERE type = 'table' AND name = 'pubmed_articles'
        """).fetchone()
        if not table_exists:
            return

        existing_cols = {
            row[1] for row in self.conn.execute("PRAGMA table_info(pubmed_articles)").fetchall()
        }
        migrations = {
            "extracted": "ALTER TABLE pubmed_articles ADD COLUMN extracted TEXT NOT NULL DEFAULT 'No'",
            "extraction_status": "ALTER TABLE pubmed_articles ADD COLUMN extraction_status TEXT",
            "extracted_at": "ALTER TABLE pubmed_articles ADD COLUMN extracted_at TEXT",
        }
        for col, sql in migrations.items():
            if col not in existing_cols:
                self.conn.execute(sql)

        self.conn.execute("""
            CREATE INDEX IF NOT EXISTS idx_pubmed_articles_extracted
            ON pubmed_articles (extracted)
        """)

    def append_result(self, result):
        now = datetime.now().isoformat()
        rows = flatten_patients(result)
        article_id = self.article_id(result.get("pmcid", ""), result.get("pmid", ""))

        with self.conn:
            # 논문 단위 처리 상태를 먼저 기록해 재실행 시 skip 판단에 사용한다.
            self.conn.execute("""
                INSERT INTO processed_articles (
                    article_id, pmcid, pmid, gene, status, patient_count, used_supplementary,
                    files_processed, error, processed_at, supplement_note
                )
                VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                ON CONFLICT(article_id) DO UPDATE SET
                    pmcid=excluded.pmcid,
                    pmid=excluded.pmid,
                    gene=excluded.gene,
                    status=excluded.status,
                    patient_count=excluded.patient_count,
                    used_supplementary=excluded.used_supplementary,
                    files_processed=excluded.files_processed,
                    error=excluded.error,
                    processed_at=excluded.processed_at,
                    supplement_note=excluded.supplement_note
            """, (
                article_id,
                result.get("pmcid", ""),
                result.get("pmid", ""),
                result.get("gene", ""),
                result.get("status", ""),
                len(result.get("patients", [])),
                str(result.get("used_supplementary", "")),
                json.dumps(result.get("files_processed", []), ensure_ascii=False),
                result.get("error", ""),
                now,
                result.get("supplement_note", ""),
            ))
            self._update_article_extraction_status(
                result.get("pmcid", ""),
                result.get("pmid", ""),
                result.get("status", ""),
                now,
            )

            for row in rows:
                self._insert_mention(row)

            saved_links = self._upsert_external_dataset_links(result, now)
            if saved_links:
                logging.info(f"external_dataset_links 저장/갱신: {saved_links}개")

    def _insert_mention(self, row):
        columns = self.MENTION_COLUMNS
        placeholders = ", ".join(["?"] * len(columns))
        sql = f"""
            INSERT INTO patient_mentions ({", ".join(columns)})
            VALUES ({placeholders})
        """
        self.conn.execute(sql, [row.get(col, "") for col in columns])

    def _upsert_external_dataset_links(self, result, timestamp):
        pmid = clean_text(result.get("pmid"))
        pmcid = clean_text(result.get("pmcid"))
        gene = clean_text(result.get("gene"))
        rows = []

        for link in result.get("external_dataset_links", []) or []:
            if not isinstance(link, dict):
                continue
            if not should_keep_external_dataset_link(link, gene):
                continue

            url = clean_text(link.get("url"))
            doi = clean_text(link.get("doi"))
            if not url and doi:
                url = f"https://doi.org/{doi}"
            url = canonical_dataset_url(url, doi)
            if not url:
                continue

            access_status = clean_text(link.get("access_status")) or "unknown"
            manual_required = normalize_bool(link.get("manual_download_required"))
            if access_status in {"manual_download_needed", "token_required", "forbidden"}:
                manual_required = 1

            rows.append({
                "pmid": pmid,
                "pmcid": pmcid,
                "gene": gene,
                "source_type": clean_text(link.get("source_type")) or "external_dataset",
                "repository": clean_text(link.get("repository")),
                "url": url,
                "doi": doi,
                "link_text": clean_text(link.get("link_text")),
                "evidence_location": clean_text(link.get("evidence_location")),
                "expected_content": clean_text(link.get("expected_content")),
                "access_status": access_status,
                "manual_download_required": manual_required,
                "manual_download_reason": clean_text(link.get("manual_download_reason")),
                "extraction_status": clean_text(link.get("extraction_status")) or "not_processed",
                "local_file_path": clean_text(link.get("local_file_path")),
                "notes": clean_text(link.get("notes")),
                "discovered_at": timestamp,
                "updated_at": timestamp,
            })

        for row in rows:
            self.conn.execute("""
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
                    local_file_path = excluded.local_file_path,
                    notes = excluded.notes,
                    updated_at = excluded.updated_at
            """, row)

        return len(rows)

    def not_specified_mentions_for_article(self, pmcid="", pmid="", before_timestamp=None):
        """variant_c_gnomad와 variant_p가 모두 Not specified인 기존 row를 반환한다."""
        if pmid:
            where_sql = "pmid = ?"
            params = [pmid]
        elif pmcid:
            where_sql = "pmcid = ?"
            params = [pmcid]
        else:
            return []

        timestamp_sql = ""
        if before_timestamp:
            timestamp_sql = """
              AND (
                    extraction_timestamp IS NULL
                    OR TRIM(extraction_timestamp) = ''
                    OR extraction_timestamp < ?
              )
            """
            params.append(before_timestamp)

        sql = f"""
            SELECT *
            FROM patient_mentions
            WHERE {where_sql}
              AND COALESCE(TRIM(variant_c_gnomad), '') = 'Not specified'
              AND COALESCE(TRIM(variant_p), '') = 'Not specified'
              {timestamp_sql}
            ORDER BY mention_id
        """
        return [dict(row) for row in self.conn.execute(sql, tuple(params)).fetchall()]

    def replace_mentions(self, result):
        """재추출 결과를 기존 mention_id 행에 UPDATE한다. 새 행은 만들지 않는다."""
        rows = flatten_patients(result)
        updated = 0
        skipped = 0
        columns = self.MENTION_COLUMNS
        set_clause = ", ".join(f"{col} = ?" for col in columns)

        with self.conn:
            for row in rows:
                mention_id = row.get("_mention_id")
                if mention_id in (None, ""):
                    skipped += 1
                    continue
                try:
                    mention_id = int(mention_id)
                except (TypeError, ValueError):
                    skipped += 1
                    continue

                existing = self.conn.execute(
                    "SELECT mention_id FROM patient_mentions WHERE mention_id = ?",
                    (mention_id,),
                ).fetchone()
                if not existing:
                    skipped += 1
                    continue

                self.conn.execute(
                    f"UPDATE patient_mentions SET {set_clause} WHERE mention_id = ?",
                    [row.get(col, "") for col in columns] + [mention_id],
                )
                updated += 1

        return {"updated": updated, "skipped": skipped, "returned": len(rows)}

    def delete_article_mentions(self, pmcid="", pmid=""):
        """논문 단위 재추출 성공 후 기존 patient_mentions 행만 교체하기 위해 삭제한다."""
        if pmid:
            where_sql = "pmid = ?"
            params = (pmid,)
        elif pmcid:
            where_sql = "pmcid = ?"
            params = (pmcid,)
        else:
            return 0

        with self.conn:
            cursor = self.conn.execute(
                f"DELETE FROM patient_mentions WHERE {where_sql}",
                params,
            )
        return cursor.rowcount

    def _update_article_extraction_status(self, pmcid, pmid, status, timestamp):
        table_exists = self.conn.execute("""
            SELECT 1
            FROM sqlite_master
            WHERE type = 'table' AND name = 'pubmed_articles'
        """).fetchone()
        if not table_exists or not (pmcid or pmid):
            return

        extracted = "Yes" if status == "success" else "No"
        if pmcid:
            self.conn.execute("""
                UPDATE pubmed_articles
                SET extracted = ?,
                    extraction_status = ?,
                    extracted_at = ?
                WHERE pmcid = ?
            """, (extracted, status, timestamp, pmcid))
        else:
            self.conn.execute("""
                UPDATE pubmed_articles
                SET extracted = ?,
                    extraction_status = ?,
                    extracted_at = ?
                WHERE pmid = ?
            """, (extracted, status, timestamp, pmid))


class PhenotypeExtractor:
    """환자별 변이 + Phenotype 추출 클래스"""

    def __init__(self, api_key, gene_name):
        if OpenAI is None:
            raise ImportError("openai 패키지가 필요합니다. `pip install openai` 후 다시 실행하세요.")
        if not api_key:
            raise ValueError("OPENAI_API_KEY 환경변수가 필요합니다.")
        # 병렬 실행 시 속도 제한(429)·일시 오류를 SDK가 기다렸다 재시도 (이전 기본 0)
        self.client = OpenAI(api_key=api_key, max_retries=4)
        self.gene_name = gene_name
        self.clinvar_cache = load_cache(CLINVAR_CACHE_PATH) if load_cache else {}
        self.clinvar_session = requests.Session() if requests else None
        self.stats = {
            'total_articles': 0,
            'successful_articles': 0,
            'failed_articles': 0,
            'total_patients': 0,
            'external_links_found': 0,
            'api_errors': 0,
        }
        # 큰 환자 표는 행 범위로 나눠 여러 번 호출한다 (None이면 한 번에).
        self.chunk_size = None
        # 입력이 너무 크면 main 파일만으로 재시도할지. 재추출(supplement가 핵심)에서는 끈다.
        self.allow_main_only_retry = True
        # ZIP·확장자 없는 supplementary를 풀어 둔 임시 폴더 (다음 논문 처리 시·종료 시 삭제)
        self._packed_temp_dirs = []
        atexit.register(lambda: self._cleanup_temp_dirs(self._packed_temp_dirs))

    @staticmethod
    def _split_guidance_items(value):
        return [item.strip() for item in str(value or "").split("|") if item.strip()]

    @classmethod
    def load_gene_review_guidance(cls):
        """유전자별 phenotype review 항목을 config CSV에서 읽는다."""
        guidance = {}
        if not GENE_REVIEW_GUIDANCE_PATH.exists():
            return guidance
        with GENE_REVIEW_GUIDANCE_PATH.open(newline="", encoding="utf-8-sig") as f:
            for row in csv.DictReader(f):
                gene = clean_text(row.get("gene")).upper()
                if not gene:
                    continue
                guidance[gene] = {
                    "onset_unit": clean_text(row.get("onset_unit")),
                    "followup_unit": clean_text(row.get("followup_unit")),
                    "phenotype_categories": cls._split_guidance_items(row.get("phenotype_categories")),
                    "response_focus": cls._split_guidance_items(row.get("response_focus")),
                    "extraction_hints": cls._split_guidance_items(row.get("extraction_hints")),
                }
        return guidance

    @classmethod
    def get_gene_review_guidance(cls, target_gene):
        """유전자별 phenotype review 항목"""
        guidance = cls.load_gene_review_guidance()
        return guidance.get(target_gene.upper(), {
            "onset_unit": "Record seizure onset in the unit used by the paper; prefer months if possible.",
            "followup_unit": "Record last follow-up age in years when possible.",
            "phenotype_categories": [
                "Epilepsy",
                "Other epilepsy",
                "Non-seizure / seizure never reported",
            ],
            "response_focus": [
                "Any reported anti-seizure medication or treatment response",
            ],
            "extraction_hints": [
                "Do not infer treatment response unless explicitly stated.",
            ],
        })

    def create_extraction_prompt(self, target_gene):
        """환자별 변이 + phenotype 추출 프롬프트"""

        gene_guidance = self.get_gene_review_guidance(target_gene)
        # 유전자별 관심 phenotype/ASM 항목만 프롬프트에 끼워 넣는다.
        phenotype_categories = "\n".join(f"- {item}" for item in gene_guidance["phenotype_categories"])
        response_focus = "\n".join(f"- {item}" for item in gene_guidance["response_focus"])
        extraction_hints = "\n".join(f"- {item}" for item in gene_guidance["extraction_hints"])

        return f"""You are an expert clinical geneticist extracting structured patient mentions from papers.
Analyze ALL attached files and extract patient-level information for the **{target_gene}** gene.

Work conservatively:
- Prefer exact text from tables, case descriptions, supplementary files, and figure legends.
- Do not fill missing values from general medical knowledge.
- When information conflicts, prefer the most patient-specific source, then tables, then main text summaries.
- Keep one row per patient-variant mention. Do not merge different patients, families, or cited sources into one row.

## TASK:
Extract information for EACH PATIENT with a {target_gene} variant. Each patient should have:
1. Genetic variant information (variant_c_original, variant_p)
2. Clinical phenotype information (seizures, age of onset, response to treatment, etc.)
3. External dataset links, if the paper points readers to a repository or dataset outside the attached files.

## OUTPUT FORMAT (JSON only):

{{
  "patients": [
    {{
      "patient_id": "Patient 1 / Case 1 / ID from paper",
	      "mention_source_type": "original_patient / cited_patient / aggregate_only / variant_only / previously_reported_patient / reanalyzed_previous_patient / unclear",
	      "originality_evidence": "our patient / this study / previously reported by ... / table adapted from ... / Not specified",
	      "cited_source": "If cited, original paper author/year/PMID/PMCID/DOI; otherwise 'Not specified'",
	      "sex": "Male / Female / Not specified",
	      "korean_patient_source": "Korean multicenter cohort / Korean patient from literature / Not Korean / Not specified",
	      "korean_patient_evidence": "Evidence for Korean source classification, or 'Not specified'",
	      "nm_transcript": "NM_001165963.4 or 'Not specified' — the transcript/NM accession number stated in the paper",
	      "variant_c_original": "HGVS cDNA notation from the paper, without NM prefix; c.123A>G or 'Not specified'",
	      "variant_p": "HGVS protein notation using 3-letter amino acid names, e.g. p.Arg41His or 'Not specified'",
	      "variant_type": "sequence_variant_or_small_indel / exon_level_cnv / large_cnv / unclear",
	      "variant_classification_in_paper": "pathogenic / likely pathogenic / VUS / likely benign / benign / polymorphism / not stated",
	      "variant_evidence_source": "table / supplementary_table / main_text / figure_only",
	      "variant_phase": "single / cis / trans / unknown_phase",
      "phenotype": {{
        "seizure_presence": "Yes / No / Not specified",
        "age_of_onset": "6 months / 5 days / Not specified",
        "last_followup_age": "3 years / 12 years / Not specified",
        "seizure_types": "Focal motor / Epileptic spasms / Generalized tonic-clonic / Not specified",
        "phenotype_review_category": "One of the gene-specific phenotype categories below, or 'Not specified'",
        "response_to_asm": "Controlled / Refractory / Unknown / Not specified",
        "asm_used": "VPA, LEV / Not specified",
        "responsive_asm": "VPA / Not specified",
        "targeted_treatment_response": "Gene-specific medication/treatment response requested below, or 'Not specified'",
        "epilepsy_type": "DEE / Focal epilepsy / Not specified",
        "other_phenotypes": "Developmental delay, intellectual disability, etc."
      }},
	      "location_in_paper": "Table 1, row 3 / Case report section / Supplementary Table S1",
	      "confidence": "High / Medium / Low"
	    }}
	  ],
	  "missing_supplementary_note": "If the main text points to a patient-level supplementary table/file that is NOT among the attached files, name it here (e.g. 'Table S2 with individual genotypes not attached'); otherwise 'Not specified'",
	  "external_dataset_links": [
	    {{
	      "source_type": "external_dataset / controlled_access_repository / registry / publisher_external_data / other",
	      "repository": "Dryad / Figshare / Zenodo / dbGaP / EGA / publisher / unknown",
	      "url": "full URL if available",
	      "doi": "dataset DOI if available, otherwise empty string",
	      "link_text": "the visible link text or dataset title",
	      "evidence_location": "Data availability / Key resources table / Methods / Table note / page or file name",
	      "expected_content": "patient-level genotype/phenotype table / variants / sequencing data / unclear",
	      "access_status": "open / controlled_access / manual_download_needed / token_required / forbidden / unknown",
	      "manual_download_required": true,
	      "manual_download_reason": "why this should be manually downloaded or empty if not needed",
	      "extraction_status": "not_processed",
	      "local_file_path": "",
	      "notes": "short conservative note"
	    }}
	  ]
	}}

## IMPORTANT NOTES:
- Extract patient-level mentions related to {target_gene} from the attached files.
- Include both this paper's own patients and clearly cited/review patients, but label `mention_source_type` accurately.
- Each row should represent one patient-variant mention.
- If phenotype information is not available for a specific field, use "Not specified"
- Do not infer missing clinical details from syndrome names, gene knowledge, or typical disease course. Extract only what is explicitly stated in the attached files.
- Only combine variant and phenotype information when they clearly refer to the same patient ID, case number, family member, or table row.
- Do not create patient-level clinical rows from a paper that lists individual variants but gives phenotype only as cohort-level summary statistics. In that case, use `mention_source_type="variant_only"` for individual variant rows and keep patient-level phenotype fields as "Not specified".
- If variant-phenotype linkage is uncertain, use `mention_source_type="variant_only"` or `unclear`, set confidence to "Low", and describe the uncertainty in `location_in_paper` or `originality_evidence`.
- For **mention_source_type**:
  * Use "original_patient" only when the patient is newly reported by this paper or clearly belongs to this paper's own cohort/case series.
  * Use "cited_patient" when the patient is copied, summarized, re-tabulated, or discussed from prior literature.
  * Use "previously_reported_patient" when the paper includes a patient/cohort explicitly described as previously reported, previously described, previously published, reported by another paper, or reused from a prior cohort.
  * Use "reanalyzed_previous_patient" when the paper re-analyzes an already reported patient/cohort.
  * Use "variant_only" when individual genotype/variant rows are present but patient-level phenotype cannot be directly matched to that individual.
  * Do not mark a patient as "original_patient" unless the paper clearly presents the patient as part of its own cohort, case report, or newly analyzed samples.
  * Use "aggregate_only" when the paper reports only grouped statistics without individual patient data.
  * Use "unclear" when originality cannot be determined.
- Before assigning "original_patient", scan Methods, Participants, Data availability, table notes, and Supplementary files for cohort-source statements. If a mixed cohort says, for example, that some patients were previously described and only a subset is new, label only the explicitly new subset as "original_patient".
- For **originality_evidence**: quote or summarize the evidence that supports the role, such as "our patient", "this study", "previously reported by X et al.", or "adapted from".
- For **cited_source**: if the patient is cited from prior literature, record the cited source if available.
- For **korean_patient_source**: use "Korean multicenter cohort" only for the project's Korean cohort or an explicitly described Korean multicenter cohort; use "Korean patient from literature" for Korean patients reported in publications; otherwise use "Not Korean / Not specified".
- For **variant_type**:
  * Use "sequence_variant_or_small_indel" for SNV, missense, nonsense, splice-site, frameshift, or small insertion/deletion variants.
  * Use "exon_level_cnv" for exon-level deletion/duplication or CNV affecting one or more named exons.
  * Use "large_cnv" for whole-gene, multi-gene, chromosomal, microdeletion, microduplication, or other large structural CNV.
  * Use "unclear" when the variant class cannot be confidently determined from the paper.
- For **nm_transcript**: record the RefSeq transcript accession (e.g. NM_001165963.4) if stated anywhere in the paper (Methods, Table header, variant notation). Use "Not specified" if not mentioned.
- For **variant_c_original**: record only the cDNA part WITHOUT the NM prefix (e.g. "c.2528A>G"), even if the paper writes "NM_001165963.4:c.2528A>G".
- Use standard HGVS cDNA formatting for **variant_c_original** when the paper uses an obvious nonstandard arrow notation: write "c.208G>A", not "c.208G->A"; write "c.364G>A", not "c364G->A".
- Use valid HGVS separators and spelling for **variant_c_original**: no spaces inside the variant, use "_" for ranges (not "-"), and use "del" (not OCR errors such as "det").
- Prefer valid cDNA HGVS in **variant_c_original**. If the paper gives a nonstandard cDNA-like, genomic-coordinate, CNV, exon-level deletion, or free-text deletion genotype that cannot be normalized confidently, preserve that exact genotype text in **variant_c_original** rather than dropping the genotype evidence.
- Do not reinterpret or invent variants. Preserve nonstandard genotype text only when it is explicitly stated for that patient.
- For **variant_p**: normalize protein notation to one HGVS-style 3-letter amino acid value when clear. Examples: p.E122K -> p.Glu122Lys; p.R687X, p.Arg687*, or p.(R687*) -> p.Arg687Ter. Do not keep uncertainty parentheses in the saved value. If a splice variant has multiple predicted transcript/protein products, use "Not specified" instead of listing multiple `p.` values. Use "Not specified" if protein consequence is not stated or cannot be read confidently.
- Do NOT output `variant_c_gnomad`, `gnomad_transcript`, or `gnomad_conversion_status`; those are added later by the script.
- For **sex**: use only "Male", "Female", or "Not specified"
- For **last_followup_age**: record the latest age at which this patient's status is described, keeping the paper's precision and unit (e.g. "6 years 8 months", "11 months", "3 days"). See the ages rule below.
- For **location_in_paper**: include the most precise source available, such as table number, row label, supplementary filename, sheet name, page, figure, or case section.
- Use `confidence` as:
  * "High" when patient, variant, and phenotype are directly linked in the same row/case.
  * "Medium" when the link is clear but spread across nearby text/tables.
  * "Low" when the link is plausible but uncertain.
- Look for information in:
  * Main text (case reports, patient descriptions)
  * Tables (patient summary tables, variant tables)
  * Supplementary files (all attached files; files unpacked from a ZIP are named "<zip>__<file>")
  * Figures (pedigrees, clinical timelines) — for phenotype context only; see the figure rule below

## Additional rules (from manual review against the clinicians' curated database):
- **First genotype report**: if a patient's clinical features were published earlier (e.g. a linkage or clinical-series paper) but THIS paper is the first to report that patient's {target_gene} variant, label the patient "original_patient", put the earlier clinical paper in `cited_source`, and write "first genotype report" plus the earlier citation in `originality_evidence`. If an earlier paper already reported the same patient's variant, use "previously_reported_patient".
- **Every row of a patient table**: when a table or supplementary table lists individual patients, extract every row for {target_gene}, including long tables that continue across pages or sheets. Do not stop after the first rows and do not summarize a long patient table as aggregate_only.
- **Review papers that also add new patients**: extract each new patient (own ID, "not previously reported", "novel", "this study", no citation) as an individual "original_patient" row, and extract literature patients separately as "previously_reported_patient"/"cited_patient" with the original citation in `cited_source`. Never collapse such a paper into one aggregate_only row.
- **Missing supplementary files**: if the text says patient-level details are in a supplementary table that is not attached, still extract the individual patients that are available in the attached files, fill `missing_supplementary_note`, and do not use aggregate_only only because the table is missing.
- **Figures**: do not read variant notations from figures (protein schematics, lollipop plots, pedigree images) when a table or text lists the same patients. If a patient's variant appears only in a figure that links it to an identifiable individual (pedigree member, case photo/timeline, labelled patient), extract it with `variant_evidence_source="figure_only"` and confidence "Low". Do not create rows for variants that appear only on a protein schematic, lollipop plot or mutation map without an individual patient.
- **patient_id**: copy the paper's own patient/case/family-member label exactly (e.g. "Patient 12", "Family 3 II-1"). Do not add suffixes such as "- variant 1" to separate variants; one patient with two {target_gene} variants keeps the same patient_id in two rows.
- **Two or more variants in one patient**: output one row per variant (never join two variants into one string) and set `variant_phase` to "cis" when the paper states the variants are on the same allele/complex allele, "trans" when on different alleles, "unknown_phase" when not stated, and "single" when the patient has only one {target_gene} variant.
- **Target gene only**: `variant_c_original`/`variant_p` must be the {target_gene} variant. Do not put a variant from another gene (e.g. SCN1A in an SCN2A row, CNTNAP2 in a CDKL5 row) into these fields; mention a co-occurring variant in another gene in `other_phenotypes` only as "Additional variant: GENE c.xxx".
- **Paper's variant classification**: record the paper's own classification of the {target_gene} variant in `variant_classification_in_paper` (ACMG class, "polymorphism", "benign", "non-pathogenic", "VUS"); use "not stated" when the paper gives none. Still extract VUS, benign, and polymorphic variants that are listed for a patient, but label them accurately.
- **Nonstandard genotype text**: if a patient-specific genotype is written in any nonstandard form (exon deletion, "del exon 1-4", IVS notation, "c.?", breakpoints with "?", genomic coordinates, karyotype/array notation), copy that text exactly into `variant_c_original` instead of "Not specified".
- **Ages and units**: keep the unit written in the paper (days, weeks, months, years). Do not round days/weeks to 0 or 1 month. `age_of_onset` is age at first seizure (or first symptom only if the paper defines onset that way); never use age at examination, diagnosis, enrolment or genetic testing as onset. `last_followup_age` is the LATEST age at which the patient is described: age at last follow-up/last evaluation, current age ("currently 21 months", "a 3-year-old boy" in a case report), age at re-examination, age at enrolment or examination when no later age is given, or age at death (write e.g. "1 month 7 days (died)"). Do not leave it "Not specified" when any such age is stated for the patient, and do not put that age only into `other_phenotypes`.
- **No cohort values in patient rows**: never copy cohort medians/ranges, percentages, or general statements from the Introduction/Discussion into an individual patient's fields.
- **Explicit negatives**: when a table cell for this patient says "None", "No", "-" with a legend meaning absent, or "absent", write the negative finding descriptively (e.g. "No endocrine disease", "No movement disorder") instead of "Not specified"; use "Not specified" only when the information is missing.
- **Table rows**: in tables with merged cells or multi-line rows, keep each value with its own patient row; re-check that variant, sex, and age come from the same row before output.

## Additional rules (from the full source check of 2026-09-27, where these were the most common errors):
- **Relatives**: do not create rows for relatives who were not genotyped (pedigree members shown as affected but not tested), and do not create rows for genotyped carriers who have no seizures and no neurological or developmental phenotype (unaffected carriers, including unaffected mosaic parents). Mention them in the proband's `other_phenotypes` only if relevant (e.g. "Unaffected carrier father").
- **Originality wording**: start `originality_evidence` with exactly one of: "First reported in this paper:", "First genotype report:", "Previously reported in <citation>:", "Cited from <citation>:", followed by the supporting quote or location. For a mixed cohort, decide this per patient.
- **Same patient in several papers from one group**: if the paper says the cohort overlaps an earlier publication (same centre, "previously described", "an extended cohort"), label the overlapping patients "previously_reported_patient" with that citation, even when individual IDs changed.
- **Legacy transcript numbering**: copy the cDNA and protein exactly as the paper writes them and record the transcript in `nm_transcript` if stated. Do not renumber variants yourself; renumbering to MANE is done later.
- **Large deletions/duplications**: in `variant_c_original`, keep the interval as written and list every other gene the paper says is inside the interval (e.g. "20q13.33 del 0.72 Mb incl. KCNQ2, CHRNA4, EEF1A2"), so that multigene CNVs can be recognised.
- **Not patients**: never create rows for animal models, cell lines without a described patient, in-vitro constructs, population/control variants, or preprint duplicates of this paper.
- **Individual response to treatment**: record a drug response only when it is stated for this patient. A sentence about a group ("5 of 13 worsened on carbamazepine") is cohort-level and must not be copied into each patient; "no effect" is not "worsening"; improvement on a combination is not the effect of one drug.
- **Column alignment**: in wide tables, check each value against the patient's own column/row header before output; a value shifted from the neighbouring patient is the most common table error.
- **Seizure fields must agree**: if `phenotype_review_category` says seizures were never reported or absent, set `seizure_presence` to "No" only when the paper states no seizures, otherwise "Not specified"; if the category is an epilepsy category, `seizure_presence` must be "Yes".
- **Cited patients need the target-gene variant**: create a cited/previously reported row only when this paper states the patient's {target_gene} variant (or clearly states the patient carries a {target_gene} variant). Do not infer the gene from the title of a cited reference.
- **Main text vs supplementary table**: when they disagree about a patient's variant or clinical value, prefer the patient-level table and note the conflict in `location_in_paper`.

## External dataset link extraction:
- Also scan the article for URLs, DOIs, accession links, or repository pages that may contain patient-level genotype, variant, phenotype, cohort, or sequencing data for this paper.
- These links are usually in Data availability, Supplementary Data, Methods, Results, Key resources table, footnotes, or table notes.
- Include dataset/repository links such as Dryad, Figshare, Zenodo, Mendeley Data, OSF, dbGaP, EGA, SRA/BioProject, or publisher-hosted external data pages when they point to this paper's patient/variant data.
- Include DOI links only when they point to a dataset, not when they are ordinary article citations.
- Do not include ordinary PubMed/PMC/journal article pages, author emails, ORCID links, generic database homepages, recruitment/project homepages, code/software-only links, ClinVar/LOVD variant pages, RNA-seq/expression-only datasets, animal/preclinical datasets, or ordinary publisher supplementary-material anchor pages.
- Do not list supplementary files already attached locally unless the paper points to a separate external repository URL/DOI outside the journal package.
- If no external patient/variant dataset links are found, return "external_dataset_links": [].

## Gene-specific phenotype review items for {target_gene}:

Onset unit rule:
- {gene_guidance["onset_unit"]}

Follow-up age rule:
- {gene_guidance["followup_unit"]}

Phenotype category:
Choose the best matching category for `phenotype_review_category`.
{phenotype_categories}

Treatment/ASM response focus:
Put the requested medication/treatment response details in `targeted_treatment_response`.
{response_focus}

Gene-specific extraction hints:
{extraction_hints}

## Common phenotype fields to look for:
- Seizure types: focal, generalized, epileptic spasms, absence, myoclonic, etc.
- Age of onset: when seizures first appeared
- Response to ASM (Anti-Seizure Medications): controlled, refractory, partial response
- ASM used: VPA, LEV, LTG, OXC, TPM, etc.
- Epilepsy syndrome: DEE, West syndrome, Dravet syndrome, etc.
- Developmental outcomes: normal, delayed, intellectual disability
- Other features: autism, movement disorders, etc.

If NO patients are found, still return any external dataset links found:
{{"patients": [], "external_dataset_links": [...]}}

Return ONLY JSON, no additional text.
"""

    def create_reextract_prompt(self, target_gene, target_mentions):
        """기존 Not specified row만 재추출하기 위한 프롬프트"""
        target_payload = []
        for row in target_mentions:
            target_payload.append({
                "mention_id": row.get("mention_id"),
                "patient_id": row.get("patient_id"),
                "mention_source_type": row.get("mention_source_type"),
                "cited_source": row.get("cited_source"),
                "sex": row.get("sex"),
                "variant_type": row.get("variant_type"),
                "korean_patient_source": row.get("korean_patient_source"),
                "age_of_onset": row.get("age_of_onset"),
                "seizure_types": row.get("seizure_types"),
                "epilepsy_type": row.get("epilepsy_type"),
                "location_in_paper": row.get("location_in_paper"),
                "originality_evidence": row.get("originality_evidence"),
                "old_variant_c_gnomad": row.get("variant_c_gnomad"),
                "old_variant_p": row.get("variant_p"),
            })

        targets_json = json.dumps(target_payload, ensure_ascii=False, indent=2)
        gene_guidance = self.get_gene_review_guidance(target_gene)
        phenotype_categories = "\n".join(f"- {item}" for item in gene_guidance["phenotype_categories"])
        response_focus = "\n".join(f"- {item}" for item in gene_guidance["response_focus"])
        extraction_hints = "\n".join(f"- {item}" for item in gene_guidance["extraction_hints"])

        return f"""You are re-extracting previously incomplete patient-level genetics rows.
Analyze ALL attached files for the **{target_gene}** gene.

Your job is to replace ONLY the existing rows listed below. Return exactly one patient object
for each `mention_id`. Keep the same `mention_id` in the output. Do not create new rows.

Previously incomplete rows:
{targets_json}

For each row, find the matching patient/case/table row in the attached files and re-extract:
- variant_c_original
- variant_p
- nm_transcript
- variant_type
- Korean patient source classification
- patient-level phenotype fields
- source/originality evidence

If the row is an aggregate-only statement or the paper truly does not provide patient-specific
cDNA/protein variant information, keep `variant_c_original` and `variant_p` as "Not specified".
Do not invent variants from gene knowledge, ClinVar, or typical disease facts.

Return strict JSON only:
{{
  "patients": [
    {{
      "mention_id": 123,
      "patient_id": "same patient/case identifier when possible",
      "mention_source_type": "original_patient / cited_patient / aggregate_only / variant_only / previously_reported_patient / reanalyzed_previous_patient / unclear",
      "originality_evidence": "evidence for source type",
      "cited_source": "If cited, original source; otherwise 'Not specified'",
      "sex": "Male / Female / Not specified",
      "korean_patient_source": "Korean multicenter cohort / Korean patient from literature / Not Korean / Not specified",
      "korean_patient_evidence": "Evidence for Korean source classification, or 'Not specified'",
      "nm_transcript": "NM accession stated in paper or 'Not specified'",
      "variant_c_original": "c.123A>G or 'Not specified'",
      "variant_p": "p.Arg41His or 'Not specified'",
      "variant_type": "sequence_variant_or_small_indel / exon_level_cnv / large_cnv / unclear",
      "phenotype": {{
        "seizure_presence": "Yes / No / Not specified",
        "age_of_onset": "6 months / Not specified",
        "last_followup_age": "3 years / Not specified",
        "seizure_types": "Focal motor / Not specified",
        "phenotype_review_category": "One category below or 'Not specified'",
        "response_to_asm": "Controlled / Refractory / Unknown / Not specified",
        "asm_used": "VPA, LEV / Not specified",
        "responsive_asm": "VPA / Not specified",
        "targeted_treatment_response": "requested treatment response or 'Not specified'",
        "epilepsy_type": "DEE / Focal epilepsy / Not specified",
        "other_phenotypes": "Developmental delay, intellectual disability, etc."
      }},
      "location_in_paper": "precise table/case/supplement location",
      "confidence": "High / Medium / Low"
    }}
  ]
}}

Rules:
- Return exactly one object per requested `mention_id`.
- Preserve the `mention_id` exactly.
- Do not output `variant_c_gnomad`, `gnomad_transcript`, or `gnomad_conversion_status`; the script adds them.
- For `variant_c_original`, record only valid cDNA HGVS beginning with `c.` and without NM prefix when available. If the requested row is an explicit CNV/exon-level genotype and no cDNA HGVS exists, preserve the explicit genotype text and classify it in `variant_type`.
- For `variant_p`, normalize to one HGVS-style 3-letter protein value, e.g. p.E122K -> p.Glu122Lys.
- If a variant cannot be confidently linked to that specific patient/mention, use "Not specified".
- Prefer supplementary tables and patient-level tables over narrative summaries.
- Do not create patient-level phenotype from cohort-level summaries. If genotype and phenotype cannot be directly matched to the same patient/case/table row, use `mention_source_type="variant_only"` or `unclear`, keep phenotype fields as "Not specified", and explain why in `originality_evidence` or `location_in_paper`.
- Before assigning "original_patient", check Methods, Participants, table notes, and Supplementary files for previous-cohort statements such as previously reported, previously described, reported by, adapted from, or re-analyzed.
- For `variant_type`, classify sequence/small indel vs exon-level CNV vs large CNV conservatively.
- For `korean_patient_source`, distinguish project/Korean multicenter cohort patients from Korean patients found in publications; otherwise use "Not Korean / Not specified".

Gene-specific phenotype review items for {target_gene}:

Onset unit rule:
- {gene_guidance["onset_unit"]}

Follow-up age rule:
- {gene_guidance["followup_unit"]}

Phenotype category:
{phenotype_categories}

Treatment/ASM response focus:
{response_focus}

Gene-specific extraction hints:
{extraction_hints}

Return ONLY JSON, no additional text.
"""

    @staticmethod
    def _is_supplementary_pdf(filename):
        """PDF 파일명이 supplementary인지 판별"""
        if is_supplementary_file(filename):
            return True
        nl = filename.lower()
        # s001, s002 등 패턴 (-s001.pdf, _s001.pdf)
        if re.search(r'[-_]s\d{3}\.pdf$', nl):
            return True
        return False

    def get_files_by_type(self, article_dir):
        """논문 디렉토리에서 main article files와 supplementary files를 분리"""
        article_path = Path(article_dir)

        if not article_path.exists():
            logging.warning(f"디렉토리가 존재하지 않습니다: {article_dir}")
            return [], []

        main_pdfs = []
        supp_files = []

        for file in article_path.iterdir():
            if not file.is_file():
                continue
            invalid_reason = invalid_downloaded_html_reason(file)
            if invalid_reason:
                logging.warning(f"잘못 저장된 HTML 파일 제외: {file.name} ({invalid_reason})")
                continue
            ext = file.suffix.lower()
            if ext in HTML_EXTENSIONS:
                if is_main_fulltext_html_file(file.name):
                    main_pdfs.append(file)
                else:
                    logging.info(f"본문 후보가 아닌 HTML 파일 제외: {file.name}")
            elif ext in NXML_EXTENSIONS:
                main_pdfs.append(file)
            elif ext == '.pdf':
                if is_rendered_fulltext_pdf_file(file.name) and (article_path / "main.fulltext.txt").exists():
                    logging.info(f"TXT 본문이 있어 렌더링 PDF는 추출 입력에서 제외: {file.name}")
                elif self._is_supplementary_pdf(file.name):
                    supp_files.append(file)
                else:
                    main_pdfs.append(file)
            elif ext == '.txt' and is_main_fulltext_text_file(file.name):
                main_pdfs.append(file)
            elif ext in SUPPLEMENTARY_EXTENSIONS:
                supp_files.append(file)

        supp_files.extend(self._unpack_packed_supplements(article_path))
        return self._dedupe_inputs(sorted(main_pdfs), sorted(supp_files))

    @staticmethod
    def _dedupe_inputs(main_files, supp_files):
        """같은 내용이 여러 형식·사본으로 들어가지 않게 입력을 정리한다.

        - 바이트가 같은 파일은 하나만
        - 본문 XML(nxml/PMC xml)은 하나만 (.nxml 우선)
        - 본문 PDF가 여럿이면 PMC 사본(PMC*.pdf)은 제외
        - 본문 XML이 있으면 PMC 텍스트 사본(PMC*.txt)은 제외
        """
        seen, kept_main, kept_supp = set(), [], []

        def digest(path):
            import hashlib
            return hashlib.sha1(path.read_bytes()).hexdigest()

        xml_mains = [f for f in main_files if f.suffix.lower() in NXML_EXTENSIONS]
        keep_xml = sorted(xml_mains, key=lambda f: (f.suffix.lower() != ".nxml", f.name))[:1]
        pdf_mains = [f for f in main_files if f.suffix.lower() == ".pdf"]
        pmc_copy = re.compile(r"^PMC\d+(\.\d+)?\.", re.I)
        drop_pdf = {f for f in pdf_mains if pmc_copy.match(f.name)} if len(pdf_mains) > 1 else set()
        if drop_pdf == set(pdf_mains):
            drop_pdf = set(sorted(drop_pdf)[1:])
        for f in main_files:
            if f.suffix.lower() in NXML_EXTENSIONS and f not in keep_xml:
                logging.info(f"   중복 본문 XML 제외: {f.name}")
                continue
            if f in drop_pdf:
                logging.info(f"   PMC 사본 PDF 제외: {f.name}")
                continue
            h = digest(f)
            if h in seen:
                logging.info(f"   같은 내용 파일 제외: {f.name}")
                continue
            seen.add(h)
            kept_main.append(f)
        for f in supp_files:
            if keep_xml and f.suffix.lower() == ".txt" and pmc_copy.match(f.name):
                logging.info(f"   PMC 텍스트 사본 제외: {f.name}")
                continue
            h = digest(f)
            if h in seen:
                logging.info(f"   같은 내용 파일 제외: {f.name}")
                continue
            seen.add(h)
            kept_supp.append(f)
        return kept_main, kept_supp

    def _filter_spreadsheet_for_gene(self, file_path):
        """큰 다유전자 표에서 대상 유전자 행 + 머리 행만 텍스트로 남긴다. 걸러낼 필요 없으면 None."""
        suffix = file_path.suffix.lower()
        size = file_path.stat().st_size
        big_text = suffix == ".txt" and size >= LARGE_TEXT_FILTER_MIN_BYTES and not is_main_fulltext_text_file(file_path.name)
        if (suffix not in {".xlsx", ".xls", ".csv", ".tsv"} and not big_text) or size < SPREADSHEET_FILTER_MIN_BYTES:
            return None
        sheets = []
        try:
            if suffix == ".xlsx":
                import openpyxl
                wb = openpyxl.load_workbook(file_path, read_only=True, data_only=True)
                for ws in wb.worksheets:
                    sheets.append((ws.title, [["" if c is None else str(c) for c in row] for row in ws.iter_rows(values_only=True)]))
            elif suffix == ".xls":
                import xlrd
                wb = xlrd.open_workbook(str(file_path))
                for i in range(wb.nsheets):
                    sh = wb.sheet_by_index(i)
                    sheets.append((sh.name, [[str(v) for v in sh.row_values(r)] for r in range(sh.nrows)]))
            elif big_text:
                # 큰 텍스트 supplement(유전자 목록·원자료 등)는 줄 단위로 대상 유전자 줄만 남긴다
                sheets.append((file_path.name, [[line] for line in file_path.read_text(encoding="utf-8", errors="replace").splitlines()]))
            else:
                sep = "\t" if suffix == ".tsv" else ","
                rows = list(csv.reader(file_path.read_text(encoding="utf-8", errors="replace").splitlines(), delimiter=sep))
                sheets.append((file_path.name, rows))
        except Exception as e:
            logging.warning(f"표 읽기 실패로 원본 사용: {file_path.name} ({e})")
            return None
        total = sum(len(rows) for _, rows in sheets)
        gene_re = re.compile(rf"(?<![A-Za-z0-9]){re.escape(self.gene_name)}(?![A-Za-z0-9])", re.I)
        hit = sum(1 for _, rows in sheets for row in rows if any(gene_re.search(c) for c in row))
        if total < SPREADSHEET_FILTER_MIN_ROWS or hit > total * 0.5:
            return None
        lines = [f"Filtered copy of supplementary file {file_path.name}: header rows plus only the rows that mention "
                 f"{self.gene_name} ({hit} of {total} rows). Rows for other genes were removed before upload.", ""]
        for title, rows in sheets:
            gene_rows = [row for row in rows[SPREADSHEET_HEADER_ROWS:] if any(gene_re.search(c) for c in row)]
            if not gene_rows:
                continue
            lines.append(f"## Sheet: {title}")
            lines += ["\t".join(row) for row in rows[:SPREADSHEET_HEADER_ROWS]]
            lines += ["\t".join(row) for row in gene_rows]
            lines.append("")
        return "\n".join(lines)

    def _unpack_packed_supplements(self, article_path):
        """ZIP과 확장자 없는 supplementary를 임시 폴더에 업로드 가능한 파일로 풀어 반환한다.

        - 확장자 없는 PDF/DOC/XLS(OLE) → 확장자를 붙여 복사
        - 확장자 없는 DOCX/XLSX/PPTX(OOXML ZIP) → 확장자를 붙여 복사
        - 그 외 ZIP → 안의 PDF·문서·표 파일을 꺼냄 (폴더에 이미 풀려 있는 같은 이름, 심사 기록 제외)
        """
        self._cleanup_temp_dirs(self._packed_temp_dirs)
        self._packed_temp_dirs = []
        files = sorted(f for f in article_path.iterdir() if f.is_file())
        loose = {f.name.lower() for f in files}
        unpacked = []
        temp_dir = None

        def write(name, data):
            nonlocal temp_dir
            if temp_dir is None:
                temp_dir = Path(tempfile.mkdtemp(prefix="supp_unpack_"))
                self._packed_temp_dirs.append(temp_dir)
            path = temp_dir / name
            n = 1
            while path.exists():
                path = temp_dir / f"{Path(name).stem}_{n}{Path(name).suffix}"
                n += 1
            path.write_bytes(data)
            unpacked.append(path)
            logging.info(f"   [UNPACK] {name}")

        for file in files:
            if file.suffix.lower() not in ("", ".zip", ".bin"):
                continue
            kind = file_signature(file)
            try:
                if kind == "pdf":
                    write(f"{file.name}.pdf", file.read_bytes())
                elif kind == "ole":
                    ext = ole_document_extension(file.read_bytes())
                    if ext:
                        write(f"{file.name}{ext}", file.read_bytes())
                elif kind == "zip":
                    with zipfile.ZipFile(file) as archive:
                        ext = ooxml_extension(archive.namelist())
                        if ext:
                            write(f"{file.name}{ext}", file.read_bytes())
                            continue
                        members = [
                            info for info in archive.infolist()
                            if not info.is_dir()
                            and "__MACOSX" not in info.filename
                            and not Path(info.filename).name.startswith(".")
                            and Path(info.filename).suffix.lower() in (PDF_EXTENSIONS | SUPPLEMENTARY_EXTENSIONS)
                            and Path(info.filename).name.lower() not in loose
                        ]
                        skipped = [m for m in members if PEER_REVIEW_FILE_PATTERN.search(Path(m.filename).name)]
                        members = [m for m in members if m not in skipped]
                        if skipped:
                            logging.info(f"   ZIP 심사 기록 {len(skipped)}개 제외: {file.name}")
                        if len(members) > MAX_ZIP_MEMBERS:
                            logging.warning(f"ZIP 문서가 많아 앞 {MAX_ZIP_MEMBERS}개만 사용: {file.name} ({len(members)}개)")
                            members = members[:MAX_ZIP_MEMBERS]
                        for member in members:
                            write(f"{file.stem}__{Path(member.filename).name}", archive.read(member))
            except Exception as e:
                logging.warning(f"supplementary 풀기 실패로 제외: {file.name} ({e})")
        return unpacked

    @staticmethod
    def _normalize_extracted_text(text):
        text = html.unescape(text or "")
        text = re.sub(r"\s+", " ", text)
        return text.strip()

    @staticmethod
    def _find_first_by_local_name(root, tag_name):
        for elem in root.iter():
            local = elem.tag.rsplit("}", 1)[-1] if isinstance(elem.tag, str) else elem.tag
            if local == tag_name:
                return elem
        return None

    @classmethod
    def _text_from_element(cls, elem):
        if elem is None:
            return ""
        return cls._normalize_extracted_text(" ".join(elem.itertext()))

    def _nxml_to_plain_text(self, file_path):
        """JATS/NXML 본문을 모델 입력용 plain text로 변환한다."""
        try:
            root = ET.parse(file_path).getroot()
            sections = []

            title = self._text_from_element(self._find_first_by_local_name(root, "article-title"))
            if title:
                sections.append(f"TITLE: {title}")

            abstract = self._text_from_element(self._find_first_by_local_name(root, "abstract"))
            if abstract:
                sections.append(f"ABSTRACT: {abstract}")

            body = self._text_from_element(self._find_first_by_local_name(root, "body"))
            if body:
                sections.append(f"BODY: {body}")

            floats = self._text_from_element(self._find_first_by_local_name(root, "floats-group"))
            if floats and floats not in body:
                sections.append(f"FLOATS_AND_TABLES: {floats}")

            text = "\n\n".join(sections) or self._normalize_extracted_text(" ".join(root.itertext()))
        except Exception as e:
            logging.warning(f"NXML 파싱 실패, 태그 제거 fallback 사용: {file_path.name} ({e})")
            raw = file_path.read_text(encoding="utf-8", errors="replace")
            text = self._normalize_extracted_text(re.sub(r"<[^>]+>", " ", raw))

        if len(text) > MAX_NXML_TEXT_CHARS:
            logging.warning(
                f"NXML 텍스트가 길어 앞부분만 사용: {file_path.name} "
                f"({len(text)} chars -> {MAX_NXML_TEXT_CHARS})"
            )
            text = text[:MAX_NXML_TEXT_CHARS]

        return (
            f"Converted full-text XML source file: {file_path.name}\n"
            "Use this text as the main article body when extracting patient-level variants.\n\n"
            f"{text}\n"
        )

    def _html_to_plain_text(self, file_path):
        """PMC 본문 HTML을 모델 입력용 plain text로 변환한다."""
        raw = file_path.read_text(encoding="utf-8", errors="replace")
        parser = HTMLTextExtractor()
        parser.feed(raw)
        text = self._normalize_extracted_text(parser.text())
        if len(text) > MAX_NXML_TEXT_CHARS:
            logging.warning(
                f"HTML 텍스트가 길어 앞부분만 사용: {file_path.name} "
                f"({len(text)} chars -> {MAX_NXML_TEXT_CHARS})"
            )
            text = text[:MAX_NXML_TEXT_CHARS]

        return (
            f"Converted full-text HTML source file: {file_path.name}\n"
            "Use this text as the main article body when extracting patient-level variants.\n\n"
            f"{text}\n"
        )

    def _prepare_upload_files(self, file_paths):
        """XML/HTML은 임시 .txt로 변환하고, 나머지는 원본 경로를 그대로 반환한다."""
        prepared = []
        temp_dirs = []
        for file_path in file_paths:
            suffix = file_path.suffix.lower()
            filtered = self._filter_spreadsheet_for_gene(file_path)
            if filtered is not None:
                temp_dir = Path(tempfile.mkdtemp(prefix="gene_filtered_"))
                temp_path = temp_dir / f"{file_path.stem}.{self.gene_name}_rows.txt"
                temp_path.write_text(filtered, encoding="utf-8")
                prepared.append(temp_path)
                temp_dirs.append(temp_dir)
                logging.info(f"   [표 필터] {file_path.name} -> {temp_path.name}")
                continue
            if suffix not in (NXML_EXTENSIONS | HTML_EXTENSIONS):
                prepared.append(file_path)
                continue

            try:
                temp_dir = Path(tempfile.mkdtemp(prefix="nxml_extract_"))
                temp_path = temp_dir / f"{file_path.name}.txt"
                if suffix in NXML_EXTENSIONS:
                    text = self._nxml_to_plain_text(file_path)
                    label = "NXML->TXT"
                else:
                    text = self._html_to_plain_text(file_path)
                    label = "HTML->TXT"
                temp_path.write_text(text, encoding="utf-8")
                if temp_path.stat().st_size > MAX_UPLOAD_FILE_BYTES:
                    mb = temp_path.stat().st_size / (1024 * 1024)
                    logging.warning(f"텍스트 변환 크기 초과로 제외: {file_path.name} ({mb:.1f} MB)")
                    shutil.rmtree(temp_dir, ignore_errors=True)
                    continue
                prepared.append(temp_path)
                temp_dirs.append(temp_dir)
                logging.info(f"   [{label}] {file_path.name} -> {temp_path.name}")
            except Exception as e:
                logging.warning(f"텍스트 변환 실패로 제외: {file_path.name} ({e})")

        return prepared, temp_dirs

    @staticmethod
    def _cleanup_temp_dirs(temp_dirs):
        for temp_dir in temp_dirs:
            shutil.rmtree(temp_dir, ignore_errors=True)

    def _upload_files(self, file_paths):
        """파일 업로드 후 [{file_id, filename}] 반환"""
        uploaded = []
        for file_path in file_paths:
            upload_path = file_path
            temp_path = None
            try:
                if file_path.suffix != file_path.suffix.lower():
                    temp = tempfile.NamedTemporaryFile(
                        suffix=file_path.suffix.lower(),
                        delete=False
                    )
                    temp.close()
                    temp_path = Path(temp.name)
                    shutil.copyfile(file_path, temp_path)
                    upload_path = temp_path

                logging.info(f"📤 Uploading {file_path.name}...")
                with open(upload_path, "rb") as f:
                    file_obj = self.client.files.create(
                        file=f,
                        purpose="user_data"
                    )
                uploaded.append({
                    'file_id': file_obj.id,
                    'filename': file_path.name
                })
                logging.info(f"   ✅ {file_obj.id}")
            except Exception as e:
                if "invalid_api_key" in str(e) or "Incorrect API key" in str(e):
                    raise InvalidApiKeyError("OpenAI API key가 유효하지 않습니다.") from e
                # 일부 파일이 실패해도 업로드된 파일이 있으면 분석은 계속 진행한다.
                logging.error(f"   ❌ 업로드 실패: {e}")
            finally:
                if temp_path and temp_path.exists():
                    temp_path.unlink()
        return uploaded

    def _filter_upload_candidates(self, file_paths):
        """OpenAI 파일 제한을 넘는 자료는 분석 후보에서 제외한다."""
        kept = []
        for file_path in file_paths:
            try:
                file_size = file_path.stat().st_size
            except OSError:
                logging.warning(f"파일 크기를 확인할 수 없어 제외: {file_path.name}")
                continue
            if file_size > MAX_UPLOAD_FILE_BYTES:
                mb = file_size / (1024 * 1024)
                logging.warning(f"파일 크기 초과로 제외: {file_path.name} ({mb:.1f} MB)")
                continue
            invalid_reason = invalid_downloaded_html_reason(file_path)
            if invalid_reason:
                logging.warning(f"잘못 저장된 HTML 파일 제외: {file_path.name} ({invalid_reason})")
                continue
            kept.append(file_path)
        return kept

    def _cleanup_files(self, uploaded_files):
        """업로드한 파일 삭제"""
        for file_info in uploaded_files:
            try:
                self.client.files.delete(file_info['file_id'])
            except Exception:
                pass

    def _call_api(self, uploaded_files, prompt_text):
        """Responses API 호출 후 response_text 반환"""
        content = []
        for file_info in uploaded_files:
            content.append({
                "type": "input_file",
                "file_id": file_info['file_id']
            })
        content.append({
            "type": "input_text",
            "text": prompt_text
        })

        response = self.client.responses.create(
            model=EXTRACTION_MODEL,
            input=[{
                "role": "user",
                "content": content
            }],
            timeout=RESPONSE_API_TIMEOUT_SECONDS,
        )
        if getattr(response, "status", None) == "incomplete":
            details = getattr(response, "incomplete_details", None)
            reason = getattr(details, "reason", None) or str(details)
            raise OutputTruncatedError(f"Responses API output incomplete ({reason}); use --chunk-size")
        if getattr(response, "output_text", None):
            return response.output_text

        for item in getattr(response, "output", []) or []:
            for part in getattr(item, "content", []) or []:
                text = getattr(part, "text", None)
                if text:
                    return text

        raise RuntimeError("Responses API returned no text output")

    def _extract_json(self, uploaded_files, prompt):
        """한 번에, 또는 chunk_size가 있으면 환자 행 범위를 나눠 여러 번 호출해 결과를 합친다."""
        if not self.chunk_size:
            return self._parse_json(self._call_api(uploaded_files, prompt))
        merged, seen, start, total = None, set(), 1, None
        for _ in range(MAX_OUTPUT_CHUNKS):
            end = start + self.chunk_size - 1
            chunk_prompt = prompt + f"""

## Chunked output (this paper has many patients)
Enumerate every individual {self.gene_name} patient-variant row in a fixed order: main-text cases first, then main tables in order, then supplementary tables/sheets in order (row order within each table).
Output ONLY rows {start} to {end} of that enumeration (1-based) in `patients`, and add "total_patient_rows": <count of all rows in the enumeration>.
If {start} is greater than the total, return "patients": [].
"""
            logging.info(f"   chunk rows {start}-{end}")
            part = self._parse_json(self._call_api(uploaded_files, chunk_prompt))
            patients = part.get("patients", []) or []
            if merged is None:
                merged = {**part, "patients": []}
            try:
                total = int(part.get("total_patient_rows")) if part.get("total_patient_rows") is not None else total
            except (TypeError, ValueError):
                pass
            for patient in patients:
                key = (str(patient.get("patient_id", "")).strip().lower(),
                       str(patient.get("variant_c_original", "")).strip().lower(),
                       str(patient.get("variant_p", "")).strip().lower())
                if key not in seen:
                    seen.add(key)
                    merged["patients"].append(patient)
            if not patients or (total is not None and end >= total):
                break
            start = end + 1
        merged["total_patient_rows"] = total
        logging.info(f"   chunk 합계 {len(merged['patients'])}행 (모델 집계 {total})")
        return merged

    @staticmethod
    def _should_retry_main_only(error):
        """supplementary까지 넣어 입력이 너무 커진 경우 main PDF만 재시도한다."""
        message = str(error)
        return any(snippet in message for snippet in CONTEXT_RETRY_ERROR_SNIPPETS)

    def _parse_json(self, response_text):
        """응답에서 JSON 파싱"""
        # 모델이 코드블록을 붙여도 JSON 객체 부분만 파싱한다.
        clean = response_text.replace("```json", "").replace("```", "").strip()
        start = clean.find('{')
        end = clean.rfind('}')
        if start != -1 and end != -1:
            return json.loads(clean[start:end+1])
        return json.loads(clean)

    def _normalize_patients_gnomad(self, patients):
        """
        환자 목록의 variant_c_original을 gnomAD 기준 transcript로 정규화.
        GPT가 추출한 variant_c_original은 논문 원문 c. 표기로 유지한다.
        정규화 결과는 variant_c_gnomad / gnomad_transcript / gnomad_conversion_status에 저장.
        """
        normalized = []
        for patient in patients:
            p = dict(patient)
            p['variant_p'] = normalize_variant_p(p.get('variant_p', ''))
            variant_c = p.get('variant_c_original', '')
            nm_source = p.get('nm_transcript', '') or ''
            if nm_source in ('Not specified', ''):
                nm_source = None

            gnomad_c, status = normalize_to_gnomad(variant_c, self.gene_name, nm_source)

            p['variant_c_gnomad'] = gnomad_c
            p['gnomad_transcript'] = GNOMAD_TRANSCRIPTS.get(self.gene_name.upper(), '')
            p['gnomad_conversion_status'] = status

            if status not in ('already_gnomad', 'assumed_gnomad', 'not_specified'):
                logging.info(f"   gnomAD 정규화 [{status}]: {variant_c} → {gnomad_c}")

            normalized.append(p)
            time.sleep(0.2)  # Mutalyzer rate limit

        return normalized

    def _correct_api_error_variant(self, patient):
        """
        Mutalyzer api_error가 난 cDNA 표기만 GPT로 재검토한다.
        확실한 형식 오류는 보정하고, 복구 불가능한 표기는 invalid_hgvs로 분류한다.
        """
        gnomad_transcript = GNOMAD_TRANSCRIPTS.get(self.gene_name.upper(), '')
        prompt = f"""
You are correcting a single HGVS cDNA variant for a patient-level genetics database.

Gene: {self.gene_name}
Target transcript for comparison: {gnomad_transcript}

Use only the information below. Do not invent a variant. If the cDNA notation cannot be
confidently corrected into valid HGVS cDNA syntax, return "Not specified".

Rules:
- Return the cDNA part only, without NM prefix.
- Correct obvious extraction/OCR mistakes only.
- Valid examples: c.2059C>T, c.333delA, c.435_447dup, c.1676+5G>A, c.190-2A>G.
- Invalid examples that should become Not specified unless clearly fixable:
  c.1922S>X, c.3244Q>X, c.322K>X, c.638K>X, c.1167_1168G>QfsX26.
- For frameshift protein notation mixed into cDNA, keep only the valid cDNA event if clear.
  Example: c.1329_1333delCAAGG with p.Lys444Glyfs*27 stays c.1329_1333delCAAGG.

Input:
- variant_c_original: {patient.get("variant_c_original", "")}
- variant_p: {patient.get("variant_p", "")}
- nm_transcript stated in paper: {patient.get("nm_transcript", "")}
- patient_id: {patient.get("patient_id", "")}
- location_in_paper: {patient.get("location_in_paper", "")}
- originality_evidence: {patient.get("originality_evidence", "")}

Return strict JSON only:
{{
  "corrected_variant_c_original": "c.... or Not specified",
  "correction_reason": "short reason",
  "confidence": "High or Medium or Low"
}}
"""
        try:
            response = self.client.responses.create(
                model=EXTRACTION_MODEL,
                input=prompt,
                text={"format": {"type": "json_object"}},
            )
            correction = json.loads(response.output_text)
        except Exception as e:
            logging.warning(f"variant QC 실패: {e}")
            patient['gnomad_conversion_status'] = 'mutalyzer_rejected'
            return patient

        corrected = (correction.get("corrected_variant_c_original") or "Not specified").strip()
        reason = correction.get("correction_reason", "")
        confidence = correction.get("confidence", "")
        patient['variant_qc'] = {
            "corrected_by": "batch_extract_patient_mentions.py",
            "correction_reason": reason,
            "correction_confidence": confidence,
        }

        if not looks_like_cdna(corrected):
            patient['variant_c_original'] = "Not specified"
            patient['variant_c_gnomad'] = "Not specified"
            patient['gnomad_transcript'] = gnomad_transcript
            patient['gnomad_conversion_status'] = 'invalid_hgvs'
            return patient

        nm_source = patient.get('nm_transcript', '') or ''
        if nm_source in ('Not specified', ''):
            nm_source = None
        gnomad_c, status = normalize_to_gnomad(corrected, self.gene_name, nm_source)
        if status == 'api_error':
            status = 'mutalyzer_rejected'

        patient['variant_c_original'] = corrected
        patient['variant_c_gnomad'] = gnomad_c
        patient['gnomad_transcript'] = gnomad_transcript
        patient['gnomad_conversion_status'] = status
        return patient

    def _qc_api_error_variants(self, patients):
        """정규화 후 api_error로 남은 patient만 저장 전에 보정한다."""
        corrected = []
        for patient in patients:
            if patient.get('gnomad_conversion_status') == 'api_error':
                before = patient.get('variant_c_original', '')
                patient = self._correct_api_error_variant(dict(patient))
                after = patient.get('variant_c_original', '')
                status = patient.get('gnomad_conversion_status', '')
                logging.info(f"   variant QC [{status}]: {before} → {after}")
                time.sleep(0.2)
            corrected.append(patient)
        return corrected

    def _annotate_patients_clinvar(self, patients):
        """variant_c_gnomad 기준으로 ClinVar ID와 pathogenicity를 붙인다."""
        if not patients:
            return patients
        if not all([query_clinvar, save_cache, usable_variant, normalize_cache_entry, self.clinvar_session]):
            logging.warning("ClinVar 조회 모듈을 사용할 수 없어 ClinVar 컬럼을 비워 저장합니다.")
            return [
                {
                    **dict(patient),
                    "clinvar": dict(patient).get("clinvar") or "not_queried",
                    "clinvar_pathogenicity": dict(patient).get("clinvar_pathogenicity", ""),
                    "clinvar_id": dict(patient).get("clinvar_id", ""),
                }
                for patient in patients
            ]

        variants = sorted({
            usable_variant(patient.get("variant_c_gnomad", ""))
            for patient in patients
        })
        variants = [variant for variant in variants if variant]

        for variant in variants:
            cached = normalize_cache_entry(self.clinvar_cache.get(variant, ""))
            if needs_detailed_refresh and needs_detailed_refresh(cached):
                try:
                    if fetch_summaries_by_uids and summarize_summary_items:
                        summaries = fetch_summaries_by_uids(
                            self.clinvar_session,
                            cached.get("uids", []),
                            email=NCBI_EMAIL,
                        )
                        items = [
                            summaries[uid]
                            for uid in cached.get("uids", [])
                            if uid in summaries
                        ]
                        if items:
                            updated = summarize_summary_items(items, fallback=cached)
                            updated["uids"] = cached.get("uids", [])
                            self.clinvar_cache[variant] = updated
                        else:
                            self.clinvar_cache[variant] = cached
                    else:
                        self.clinvar_cache[variant] = query_clinvar(
                            variant,
                            self.clinvar_session,
                            email=NCBI_EMAIL,
                        )
                except Exception as exc:
                    logging.warning(f"ClinVar legacy 상세 조회 실패 ({variant}): {exc}")
                    self.clinvar_cache[variant] = cached
                save_cache(CLINVAR_CACHE_PATH, self.clinvar_cache)
                time.sleep(0.34)
                cached = normalize_cache_entry(self.clinvar_cache.get(variant, ""))

            if cached.get("status") not in {"query_error", "unknown", "not_queried"}:
                continue
            try:
                self.clinvar_cache[variant] = query_clinvar(
                    variant,
                    self.clinvar_session,
                    email=NCBI_EMAIL,
                )
            except Exception as exc:
                self.clinvar_cache[variant] = {
                    "status": "query_error",
                    "clinvar": f"query_error: {type(exc).__name__}",
                    "clinvar_pathogenicity": "",
                    "clinvar_id": "",
                    "classification_raw": "",
                    "uids": [],
                }
                logging.warning(f"ClinVar 조회 실패 ({variant}): {exc}")
            save_cache(CLINVAR_CACHE_PATH, self.clinvar_cache)
            time.sleep(0.34)

        annotated = []
        for patient in patients:
            p = dict(patient)
            variant = usable_variant(p.get("variant_c_gnomad", ""))
            entry = normalize_cache_entry(self.clinvar_cache.get(variant, "")) if variant else normalize_cache_entry("not_queried")
            p["clinvar"] = entry.get("clinvar", "not_queried")
            p["clinvar_pathogenicity"] = entry.get("clinvar_pathogenicity", "")
            p["clinvar_id"] = entry.get("clinvar_id", "")
            # P/LP가 아닌 분류(VUS·상충·benign)도 구분할 수 있게 ClinVar 원 분류를 보존한다.
            p["clinvar_classification"] = entry.get("classification_raw", "")
            annotated.append(p)
        return annotated

    def extract_from_article(self, article_dir, pmcid="", pmid=""):
        """단일 논문 폴더에서 환자별 정보 추출 (main PDF + supplementary 전체 분석)"""
        self.stats['total_articles'] += 1
        article_label = pmcid or f"PMID:{pmid}"

        logging.info(f"\n{'='*60}")
        logging.info(f"Processing {article_label}")
        logging.info(f"Directory: {article_dir}")
        logging.info(f"{'='*60}")

        # 1. PDF와 supplementary 파일 분리
        pdf_files, supp_files = self.get_files_by_type(article_dir)

        if not pdf_files and not supp_files:
            logging.warning(f"{article_label}: 분석 가능한 파일 없음")
            self.stats['failed_articles'] += 1
            return {
                'pmcid': pmcid,
                'pmid': pmid,
                'gene': self.gene_name,
                'status': 'no_files',
                'files_processed': [],
                'used_supplementary': False,
                'patients': [],
                'error': 'No supported files found'
            }

        logging.info(f"📁 Main files: {len(pdf_files)}개, Supplementary: {len(supp_files)}개")
        for f in pdf_files:
            if f.suffix.lower() in NXML_EXTENSIONS:
                label = "NXML"
            elif f.suffix.lower() in HTML_EXTENSIONS:
                label = "HTML"
            else:
                label = "PDF"
            logging.info(f"   [{label}] {f.name}")
        for f in supp_files:
            logging.info(f"   [SUP] {f.name}")

        all_files = self._filter_upload_candidates(pdf_files + supp_files)
        temp_dirs = []
        all_files, temp_dirs = self._prepare_upload_files(all_files)
        uploaded_files = self._upload_files(all_files)
        # DB에는 후보 파일이 아니라 실제 API 분석에 들어간 파일만 기록한다.
        uploaded_filenames = [file_info['filename'] for file_info in uploaded_files]
        if not uploaded_files:
            self._cleanup_temp_dirs(temp_dirs)
            self.stats['failed_articles'] += 1
            return {
                'pmcid': pmcid,
                'pmid': pmid,
                'gene': self.gene_name,
                'status': 'upload_failed',
                'files_processed': [],
                'used_supplementary': len(supp_files) > 0,
                'patients': [],
                'error': 'File upload failed'
            }

        try:
            # 한 번의 Responses API 호출에 main/supplementary 전체를 같이 넣는다.
            prompt = self.create_extraction_prompt(self.gene_name)

            logging.info(f"🤖 파일 {len(uploaded_files)}개 분석 중...")
            try:
                result = self._extract_json(uploaded_files, prompt)
            except Exception as e:
                if not supp_files or not self.allow_main_only_retry or not self._should_retry_main_only(e):
                    raise
                logging.warning(f"입력이 너무 커서 main PDF만 재시도: {e}")
                self._cleanup_files(uploaded_files)
                self._cleanup_temp_dirs(temp_dirs)
                main_only_files = self._filter_upload_candidates(pdf_files)
                main_only_files, temp_dirs = self._prepare_upload_files(main_only_files)
                uploaded_files = self._upload_files(main_only_files)
                uploaded_filenames = [file_info['filename'] for file_info in uploaded_files]
                if not uploaded_files:
                    raise
                logging.info(f"🤖 main file {len(uploaded_files)}개로 재분석 중...")
                result = self._extract_json(uploaded_files, prompt)

            patients = result.get('patients', [])
            supplement_note = clean_text(result.get('missing_supplementary_note', ''))
            if supplement_note.lower() in {"not specified", "none", "n/a"}:
                supplement_note = ""
            if supplement_note:
                logging.warning(f"   본문이 가리키는 supplementary 미확보: {supplement_note}")
            external_dataset_links = result.get('external_dataset_links', [])
            if not isinstance(external_dataset_links, list):
                external_dataset_links = []
            patients = self._normalize_patients_gnomad(patients)
            patients = self._qc_api_error_variants(patients)
            patients = self._annotate_patients_clinvar(patients)

            logging.info(f"   환자 {len(patients)}명 발견")
            logging.info(f"   외부 dataset 링크 후보 {len(external_dataset_links)}개 발견")

        except (InvalidApiKeyError, FatalApiQuotaError):
            raise
        except Exception as e:
            if is_quota_error(e):
                raise FatalApiQuotaError("OpenAI API quota exceeded; stopping batch.") from e
            logging.error(f"❌ 분석 실패: {e}")
            self.stats['failed_articles'] += 1
            self.stats['api_errors'] += 1
            return {
                'pmcid': pmcid,
                'pmid': pmid,
                'gene': self.gene_name,
                'status': 'api_error',
                'files_processed': uploaded_filenames,
                'used_supplementary': len(supp_files) > 0,
                'patients': [],
                'error': str(e)
            }
        finally:
            self._cleanup_files(uploaded_files)
            self._cleanup_temp_dirs(temp_dirs)

        # ── 결과 구성 ──
        self.stats['total_patients'] += len(patients)
        self.stats['external_links_found'] += len(external_dataset_links)
        self.stats['successful_articles'] += 1

        files_processed = uploaded_filenames
        supp_names = {f.name for f in supp_files}
        used_supplementary = any(name in supp_names for name in files_processed)

        logging.info(f"✅ 최종: {len(patients)}명 (supplementary 사용: {used_supplementary})")

        return {
            'pmcid': pmcid,
            'pmid': pmid,
            'gene': self.gene_name,
            'status': 'success',
            'files_processed': files_processed,
            'used_supplementary': used_supplementary,
            'supp_files_available': [f.name for f in supp_files],
            'supplement_note': supplement_note,
            'patients': patients,
            'external_dataset_links': external_dataset_links,
            'timestamp': datetime.now().isoformat()
        }

    def reextract_mentions_from_article(self, article_dir, target_mentions, pmcid="", pmid=""):
        """기존 mention_id 목록만 재추출한다. 결과는 append가 아니라 UPDATE용이다."""
        self.stats['total_articles'] += 1
        article_label = pmcid or f"PMID:{pmid}"

        logging.info(f"\n{'='*60}")
        logging.info(f"Re-extracting {len(target_mentions)} incomplete mentions from {article_label}")
        logging.info(f"Directory: {article_dir}")
        logging.info(f"{'='*60}")

        pdf_files, supp_files = self.get_files_by_type(article_dir)

        if not pdf_files and not supp_files:
            logging.warning(f"{article_label}: 분석 가능한 파일 없음")
            self.stats['failed_articles'] += 1
            return {
                'pmcid': pmcid,
                'pmid': pmid,
                'gene': self.gene_name,
                'status': 'no_files',
                'files_processed': [],
                'used_supplementary': False,
                'patients': [],
                'error': 'No supported files found',
                'timestamp': datetime.now().isoformat(),
            }

        logging.info(f"📁 Main files: {len(pdf_files)}개, Supplementary: {len(supp_files)}개")
        all_files = self._filter_upload_candidates(pdf_files + supp_files)
        temp_dirs = []
        all_files, temp_dirs = self._prepare_upload_files(all_files)
        uploaded_files = self._upload_files(all_files)
        uploaded_filenames = [file_info['filename'] for file_info in uploaded_files]
        if not uploaded_files:
            self._cleanup_temp_dirs(temp_dirs)
            self.stats['failed_articles'] += 1
            return {
                'pmcid': pmcid,
                'pmid': pmid,
                'gene': self.gene_name,
                'status': 'upload_failed',
                'files_processed': [],
                'used_supplementary': len(supp_files) > 0,
                'patients': [],
                'error': 'File upload failed',
                'timestamp': datetime.now().isoformat(),
            }

        try:
            prompt = self.create_reextract_prompt(self.gene_name, target_mentions)

            logging.info(f"🤖 파일 {len(uploaded_files)}개로 기존 mention 재추출 중...")
            try:
                response_text = self._call_api(uploaded_files, prompt)
            except Exception as e:
                if not supp_files or not self._should_retry_main_only(e):
                    raise
                logging.warning(f"입력이 너무 커서 main PDF만 재시도: {e}")
                self._cleanup_files(uploaded_files)
                self._cleanup_temp_dirs(temp_dirs)
                main_only_files = self._filter_upload_candidates(pdf_files)
                main_only_files, temp_dirs = self._prepare_upload_files(main_only_files)
                uploaded_files = self._upload_files(main_only_files)
                uploaded_filenames = [file_info['filename'] for file_info in uploaded_files]
                if not uploaded_files:
                    raise
                logging.info(f"🤖 main file {len(uploaded_files)}개로 재추출 중...")
                response_text = self._call_api(uploaded_files, prompt)

            result = self._parse_json(response_text)
            patients = result.get('patients', [])
            patients = self._normalize_patients_gnomad(patients)
            patients = self._qc_api_error_variants(patients)
            patients = self._annotate_patients_clinvar(patients)

            requested_ids = {int(row["mention_id"]) for row in target_mentions}
            filtered = []
            for patient in patients:
                try:
                    mention_id = int(patient.get("mention_id"))
                except (TypeError, ValueError):
                    continue
                if mention_id in requested_ids:
                    filtered.append(patient)
            patients = filtered

            logging.info(f"   재추출 mention {len(patients)}/{len(target_mentions)}개 반환")

        except (InvalidApiKeyError, FatalApiQuotaError):
            raise
        except Exception as e:
            if is_quota_error(e):
                raise FatalApiQuotaError("OpenAI API quota exceeded; stopping batch.") from e
            logging.error(f"❌ 재추출 실패: {e}")
            self.stats['failed_articles'] += 1
            self.stats['api_errors'] += 1
            return {
                'pmcid': pmcid,
                'pmid': pmid,
                'gene': self.gene_name,
                'status': 'api_error',
                'files_processed': uploaded_filenames,
                'used_supplementary': len(supp_files) > 0,
                'patients': [],
                'error': str(e),
                'timestamp': datetime.now().isoformat(),
            }
        finally:
            self._cleanup_files(uploaded_files)
            self._cleanup_temp_dirs(temp_dirs)

        self.stats['total_patients'] += len(patients)
        self.stats['successful_articles'] += 1

        files_processed = uploaded_filenames
        supp_names = {f.name for f in supp_files}
        used_supplementary = any(name in supp_names for name in files_processed)

        return {
            'pmcid': pmcid,
            'pmid': pmid,
            'gene': self.gene_name,
            'status': 'success',
            'files_processed': files_processed,
            'used_supplementary': used_supplementary,
            'supp_files_available': [f.name for f in supp_files],
            'patients': patients,
            'timestamp': datetime.now().isoformat(),
        }

    def print_stats(self):
        """통계 출력"""
        logging.info(f"\n{'='*60}")
        logging.info("FINAL STATISTICS")
        logging.info(f"{'='*60}")
        logging.info(f"Total articles processed: {self.stats['total_articles']}")
        logging.info(f"Successful: {self.stats['successful_articles']}")
        logging.info(f"Failed: {self.stats['failed_articles']}")
        logging.info(f"Total patients found: {self.stats['total_patients']}")
        logging.info(f"External dataset link candidates found: {self.stats['external_links_found']}")
        logging.info(f"API errors: {self.stats['api_errors']}")
        logging.info(f"")


def load_article_csv(path):
    import pandas as pd
    df = pd.read_csv(path, dtype=str)
    if "downloaded" in df.columns:
        df = df[df["downloaded"].fillna("Yes") == "Yes"]
    rows = []
    for _, row in df.iterrows():
        pmcid = "" if pd.isna(row.get("pmcid")) else str(row.get("pmcid")).strip()
        pmid = "" if pd.isna(row.get("pmid")) else str(row.get("pmid")).strip()
        if pmcid or pmid:
            rows.append({"pmcid": pmcid, "pmid": pmid})
    return rows


def run_staged(extractor, args, base_dir, stage_dir):
    """DB 없이 추출만 하고 결과 JSON을 저장한다. 이미 저장된 JSON은 건너뛴다 (재실행 안전)."""
    gene = args.gene.upper()
    articles = [{"pmcid": "", "pmid": args.pmid}] if args.pmid else load_article_csv(args.article_csv)
    if args.max_articles:
        articles = articles[:args.max_articles]
    logging.info(f"[stage] {gene}: {len(articles)}개 논문 -> {stage_dir}")
    for i, article in enumerate(articles, 1):
        out = stage_dir / f"{gene}_{article['pmid'] or article['pmcid']}.json"
        if out.exists() and not args.force:
            logging.info(f"[{i}/{len(articles)}] {out.name} 이미 있음 (스킵)")
            continue
        logging.info(f"[{i}/{len(articles)}] PMID:{article['pmid']}")
        result = extractor.extract_from_article(base_dir / (article["pmid"] or article["pmcid"]),
                                                article["pmcid"], article["pmid"])
        result["stage_settings"] = {"chunk_size": extractor.chunk_size,
                                    "allow_main_only_retry": extractor.allow_main_only_retry}
        out.write_text(json.dumps(result, ensure_ascii=False, indent=1), encoding="utf-8")
        logging.info(f"   저장: {out.name} (status {result.get('status')}, patients {len(result.get('patients', []))})")
    extractor.print_stats()


def main():
    parser = argparse.ArgumentParser(description='환자별 변이 + Phenotype 추출')
    parser.add_argument('--pmcid', help='단일 PMCID 처리')
    parser.add_argument('--pmid', help='단일 PMID 처리 (PMCID가 없는 논문용)')
    parser.add_argument('--all', action='store_true', help='전체 논문 폴더 배치 처리')
    parser.add_argument('--gene', default=TARGET_GENE, help=f'유전자명 (기본: {TARGET_GENE})')
    parser.add_argument('--input_dir', default='article_files',
                        help='논문 폴더가 있는 디렉토리')
    parser.add_argument('--db-dir', default='databases',
                        help='유전자별 SQLite DB 저장 디렉토리')
    parser.add_argument('--max-articles', type=int, default=None, help='최대 처리 개수')
    parser.add_argument('--article-csv', default=None,
                        help='처리할 논문 목록 CSV (downloaded==Yes이고 pmcid 또는 pmid가 있는 행만 처리)')
    parser.add_argument('--force', action='store_true',
                        help='이미 성공 처리된 논문도 다시 처리')
    parser.add_argument('--replace-article', action='store_true',
                        help='일반 추출 성공 시 해당 논문의 기존 patient_mentions 행을 삭제한 뒤 새 결과로 교체')
    parser.add_argument('--replace-not-specified', action='store_true',
                        help='variant_c_gnomad와 variant_p가 모두 Not specified인 기존 mention row만 재추출해 UPDATE')
    parser.add_argument('--replace-before',
                        help='--replace-not-specified에서 이 ISO timestamp 이전 extraction_timestamp row만 대상')
    parser.add_argument('--stage-dir',
                        help='DB에 쓰지 않고 추출 결과를 {GENE}_{PMID}.json으로 저장 (이미 처리된 논문도 다시 추출). '
                             '기존 행과 대조해 빠진 환자만 추가하는 재추출용: scripts/curation/merge_staged_mentions.py')
    parser.add_argument('--chunk-size', type=int, default=None,
                        help='환자 행을 이 크기로 나눠 여러 번 호출 (큰 환자 표용, 예: 40)')
    parser.add_argument('--no-main-only-retry', action='store_true',
                        help='입력이 너무 클 때 supplement를 빼고 main만으로 재시도하지 않음 (--stage-dir이면 자동)')
    parser.add_argument('--model', default=None,
                        help=f'OpenAI 모델 (기본: OPENAI_MODEL 환경변수 또는 {EXTRACTION_MODEL})')

    args = parser.parse_args()
    if args.model:
        globals()["EXTRACTION_MODEL"] = args.model
    logging.info(f"추출 모델: {globals()['EXTRACTION_MODEL']}")

    workspace = WORKSPACE_ROOT
    # CLI 경로는 repo root 기준으로 해석한다.
    base_dir = workspace / args.input_dir
    db_dir = workspace / args.db_dir
    db_path = db_dir / f"{args.gene.upper()}.sqlite"

    # Extractor 초기화
    extractor = PhenotypeExtractor(
        api_key=OPENAI_API_KEY,
        gene_name=args.gene
    )
    extractor.chunk_size = args.chunk_size
    stage_dir = Path(args.stage_dir) if args.stage_dir else None
    if stage_dir or args.no_main_only_retry:
        extractor.allow_main_only_retry = False
    if stage_dir:
        if not stage_dir.is_absolute():
            stage_dir = workspace / stage_dir
        stage_dir.mkdir(parents=True, exist_ok=True)
        if args.pmcid or args.replace_not_specified or args.replace_article:
            parser.error("--stage-dir는 --pmid/--article-csv와만 사용 (DB를 건드리지 않음)")
        run_staged(extractor, args, base_dir, stage_dir)
        return
    patient_db = PatientGeneDatabase(db_path)

    logging.info("=" * 60)
    logging.info("Patient-level Phenotype Extraction")
    logging.info(f"Gene: {args.gene}")
    logging.info(f"Input: {base_dir}")
    logging.info(f"DB: {db_path}")
    logging.info("=" * 60)

    def article_dir_for(article):
        if article["pmid"]:
            return base_dir / article["pmid"]
        return base_dir / article["pmcid"]

    # 처리할 논문 리스트를 단일 ID, CSV, 디렉토리 전체 중 하나로 만든다.
    if args.pmcid:
        article_list = [{"pmcid": args.pmcid, "pmid": patient_db.pmid_for_pmcid(args.pmcid)}]
    elif args.pmid:
        article_list = [{"pmcid": "", "pmid": args.pmid}]
    elif args.article_csv:
        import pandas as pd
        df = pd.read_csv(args.article_csv)
        article_list = []
        for _, row in df[df['downloaded'] == 'Yes'].iterrows():
            pmcid = "" if "pmcid" not in row or pd.isna(row.get("pmcid")) else str(row.get("pmcid")).strip()
            pmid = "" if "pmid" not in row or pd.isna(row.get("pmid")) else str(row.get("pmid")).strip()
            if pmcid or pmid:
                article_list.append({"pmcid": pmcid, "pmid": pmid})
        logging.info(f"CSV에서 {len(article_list)}개 논문 ID 로드")
    elif args.all:
        article_list = []
        for d in sorted([d for d in base_dir.iterdir() if d.is_dir()]):
            if d.name.isdigit():
                article_list.append({"pmcid": "", "pmid": d.name})
    else:
        parser.error("--pmcid, --pmid, --all, 또는 --article-csv 중 하나 필요")
        return

    if args.max_articles:
        article_list = article_list[:args.max_articles]

    logging.info(f"총 {len(article_list)}개 논문 처리 예정")

    # 배치 처리
    try:
        for i, article in enumerate(article_list, 1):
            pmcid = article["pmcid"]
            pmid = article["pmid"]
            article_label = pmcid or f"PMID:{pmid}"
            logging.info(f"\n[{i}/{len(article_list)}] {article_label}")

            if args.replace_not_specified:
                target_mentions = patient_db.not_specified_mentions_for_article(
                    pmcid,
                    pmid,
                    before_timestamp=args.replace_before,
                )
                if not target_mentions:
                    logging.info(f"{article_label}: 교체 대상 Not specified mention 없음 (스킵)")
                    continue

                article_dir = article_dir_for(article)
                patient_db.refresh_article_files(article_dir, pmcid, pmid)
                result = extractor.reextract_mentions_from_article(
                    article_dir,
                    target_mentions,
                    pmcid,
                    pmid,
                )
                if result.get("status") == "success":
                    summary = patient_db.replace_mentions(result)
                    logging.info(
                        f"DB UPDATE 완료: updated={summary['updated']}, "
                        f"skipped={summary['skipped']}, returned={summary['returned']}"
                    )
                else:
                    logging.warning(f"{article_label}: 재추출 실패 상태라 DB UPDATE 건너뜀 ({result.get('status')})")
                continue

            if patient_db.has_successful_article(pmcid, pmid) and not args.force:
                logging.info(f"{article_label}: 이미 DB에 성공 처리됨 (스킵)")
                continue

            article_dir = article_dir_for(article)
            patient_db.refresh_article_files(article_dir, pmcid, pmid)
            result = extractor.extract_from_article(article_dir, pmcid, pmid)
            if args.replace_article and result.get("status") == "success":
                deleted = patient_db.delete_article_mentions(pmcid, pmid)
                logging.info(f"기존 patient_mentions 교체 삭제: {deleted}개")
            patient_db.append_result(result)
            logging.info(f"DB append 완료: {db_path}")

            # 중간 통계
            if i % 10 == 0:
                logging.info(f"\n📊 Progress: {i}/{len(article_list)}")
                logging.info(f"   성공: {extractor.stats['successful_articles']}")
                logging.info(f"   환자: {extractor.stats['total_patients']}명")
    finally:
        patient_db.close()

    extractor.print_stats()

    logging.info(f"\nDB 저장: {db_path}")
    logging.info(f"{'='*60}")


if __name__ == '__main__':
    main()
