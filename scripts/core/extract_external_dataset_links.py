#!/usr/bin/env python3
"""
논문 본문/파일에서 환자-변이정보가 들어 있을 수 있는 외부 dataset 링크만 추출한다.

기존 patient extraction 파이프라인의 파일 업로드/Responses API 흐름을 가볍게 재사용하되,
환자 row, ClinVar, Mutalyzer 처리는 하지 않는다. 결과는 유전자별 SQLite DB의
external_dataset_links 테이블에 저장한다.
"""

import argparse
import json
import logging
import os
import re
import sqlite3
import tempfile
import shutil
from datetime import datetime
from pathlib import Path

try:
    import pandas as pd
except ImportError:
    pd = None

try:
    from openai import OpenAI
except ImportError:
    OpenAI = None

try:
    from batch_extract_patient_mentions import (
        FatalApiQuotaError,
        InvalidApiKeyError,
        PatientGeneDatabase,
        RESPONSE_API_TIMEOUT_SECONDS,
        is_supplementary_file,
        is_quota_error,
        load_env_file,
    )
except ImportError:
    from scripts.core.batch_extract_patient_mentions import (
        FatalApiQuotaError,
        InvalidApiKeyError,
        PatientGeneDatabase,
        RESPONSE_API_TIMEOUT_SECONDS,
        is_supplementary_file,
        is_quota_error,
        load_env_file,
    )


logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
)

WORKSPACE = Path(__file__).resolve().parents[2]
load_env_file(WORKSPACE / ".env")
OPENAI_API_KEY = os.getenv("OPENAI_API_KEY")

DEFAULT_MODEL = os.getenv("EXTERNAL_LINK_MODEL", "gpt-5.4-mini")
MAX_UPLOAD_FILE_BYTES = int(os.getenv("EXTERNAL_LINK_MAX_UPLOAD_MB", "32")) * 1024 * 1024

PDF_EXTENSIONS = {".pdf"}


def is_supplementary_pdf(filename):
    nl = filename.lower()
    if is_supplementary_file(filename):
        return True
    supplementary_substrings = (
        "datasheet",
        "data_sheet",
        "sdata",
        "presentation_",
        "media-",
        "fig-data",
        "-fig",
        "_fig",
        "figure",
        "moesm",
        "mmc",
        "suppinfo",
        "supporting",
        "appendix",
        "review_history",
    )
    if any(token in nl for token in supplementary_substrings):
        return True
    if "-fig" in nl and "-data" in nl:
        return True
    if re.search(r"fig\d+-data\d*\.pdf$", nl):
        return True
    supplementary_patterns = (
        "-s001.pdf",
        "_s001.pdf",
        "-s002.pdf",
        "_s002.pdf",
        ".s001.pdf",
        ".s002.pdf",
        ".sapp.pdf",
        "_sapp.pdf",
        "-sapp.pdf",
    )
    if nl.endswith(supplementary_patterns):
        return True
    if re.search(r"[_\.-]app\d+\.pdf$", nl):
        return True
    return False


def clean_text(value):
    if value is None:
        return ""
    return str(value).strip()


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


def should_keep_link(link, target_gene=""):
    """외부 dataset landing/download 링크가 아닌 obvious false positive를 제외한다."""
    repository = clean_text(link.get("repository")).lower()
    source_type = clean_text(link.get("source_type")).lower()
    expected_content = clean_text(link.get("expected_content")).lower()
    evidence_location = clean_text(link.get("evidence_location")).lower()
    link_text = clean_text(link.get("link_text")).lower()
    url = clean_text(link.get("url")).lower()
    combined_text = " ".join([repository, source_type, expected_content, evidence_location, link_text, url])

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

    dataset_link_signals = (
        "dryad",
        "figshare",
        "zenodo",
        "mendeley",
        "osf.io",
        "dataverse",
        "dbgap",
        "ega-archive",
        "ebi.ac.uk/ega",
        "arrayexpress",
        "geo/query",
        "ncbi.nlm.nih.gov/geo",
        "sra",
        "bioproject",
        "dataset",
        "data set",
        "download",
        "accession",
        "repository",
    )
    if not any(signal in combined_text for signal in dataset_link_signals):
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
    has_dataset_signal = any(term in combined_text for term in dataset_terms)
    if has_dataset_signal:
        return True

    if repository in {"unknown", "", "publisher"} and expected_content in {"unclear", ""}:
        if "citation" in evidence_location or "footer" in evidence_location:
            return False
        if "doi.org/10." in url:
            return False

    return True


class ExternalDatasetLinkDB:
    """external_dataset_links 저장 전용 wrapper."""

    COLUMNS = [
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
        self.patient_db = PatientGeneDatabase(db_path)
        self.conn = self.patient_db.conn

    def close(self):
        self.patient_db.close()

    def pmid_for_pmcid(self, pmcid):
        return self.patient_db.pmid_for_pmcid(pmcid)

    def pmcid_for_pmid(self, pmid):
        row = self.conn.execute(
            "SELECT pmcid FROM pubmed_articles WHERE pmid = ?",
            (pmid,),
        ).fetchone()
        return row["pmcid"] if row and row["pmcid"] else ""

    def refresh_article_files(self, article_dir, pmcid="", pmid=""):
        self.patient_db.refresh_article_files(article_dir, pmcid, pmid)

    def has_links_for_article(self, pmcid="", pmid=""):
        if pmid:
            row = self.conn.execute(
                "SELECT 1 FROM external_dataset_links WHERE pmid = ? LIMIT 1",
                (pmid,),
            ).fetchone()
            return row is not None
        if pmcid:
            row = self.conn.execute(
                "SELECT 1 FROM external_dataset_links WHERE pmcid = ? LIMIT 1",
                (pmcid,),
            ).fetchone()
            return row is not None
        return False

    def upsert_links(self, result):
        now = datetime.now().isoformat(timespec="seconds")
        pmid = clean_text(result.get("pmid"))
        pmcid = clean_text(result.get("pmcid"))
        gene = clean_text(result.get("gene"))
        rows = []

        for link in result.get("external_dataset_links", []):
            if not should_keep_link(link, gene):
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

            row = {
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
                "discovered_at": now,
                "updated_at": now,
            }
            rows.append(row)

        with self.conn:
            for row in rows:
                self.conn.execute(
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
                        local_file_path = excluded.local_file_path,
                        notes = excluded.notes,
                        updated_at = excluded.updated_at
                    """,
                    row,
                )
        return len(rows)


class ExternalDatasetLinkExtractor:
    """Responses API로 외부 dataset 링크만 추출한다."""

    def __init__(self, api_key, gene_name, model=DEFAULT_MODEL):
        if OpenAI is None:
            raise ImportError("openai 패키지가 필요합니다. `pip install openai` 후 다시 실행하세요.")
        if not api_key:
            raise ValueError("OPENAI_API_KEY 환경변수가 필요합니다.")
        self.client = OpenAI(api_key=api_key, max_retries=0)
        self.gene_name = gene_name
        self.model = model
        self.stats = {
            "total_articles": 0,
            "successful_articles": 0,
            "failed_articles": 0,
            "links_found": 0,
            "api_errors": 0,
        }

    def get_scan_files(self, article_dir):
        article_path = Path(article_dir)
        if not article_path.exists():
            logging.warning(f"디렉토리가 존재하지 않습니다: {article_dir}")
            return []
        files = []
        for file_path in article_path.iterdir():
            if not file_path.is_file():
                continue
            if file_path.suffix.lower() in PDF_EXTENSIONS and not is_supplementary_pdf(file_path.name):
                files.append(file_path)
        return sorted(files)

    def _filter_upload_candidates(self, file_paths):
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
            kept.append(file_path)
        return kept

    def _upload_files(self, file_paths):
        uploaded = []
        for file_path in file_paths:
            upload_path = file_path
            temp_path = None
            try:
                if file_path.suffix != file_path.suffix.lower():
                    temp = tempfile.NamedTemporaryFile(suffix=file_path.suffix.lower(), delete=False)
                    temp.close()
                    temp_path = Path(temp.name)
                    shutil.copyfile(file_path, temp_path)
                    upload_path = temp_path

                logging.info(f"📤 Uploading {file_path.name}...")
                with open(upload_path, "rb") as f:
                    file_obj = self.client.files.create(file=f, purpose="user_data")
                uploaded.append({"file_id": file_obj.id, "filename": file_path.name})
            except Exception as exc:
                if "invalid_api_key" in str(exc) or "Incorrect API key" in str(exc):
                    raise InvalidApiKeyError("OpenAI API key가 유효하지 않습니다.") from exc
                logging.error(f"   ❌ 업로드 실패: {exc}")
            finally:
                if temp_path and temp_path.exists():
                    temp_path.unlink()
        return uploaded

    def _cleanup_files(self, uploaded_files):
        for file_info in uploaded_files:
            try:
                self.client.files.delete(file_info["file_id"])
            except Exception:
                pass

    def create_prompt(self, pmid="", pmcid=""):
        return f"""You are screening a genetics paper for external data links.

Gene of interest: {self.gene_name}
PMID: {pmid or "Not specified"}
PMCID: {pmcid or "Not specified"}

Task:
Find URLs, DOIs, accession links, or repository pages in the attached article files that may contain
patient-level genotype/variant/phenotype data for the paper. These are usually in sections such as
Data availability, Supplementary Data, Methods, Results, Key resources table, footnotes, or table notes.

Include links such as:
- Dryad, Figshare, Zenodo, Mendeley Data, OSF, GitHub, institutional repositories
- dbGaP/EGA/ClinVar/LOVD/DECIPHER or other accession pages when the paper says patient/variant data are there
- Publisher-hosted external data pages that are not ordinary PMC supplementary files
- DOI links that point to a dataset

Do NOT include:
- ordinary references to PubMed, PMC, journal article DOI pages, author email, ORCID, generic software homepages,
  or general database homepages unless the paper states that the paper's patient/variant data are available there
- recruitment websites, disease community websites, project homepages, institutional homepages, or patient sign-up pages
  unless the page is explicitly a dataset download/landing page with accession, DOI, or downloadable data file
- links that are only about code, software, reagents, protocols, or non-human/basic-science datasets
- ordinary publisher supplementary material pages such as a journal page ending in #supplementary-material
- ordinary supplementary filenames already attached locally, unless they are a separate repository URL/DOI outside the journal

Return strict JSON only:
{{
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

Rules:
- Be conservative. Only include links that plausibly contain patient, genotype, variant, phenotype, cohort, or raw sequencing data.
- If a dataset is outside the attached files and would need a separate browser/repository download, set manual_download_required to true.
- If the paper states all data are available from the corresponding author only, do not create a URL row unless there is an actual URL/DOI/accession.
- If no such links are found, return {{"external_dataset_links": []}}.
"""

    def _call_api(self, uploaded_files, prompt_text):
        content = [{"type": "input_file", "file_id": file_info["file_id"]} for file_info in uploaded_files]
        content.append({"type": "input_text", "text": prompt_text})
        response = self.client.responses.create(
            model=self.model,
            input=[{"role": "user", "content": content}],
            text={"format": {"type": "json_object"}},
            timeout=RESPONSE_API_TIMEOUT_SECONDS,
        )
        if getattr(response, "output_text", None):
            return response.output_text
        for item in getattr(response, "output", []) or []:
            for part in getattr(item, "content", []) or []:
                text = getattr(part, "text", None)
                if text:
                    return text
        raise RuntimeError("Responses API returned no text output")

    @staticmethod
    def _parse_json(response_text):
        clean = response_text.replace("```json", "").replace("```", "").strip()
        start = clean.find("{")
        end = clean.rfind("}")
        if start != -1 and end != -1:
            return json.loads(clean[start : end + 1])
        return json.loads(clean)

    def extract_from_article(self, article_dir, pmcid="", pmid=""):
        self.stats["total_articles"] += 1
        article_label = pmcid or f"PMID:{pmid}"
        logging.info(f"\n{'=' * 60}")
        logging.info(f"External dataset link scan: {article_label}")
        logging.info(f"Directory: {article_dir}")
        logging.info(f"{'=' * 60}")

        files = self._filter_upload_candidates(self.get_scan_files(article_dir))
        if not files:
            self.stats["failed_articles"] += 1
            return {
                "pmcid": pmcid,
                "pmid": pmid,
                "gene": self.gene_name,
                "status": "no_files",
                "files_processed": [],
                "external_dataset_links": [],
                "error": "No supported files found",
            }

        logging.info(f"메인 PDF 스캔 파일: {len(files)}개")
        for file_path in files:
            logging.info(f"   [FILE] {file_path.name}")

        uploaded_files = self._upload_files(files)
        uploaded_filenames = [file_info["filename"] for file_info in uploaded_files]
        if not uploaded_files:
            self.stats["failed_articles"] += 1
            return {
                "pmcid": pmcid,
                "pmid": pmid,
                "gene": self.gene_name,
                "status": "upload_failed",
                "files_processed": [],
                "external_dataset_links": [],
                "error": "File upload failed",
            }

        try:
            prompt = self.create_prompt(pmid=pmid, pmcid=pmcid)
            logging.info(f"🤖 {self.model}로 외부 dataset 링크 추출 중...")
            response_text = self._call_api(uploaded_files, prompt)
            parsed = self._parse_json(response_text)
            links = parsed.get("external_dataset_links", [])
            if not isinstance(links, list):
                links = []
            self.stats["successful_articles"] += 1
            self.stats["links_found"] += len(links)
            logging.info(f"   외부 링크 {len(links)}개 발견")
            return {
                "pmcid": pmcid,
                "pmid": pmid,
                "gene": self.gene_name,
                "status": "success",
                "files_processed": uploaded_filenames,
                "external_dataset_links": links,
                "timestamp": datetime.now().isoformat(timespec="seconds"),
            }
        except (InvalidApiKeyError, FatalApiQuotaError):
            raise
        except Exception as exc:
            if is_quota_error(exc):
                raise FatalApiQuotaError("OpenAI API quota exceeded; stopping batch.") from exc
            self.stats["failed_articles"] += 1
            self.stats["api_errors"] += 1
            logging.error(f"❌ 링크 추출 실패: {exc}")
            return {
                "pmcid": pmcid,
                "pmid": pmid,
                "gene": self.gene_name,
                "status": "api_error",
                "files_processed": uploaded_filenames,
                "external_dataset_links": [],
                "error": str(exc),
            }
        finally:
            self._cleanup_files(uploaded_files)

    def print_stats(self):
        logging.info(f"\n{'=' * 60}")
        logging.info("FINAL STATISTICS")
        logging.info(f"{'=' * 60}")
        logging.info(f"Total articles processed: {self.stats['total_articles']}")
        logging.info(f"Successful: {self.stats['successful_articles']}")
        logging.info(f"Failed: {self.stats['failed_articles']}")
        logging.info(f"External links found: {self.stats['links_found']}")
        logging.info(f"API errors: {self.stats['api_errors']}")


def load_article_list(args, db, base_dir):
    if args.pmcid:
        return [{"pmcid": args.pmcid, "pmid": db.pmid_for_pmcid(args.pmcid)}]
    if args.pmid:
        return [{"pmcid": db.pmcid_for_pmid(args.pmid), "pmid": args.pmid}]
    if args.article_csv:
        if pd is None:
            raise ImportError("CSV 입력에는 pandas가 필요합니다. `pip install pandas` 후 다시 실행하세요.")
        df = pd.read_csv(args.article_csv)
        article_list = []
        downloaded = df["downloaded"] == "Yes" if "downloaded" in df.columns else [True] * len(df)
        for _, row in df[downloaded].iterrows():
            pmcid = "" if "pmcid" not in row or pd.isna(row.get("pmcid")) else str(row.get("pmcid")).strip()
            pmid = "" if "pmid" not in row or pd.isna(row.get("pmid")) else str(row.get("pmid")).strip()
            if pmcid or pmid:
                article_list.append({"pmcid": pmcid, "pmid": pmid})
        return article_list
    if args.all:
        return [{"pmcid": "", "pmid": d.name} for d in sorted(base_dir.iterdir()) if d.is_dir() and d.name.isdigit()]
    raise ValueError("--pmcid, --pmid, --all, 또는 --article-csv 중 하나 필요")


def main():
    parser = argparse.ArgumentParser(description="논문에서 외부 환자/변이 dataset 링크만 추출")
    parser.add_argument("--pmcid", help="단일 PMCID 처리")
    parser.add_argument("--pmid", help="단일 PMID 처리")
    parser.add_argument("--all", action="store_true", help="전체 논문 폴더 배치 처리")
    parser.add_argument("--gene", required=True, help="유전자명")
    parser.add_argument("--model", default=DEFAULT_MODEL, help=f"OpenAI 모델 (기본: {DEFAULT_MODEL})")
    parser.add_argument("--input_dir", default="article_files", help="논문 폴더가 있는 디렉토리")
    parser.add_argument("--db-dir", default="databases", help="유전자별 SQLite DB 저장 디렉토리")
    parser.add_argument("--article-csv", default=None, help="처리할 논문 목록 CSV")
    parser.add_argument("--max-articles", type=int, default=None, help="최대 처리 개수")
    parser.add_argument("--force", action="store_true", help="이미 링크가 있는 논문도 다시 처리")
    parser.add_argument("--dry-run", action="store_true", help="DB에 저장하지 않고 결과만 로그 출력")
    args = parser.parse_args()

    base_dir = WORKSPACE / args.input_dir
    db_path = WORKSPACE / args.db_dir / f"{args.gene.upper()}.sqlite"

    link_db = ExternalDatasetLinkDB(db_path)
    extractor = ExternalDatasetLinkExtractor(
        api_key=OPENAI_API_KEY,
        gene_name=args.gene,
        model=args.model,
    )

    logging.info("=" * 60)
    logging.info("External Dataset Link Extraction")
    logging.info(f"Gene: {args.gene}")
    logging.info(f"Model: {args.model}")
    logging.info(f"Input: {base_dir}")
    logging.info(f"DB: {db_path}")
    logging.info("=" * 60)

    def article_dir_for(article):
        if article["pmid"]:
            return base_dir / article["pmid"]
        return base_dir / article["pmcid"]

    try:
        article_list = load_article_list(args, link_db, base_dir)
        if args.max_articles:
            article_list = article_list[: args.max_articles]
        logging.info(f"총 {len(article_list)}개 논문 처리 예정")

        for i, article in enumerate(article_list, 1):
            pmcid = article["pmcid"]
            pmid = article["pmid"]
            article_label = pmcid or f"PMID:{pmid}"
            logging.info(f"\n[{i}/{len(article_list)}] {article_label}")

            if not args.force and link_db.has_links_for_article(pmcid=pmcid, pmid=pmid):
                logging.info(f"{article_label}: external_dataset_links 기존 row 있음 (스킵)")
                continue

            article_dir = article_dir_for(article)
            link_db.refresh_article_files(article_dir, pmcid, pmid)
            result = extractor.extract_from_article(article_dir, pmcid=pmcid, pmid=pmid)

            if result.get("status") == "success":
                if args.dry_run:
                    logging.info(json.dumps(result.get("external_dataset_links", []), ensure_ascii=False, indent=2))
                else:
                    saved = link_db.upsert_links(result)
                    logging.info(f"DB 저장 완료: {saved}개")

            if i < len(article_list):
                import time

                time.sleep(1)

    except KeyboardInterrupt:
        logging.warning("사용자 중단")
    except FatalApiQuotaError as exc:
        logging.error(str(exc))
    finally:
        extractor.print_stats()
        link_db.close()


if __name__ == "__main__":
    main()
