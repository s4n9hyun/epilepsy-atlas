#!/usr/bin/env python3
"""논문 파일에 적힌 변이와 DB에 추출된 변이를 비교해 덜 추출된 논문을 찾는다 (API 비용 없음).

큰 표·긴 supplement에서 일부 행만 추출되거나(예: 148행 중 29행), ZIP 등으로 파일을 못 읽은 논문을 잡는다.
- 논문 폴더의 모든 문서(본문, supplement, ZIP 안 문서, 확장자 없는 파일)에서 텍스트를 뽑는다.
- 대상 유전자 기호가 있는 줄(다른 대상 유전자가 없는 문서는 전체)에서 c./p. 표기를 찾는다.
- 같은 PMID의 DB 변이(모든 mention)와 비교해 없는 변이 수를 센다.
- 없는 변이가 많은 논문을 `undercount_candidate`로 표시한다. 변이 목록 논문(기능 연구·총설)도 걸리므로
  재추출 전 사람이 확인한다.

사용:
  python3 scripts/core/audit_supplement_coverage.py --gene STXBP1 --output-dir runs/YYYYMMDD_coverage
  python3 scripts/core/audit_supplement_coverage.py --gene SCN1A --pmid 17347258 23808377
"""
import argparse
import csv
import re
import sqlite3
import subprocess
import sys
import tempfile
import zipfile
from datetime import datetime
from pathlib import Path

WORKSPACE_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(Path(__file__).resolve().parent))
from audit_variant_consistency import GNOMAD_TRANSCRIPTS, fetch_cds  # noqa: E402

GENES = ["SCN1A", "SCN2A", "KCNQ2", "CDKL5", "SYNGAP1", "EEF1A2",
         "PCDH19", "SCN8A", "GABRA1", "STXBP1", "SMC1A", "FOXG1", "KCNT1"]
DOC_EXT = {".pdf", ".doc", ".docx", ".xls", ".xlsx", ".csv", ".tsv", ".txt", ".rtf", ".html", ".htm", ".nxml", ".xml"}
SIG = {b"PK\x03\x04": "zip", b"%PDF": "pdf", b"\xd0\xcf\x11\xe0": "ole"}
AA3 = {"A": "ala", "R": "arg", "N": "asn", "D": "asp", "C": "cys", "Q": "gln", "E": "glu", "G": "gly", "H": "his",
       "I": "ile", "L": "leu", "K": "lys", "M": "met", "F": "phe", "P": "pro", "S": "ser", "T": "thr", "W": "trp",
       "Y": "tyr", "V": "val", "X": "ter", "*": "ter"}
C_RE = re.compile(r"c\.\s?[-*]?\d+(?:[+-]\d+)?(?:_[-*]?\d+(?:[+-]\d+)?)?\s?"
                  r"(?:[ACGT]\s?>\s?[ACGT]|del[ACGT]*|dup[ACGT]*|ins[ACGT]+|delins[ACGT]+)", re.I)
P3_RE = re.compile(r"p\.\(?([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2}|\*|Ter|fs|=)", re.I)
P1_RE = re.compile(r"p\.\(?([ACDEFGHIKLMNPQRSTVWY])(\d+)([ACDEFGHIKLMNPQRSTVWYX*]|fs)\b")
# 표에서 p. 없이 쓴 단백질 표기 (예: R1882Q). 오탐을 줄이려고 위치 2~4자리만 본다.
BARE_P_RE = re.compile(r"(?<![\w.])([ACDEFGHIKLMNPQRSTVWY])(\d{2,4})([ACDEFGHIKLMNPQRSTVWYX*])(?![\w])")
UNDERCOUNT_MIN_MISSING = 5
UNDERCOUNT_MIN_RATIO = 0.3


def kind(data):
    return next((v for s, v in SIG.items() if data.startswith(s)), None)


def norm_c(value):
    return re.sub(r"\s", "", value).lower()


def norm_p(aa1, pos, aa2):
    a = AA3.get(aa1.upper(), aa1.lower()) if len(aa1) == 1 else aa1.lower()
    b = AA3.get(aa2.upper(), aa2.lower()) if len(aa2) == 1 else aa2.lower()
    if b in ("*", "x"):
        b = "ter"
    return f"p.{a}{pos}{'fs' if b.startswith('fs') else b}"


def text_of(name, data):
    """파일 종류별로 텍스트를 뽑는다 (pdftotext / openpyxl / xlrd / textutil)."""
    ext = Path(name).suffix.lower()
    k = kind(data[:4])
    with tempfile.TemporaryDirectory() as tmp:
        path = Path(tmp) / ("f" + (ext or ".bin"))
        path.write_bytes(data)
        try:
            if ext == ".pdf" or k == "pdf":
                return subprocess.run(["pdftotext", "-layout", str(path), "-"], capture_output=True,
                                      text=True, timeout=120).stdout
            if ext == ".xlsx" or (k == "zip" and ext == ""):
                try:
                    import openpyxl
                    wb = openpyxl.load_workbook(path, read_only=True, data_only=True)
                    return "\n".join("\t".join("" if c is None else str(c) for c in row)
                                     for ws in wb.worksheets for row in ws.iter_rows(values_only=True))
                except Exception:
                    pass
            if ext == ".xls" or (k == "ole" and ext == ""):
                try:
                    import xlrd
                    wb = xlrd.open_workbook(str(path))
                    return "\n".join("\t".join(str(v) for v in wb.sheet_by_index(i).row_values(r))
                                     for i in range(wb.nsheets) for r in range(wb.sheet_by_index(i).nrows))
                except Exception:
                    pass
            if ext in (".doc", ".docx", ".rtf") or k in ("ole", "zip"):
                out = Path(tmp) / "out.txt"
                target = path.with_suffix(ext if ext in (".doc", ".docx", ".rtf") else (".doc" if k == "ole" else ".docx"))
                path.rename(target)
                subprocess.run(["textutil", "-convert", "txt", str(target), "-output", str(out)],
                               capture_output=True, timeout=120)
                return out.read_text(errors="ignore") if out.exists() else ""
            if ext in (".csv", ".tsv", ".txt", ".html", ".htm", ".nxml", ".xml"):
                raw = data.decode("utf-8", errors="ignore")
                return re.sub(r"<[^>]+>", "\n", raw) if ext in (".html", ".htm", ".nxml", ".xml") else raw
        except Exception:
            return ""
    return ""


def article_documents(article_dir):
    """폴더의 문서와 ZIP 안 문서를 (이름, bytes)로 돌려준다. 폴더에 이미 풀린 같은 이름은 한 번만."""
    files = [f for f in sorted(article_dir.iterdir()) if f.is_file()]
    loose = {f.name.lower() for f in files}
    for f in files:
        data = f.read_bytes()
        k = kind(data[:4])
        if k == "zip" and f.suffix.lower() in ("", ".zip"):
            try:
                archive = zipfile.ZipFile(f)
                names = archive.namelist()
            except Exception:
                continue
            if "[Content_Types].xml" in names:
                yield f.name, data
                continue
            for n in names:
                if (Path(n).suffix.lower() in DOC_EXT and Path(n).name.lower() not in loose
                        and "__MACOSX" not in n):
                    yield f"{f.name}:{n}", archive.read(n)
        elif f.suffix.lower() in DOC_EXT or (f.suffix == "" and k in ("pdf", "ole")):
            if f.name.lower().startswith("pubmed_"):
                continue
            yield f.name, data


def variants_for_gene(text, gene, max_aa=None):
    others = [g for g in GENES if g != gene and re.search(rf"\b{g}\b", text)]
    lines = text.split("\n") if others else [text]
    found = set()
    for line in lines:
        if others and not re.search(rf"\b{gene}\b", line):
            continue
        found |= {("c", norm_c(m.group(0))) for m in C_RE.finditer(line)}
        found |= {("p", norm_p(*m.groups())) for m in P3_RE.finditer(line)}
        found |= {("p", norm_p(*m.groups())) for m in P1_RE.finditer(line)}
        found |= {("p", norm_p(*m.groups())) for m in BARE_P_RE.finditer(line)}
    if max_aa:
        # 단백질 길이를 넘는 번호는 다른 열(염기 번호 등)에서 잘못 잡힌 것
        found = {(t, v) for t, v in found if t != "p" or int(re.search(r"\d+", v).group(0)) <= max_aa}
    return found


def db_variants(conn, pmid):
    have, mentions = set(), 0
    for c1, c2, p in conn.execute(
            "SELECT variant_c_original, variant_c_gnomad, variant_p FROM patient_mentions WHERE pmid=?", (pmid,)):
        mentions += 1
        for c in (c1, c2):
            have |= {("c", norm_c(m.group(0))) for m in C_RE.finditer(str(c or ""))}
        have |= {("p", norm_p(*m.groups())) for m in P3_RE.finditer(str(p or ""))}
    return have, mentions


def main():
    parser = argparse.ArgumentParser(description="논문 파일 변이 vs DB 변이 대조")
    parser.add_argument("--gene", nargs="+", required=True)
    parser.add_argument("--pmid", nargs="*", help="특정 PMID만")
    parser.add_argument("--article-root", default=str(WORKSPACE_ROOT / "article_files"))
    parser.add_argument("--output-dir", default=str(WORKSPACE_ROOT / "runs" / f"{datetime.now():%Y%m%d}_coverage"))
    args = parser.parse_args()

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / f"supplement_coverage_{datetime.now():%Y%m%d_%H%M%S}.csv"
    fields = ["gene", "pmid", "documents", "text_chars", "variants_in_files", "db_mentions",
              "missing_c", "missing_p", "undercount_candidate", "missing_examples"]
    with out_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        for gene in [g.upper() for g in args.gene]:
            try:
                max_aa = len(fetch_cds(GNOMAD_TRANSCRIPTS[gene])) // 3
            except Exception:
                max_aa = None
            conn = sqlite3.connect(WORKSPACE_ROOT / "databases" / f"{gene}.sqlite")
            has_mentions = conn.execute("SELECT 1 FROM sqlite_master WHERE name='patient_mentions'").fetchone()
            pmids = args.pmid or [str(r[0]) for r in conn.execute("SELECT pmid FROM pubmed_articles")]
            flagged = 0
            for pmid in pmids:
                article_dir = Path(args.article_root) / pmid
                if not article_dir.is_dir():
                    continue
                texts = [text_of(name, data) for name, data in article_documents(article_dir)]
                found = set().union(*(variants_for_gene(t, gene, max_aa) for t in texts)) if texts else set()
                have, mentions = db_variants(conn, pmid) if has_mentions else (set(), 0)
                missing_c = sorted(v for t, v in found if t == "c" and ("c", v) not in have)
                missing_p = sorted(v for t, v in found if t == "p" and ("p", v) not in have)
                missing = max(len(missing_c), len(missing_p))
                candidate = (missing >= UNDERCOUNT_MIN_MISSING
                             and missing >= UNDERCOUNT_MIN_RATIO * max(len(found) / 2, 1))
                flagged += candidate
                writer.writerow({
                    "gene": gene, "pmid": pmid, "documents": len(texts), "text_chars": sum(map(len, texts)),
                    "variants_in_files": len(found), "db_mentions": mentions,
                    "missing_c": len(missing_c), "missing_p": len(missing_p),
                    "undercount_candidate": "yes" if candidate else "",
                    "missing_examples": "; ".join((missing_c + missing_p)[:12]),
                })
            print(f"{gene}: articles {len(pmids)}, undercount candidates {flagged}")
    print(f"-> {out_path}")


if __name__ == "__main__":
    main()
