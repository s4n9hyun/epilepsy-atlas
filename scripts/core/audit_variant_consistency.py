#!/usr/bin/env python3
"""추출된 변이 표기를 MANE transcript 서열로 검사한다 (API 비용 없음, 자동 수정 없음).

검사 항목 (행마다 flag, 여러 개면 ';'로 연결):
  ref_mismatch        cDNA 치환의 ref 염기가 MANE CDS와 다름 (오타·다른 transcript·다른 유전자)
  protein_mismatch    cDNA를 번역한 아미노산이 variant_p와 다름 (위치 또는 결과)
  legacy_numbering    MANE에서는 안 맞지만 옛 transcript(예: SCN1A NM_006920)에서는 cDNA·protein이 맞음
                      → 값이 옛 번호 그대로 저장됨, MANE 변환 필요
  beyond_cds          cDNA 위치가 MANE CDS 길이를 넘음
  other_gene          변이 필드나 근거 문장에 다른 대상 유전자 변이가 적혀 있음
  legacy_unconverted  논문 transcript가 MANE와 다르고 변환이 실패했으며 MANE 서열로도 확인되지 않음

서열은 NCBI efetch(fasta_cds_na)로 받아 databases/reference_cds/에 캐시한다.

사용:
  python3 scripts/core/audit_variant_consistency.py --gene SCN2A --output-dir runs/YYYYMMDD_variant_qc
"""
import argparse
import csv
import re
import sqlite3
import sys
import time
import urllib.request
from collections import Counter
from datetime import datetime
from pathlib import Path

WORKSPACE_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(Path(__file__).resolve().parent))
from batch_extract_patient_mentions import AA_1_TO_3, GNOMAD_TRANSCRIPTS  # noqa: E402

CDS_CACHE = WORKSPACE_ROOT / "databases" / "reference_cds"
# 논문에서 자주 쓰는 옛 transcript (MANE에서 안 맞을 때 대조)
LEGACY_TRANSCRIPTS = {
    "SCN1A": ["NM_006920", "NM_001202435"],
    "SCN2A": ["NM_021007"],
    "KCNQ2": ["NM_004518"],
    "CDKL5": ["NM_003159"],
    "SYNGAP1": ["NM_001130066"],
    "STXBP1": ["NM_003165"],
    "KCNT1": ["NM_001272003"],
    "SCN8A": ["NM_014191", "NM_001177984"],
    "PCDH19": ["NM_001105243", "NM_020766"],
    "SMC1A": [],
    "FOXG1": [],
    "GABRA1": ["NM_000806"],
    "EEF1A2": [],
}
CONVERTED = {"already_gnomad", "assumed_gnomad", "converted", "mutalyzer_converted", "manual_corrected"}
CODON_TABLE = {
    a + b + c: aa
    for (a, b, c), aa in zip(
        ((x, y, z) for x in "TCAG" for y in "TCAG" for z in "TCAG"),
        "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    )
}
C_SUB = re.compile(r"c\.(\d+)([ACGT])>([ACGT])$", re.I)
C_POS = re.compile(r"c\.(\d+)")
P_SIMPLE = re.compile(r"p\.\(?([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2}|Ter|=)\)?$")


def fetch_cds(accession):
    CDS_CACHE.mkdir(parents=True, exist_ok=True)
    path = CDS_CACHE / f"{accession}.fa"
    if not path.exists():
        url = ("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore"
               f"&id={accession}&rettype=fasta_cds_na&retmode=text")
        with urllib.request.urlopen(url, timeout=60) as resp:
            path.write_text(resp.read().decode())
        time.sleep(0.4)
    lines = path.read_text().splitlines()
    return "".join(line.strip() for line in lines if not line.startswith(">")).upper()


def aa3(one):
    return AA_1_TO_3.get(one, "?")


def cdna(row):
    for value in (row["variant_c_gnomad"], row["variant_c_original"]):
        value = str(value or "").strip()
        if "c." in value:
            return "c." + value.split("c.", 1)[1].replace(" ", "")
    return ""


def check_substitution(cds, pos, ref, alt, p):
    """치환 변이를 한 CDS에서 확인한다. (ref 일치, protein 일치 또는 None, 번역 결과)"""
    if pos > len(cds):
        return False, None, ""
    codon_no = (pos - 1) // 3 + 1
    codon = cds[(codon_no - 1) * 3: codon_no * 3]
    mutant = codon[:(pos - 1) % 3] + alt + codon[(pos - 1) % 3 + 1:]
    expected_ref, expected_alt = aa3(CODON_TABLE.get(codon, "?")), aa3(CODON_TABLE.get(mutant, "?"))
    if expected_ref == expected_alt:
        expected_alt = "="
    expected = f"p.{expected_ref}{codon_no}{expected_alt}"
    pm = P_SIMPLE.match(p)
    if not pm:
        return cds[pos - 1] == ref, None, expected
    p_ref, p_pos, p_alt = pm.group(1), int(pm.group(2)), pm.group(3)
    same = (p_ref, p_pos) == (expected_ref, codon_no) and (
        p_alt == expected_alt or (expected_alt == "=" and p_alt == expected_ref))
    return cds[pos - 1] == ref, same, expected


def check_row(row, gene, cds, legacy_cds, gene_pattern):
    flags, notes = [], []
    c = cdna(row)
    p = str(row["variant_p"] or "").strip()
    sub = C_SUB.match(c)
    pos_match = C_POS.match(c)
    verified = False
    if sub:
        pos, ref, alt = int(sub.group(1)), sub.group(2).upper(), sub.group(3).upper()
        ref_ok, protein_ok, expected = check_substitution(cds, pos, ref, alt, p)
        verified = ref_ok and protein_ok is not False
        if not verified:
            legacy_hit = next((acc for acc, seq in legacy_cds.items()
                               if check_substitution(seq, pos, ref, alt, p)[:2] in ((True, True), (True, None))), None)
            if pos > len(cds):
                flags.append("beyond_cds")
                notes.append(f"c.{pos} > CDS {len(cds)}")
            elif legacy_hit and (protein_ok is False or not ref_ok):
                flags.append("legacy_numbering")
                notes.append(f"fits {legacy_hit}, not MANE (MANE reading {expected})")
            else:
                if not ref_ok:
                    flags.append("ref_mismatch")
                    notes.append(f"MANE c.{pos} is {cds[pos - 1]}")
                if protein_ok is False:
                    flags.append("protein_mismatch")
                    notes.append(f"cDNA gives {expected}")
    elif pos_match and int(pos_match.group(1)) > len(cds):
        flags.append("beyond_cds")
        notes.append(f"c.{pos_match.group(1)} > CDS {len(cds)}")
    text = " ".join(str(row[k] or "") for k in ("variant_c_original", "variant_p", "location_in_paper"))
    other = gene_pattern.search(text)
    if other and other.group(0) != gene:
        flags.append("other_gene")
        notes.append(f"mentions {other.group(0)}")
    stated = str(row["nm_transcript"] or "").split(".")[0]
    mane = GNOMAD_TRANSCRIPTS[gene].split(".")[0]
    if (stated.startswith("NM_") and stated != mane and not verified and "legacy_numbering" not in flags
            and str(row["gnomad_conversion_status"] or "") not in CONVERTED):
        flags.append("legacy_unconverted")
        notes.append(f"paper transcript {row['nm_transcript']}, status {row['gnomad_conversion_status']}")
    return flags, notes


def main():
    parser = argparse.ArgumentParser(description="변이 표기 MANE 서열 검사")
    parser.add_argument("--gene", nargs="+", default=list(GNOMAD_TRANSCRIPTS))
    parser.add_argument("--db-dir", default=str(WORKSPACE_ROOT / "databases"))
    parser.add_argument("--output-dir", default=str(WORKSPACE_ROOT / "runs" / f"{datetime.now():%Y%m%d}_variant_qc"))
    parser.add_argument("--include-excluded", action="store_true", help="manual_exclude 행도 검사")
    args = parser.parse_args()

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / f"variant_consistency_{datetime.now():%Y%m%d_%H%M%S}.csv"
    gene_pattern = re.compile(r"\b(" + "|".join(GNOMAD_TRANSCRIPTS) + r"|CNTNAP2|SCN3A|SCN9A|KCNQ3|MECP2)\b")
    fields = ["gene", "mention_id", "pmid", "patient_id", "mention_source_type", "variant_c", "variant_p",
              "nm_transcript", "gnomad_conversion_status", "flags", "notes"]
    with out_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        for gene in [g.upper() for g in args.gene]:
            db = Path(args.db_dir) / f"{gene}.sqlite"
            if not db.exists():
                continue
            conn = sqlite3.connect(db)
            conn.row_factory = sqlite3.Row
            if not conn.execute("SELECT 1 FROM sqlite_master WHERE name='patient_mentions'").fetchone():
                continue
            cds = fetch_cds(GNOMAD_TRANSCRIPTS[gene])
            legacy_cds = {acc: fetch_cds(acc) for acc in LEGACY_TRANSCRIPTS.get(gene, [])}
            counts = Counter()
            rows = conn.execute("SELECT * FROM patient_mentions").fetchall()
            for row in rows:
                if row["mention_source_type"] == "manual_exclude" and not args.include_excluded:
                    continue
                flags, notes = check_row(row, gene, cds, legacy_cds, gene_pattern)
                counts.update(flags)
                if flags:
                    writer.writerow({
                        "gene": gene, "mention_id": row["mention_id"], "pmid": row["pmid"],
                        "patient_id": row["patient_id"], "mention_source_type": row["mention_source_type"],
                        "variant_c": cdna(row), "variant_p": row["variant_p"], "nm_transcript": row["nm_transcript"],
                        "gnomad_conversion_status": row["gnomad_conversion_status"],
                        "flags": ";".join(flags), "notes": "; ".join(notes),
                    })
            print(f"{gene}: rows {len(rows)}, CDS {len(cds)} nt, flags {dict(counts)}")
    print(f"-> {out_path}")


if __name__ == "__main__":
    main()
