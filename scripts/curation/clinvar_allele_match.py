"""Allele-level ClinVar matching of atlas variants (reviewer comments R2-M7, R3-3, R4-M2).

1. Every distinct small variant in the canonical datasets (HGVS c. on MANE or on the alternative transcript of an exon absent
   from MANE) is converted to a GRCh38 genomic allele with NCBI Variation Services (hgvs -> contextual SPDI), cached.
2. ClinVar alleles come from the frozen variant_summary release (runs/paper/clinvar_snapshot; GRCh38 VCF fields).
3. Both are normalised the same way (trim shared bases, then left-shift indels on the GRCh38 reference of the gene region)
   and matched on (chromosome, position, ref, alt).
Writes runs/review_fixes/clinvar/allele_map.csv: gene, hgvs, chrom, pos, ref, alt, vcv, significance, review_status, status.
"""
import json
import re
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

import pandas as pd

W = Path(__import__("os").environ.get("EPILEPSY_WORKSPACE", Path(__file__).resolve().parents[2]))
sys.path.insert(0, str(W / "paper/analysis"))
from common import GENES, MANE, load_rows  # noqa: E402

OUT = W / "runs/review_fixes/clinvar"
SNAP = W / "runs/paper/clinvar_snapshot"
ALT_OK = {"SCN8A": "NM_014191.4", "CDKL5": "NM_003159.3", "SCN2A": "NM_021007.3"}
HGVS_C = re.compile(r"^c\.[-*]?\d+(?:[+-]\d+)?(?:_[-*]?\d+(?:[+-]\d+)?)?(?:[ACGT]>[ACGT]|del[ACGT]*(?:ins[ACGT]+)?|dup[ACGT]*|ins[ACGT]+|delins[ACGT]+)$")
CHR_ACC = {"1": "NC_000001.11", "2": "NC_000002.12", "6": "NC_000006.12", "9": "NC_000009.12", "12": "NC_000012.12",
           "14": "NC_000014.9", "20": "NC_000020.11", "22": "NC_000022.11", "X": "NC_000023.11", "4": "NC_000004.12", "5": "NC_000005.10"}


def get(url, tries=4):
    for k in range(tries):
        try:
            with urllib.request.urlopen(url, timeout=60) as r:
                return json.loads(r.read().decode())
        except urllib.error.HTTPError as e:
            if e.code in (400, 404, 422):
                return {"error": e.code}
            time.sleep(2 + 3 * k)
        except Exception:
            time.sleep(2 + 3 * k)
    return {"error": "failed"}


def to_spdi(hgvs, cache):
    """HGVS c. -> GRCh38 canonical SPDI [seq_id, 0-based position, deleted, inserted] (None if invalid)."""
    if hgvs in cache:
        return cache[hgvs]
    j = get("https://api.ncbi.nlm.nih.gov/variation/v0/hgvs/" + urllib.parse.quote(hgvs, safe="") + "/contextuals")
    time.sleep(0.3)
    sp = None
    try:
        t = j["data"]["spdis"][0]
        spdi = f"{t['seq_id']}:{t['position']}:{t['deleted_sequence']}:{t['inserted_sequence']}"
        c = get("https://api.ncbi.nlm.nih.gov/variation/v0/spdi/" + urllib.parse.quote(spdi, safe="") + "/canonical_representative")
        time.sleep(0.3)
        d = c["data"]
        sp = [d["seq_id"], d["position"], d["deleted_sequence"], d["inserted_sequence"]]
    except Exception:
        sp = None
    cache[hgvs] = sp
    return sp


def region_seq(acc, start, stop, cache_dir):
    f = cache_dir / f"{acc}_{start}_{stop}.txt"
    if not f.exists():
        url = (f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={acc}&rettype=fasta&retmode=text"
               f"&seq_start={start}&seq_stop={stop}")
        with urllib.request.urlopen(url, timeout=120) as r:
            f.write_text("".join(x.strip() for x in r.read().decode().splitlines() if not x.startswith(">")).upper())
        time.sleep(0.4)
    return f.read_text()


class Ref:
    """GRCh38 sequence windows around the study genes (1-based coordinates)."""

    def __init__(self, windows, cache_dir):
        self.w = {k: (s, region_seq(k[0], s, e, cache_dir)) for k, (s, e) in windows.items()}

    def base(self, acc, pos):
        for (a, _), (s, seq) in self.w.items():
            if a == acc and s <= pos < s + len(seq):
                return seq[pos - s]
        return None


def normalise(acc, pos, ref, alt, R):
    """VCF-style (1-based pos, ref, alt, possibly empty) -> trimmed, left-shifted representation."""
    while ref and alt and ref[-1] == alt[-1]:
        ref, alt = ref[:-1], alt[:-1]
    while ref and alt and ref[0] == alt[0]:
        ref, alt, pos = ref[1:], alt[1:], pos + 1
    if ref and alt:
        return (acc, pos, ref, alt)
    if not ref and not alt:  # no change (reference-only record)
        return (acc, pos, "", "")
    seq = ref or alt
    # left-shift a pure insertion/deletion: the event sits before position pos
    for _ in range(500):
        prev = R.base(acc, pos - 1)
        if prev is None or prev != seq[-1]:
            break
        seq = prev + seq[:-1]
        pos -= 1
    return (acc, pos, seq if ref else "", seq if alt else "")


def status(sig):
    s = str(sig).lower()
    if "conflicting" in s:
        return "Conflicting"
    if "pathogenic" in s and "benign" not in s and "uncertain" not in s:
        return "Pathogenic/likely pathogenic"
    if "uncertain" in s:
        return "Uncertain significance"
    if "benign" in s and "pathogenic" not in s:
        return "Benign/likely benign"
    return "Other"


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "ref").mkdir(exist_ok=True)
    rows = load_rows()
    hg = set()
    for g, tr, c in zip(rows.gene, rows.gnomad_transcript, rows.variant_c):
        c = re.sub(r"\s", "", str(c))
        if not HGVS_C.match(c):
            continue
        if tr.split(".")[0] == MANE[g].split(".")[0] or tr == ALT_OK.get(g):
            hg.add((g, f"{tr}:{c}"))
    cache_f = OUT / "spdi_cache.json"
    cache = json.loads(cache_f.read_text()) if cache_f.exists() else {}
    for k, (g, h) in enumerate(sorted(hg)):
        to_spdi(h, cache)
        if k % 200 == 0:
            cache_f.write_text(json.dumps(cache))
            print(k, len(hg), flush=True)
    cache_f.write_text(json.dumps(cache))
    vs = pd.read_csv(SNAP / "variant_summary_13genes.tsv", sep="\t", dtype=str, keep_default_na=False)
    sig_col = "ClinicalSignificance" if "ClinicalSignificance" in vs.columns else "GermlineClassification"
    vs = vs[(vs.Assembly == "GRCh38") & (vs.PositionVCF.str.match(r"^\d+$"))].copy()
    vs["acc"] = vs.Chromosome.map(CHR_ACC)
    vs = vs[vs.acc.notna()]
    # windows covering atlas alleles and ClinVar alleles, per chromosome
    pos = {}
    for (g, h) in hg:
        sp = cache.get(h)
        if sp:
            pos.setdefault(sp[0], []).append(sp[1])
    for acc, p in zip(vs.acc, vs.PositionVCF.astype(int)):
        if acc in pos:
            pos[acc].append(p)
    windows = {}
    for acc, ps in pos.items():
        ps = sorted(ps)
        start = ps[0]
        prev = ps[0]
        for p in ps[1:] + [None]:
            if p is None or p - prev > 200000:
                windows[(acc, start)] = (max(1, start - 2000), prev + 2000)
                start = p
            prev = p if p is not None else prev
    R = Ref(windows, OUT / "ref")
    cv = {}
    for acc, p, rf, al, vid, sig, rev in zip(vs.acc, vs.PositionVCF.astype(int), vs.ReferenceAlleleVCF, vs.AlternateAlleleVCF,
                                            vs.VariationID, vs[sig_col], vs.ReviewStatus):
        if not re.match(r"^[ACGTN]+$", rf) or not re.match(r"^[ACGTN]+$", al):
            continue
        cv.setdefault(normalise(acc, p, rf, al, R), []).append((vid, sig, rev))
    out = []
    for (g, h) in sorted(hg):
        sp = cache.get(h)
        if not sp:
            out.append(dict(gene=g, hgvs=h, chrom="", pos="", ref="", alt="", vcv="", significance="", review_status="",
                            status="Unmapped"))
            continue
        acc, p0, d, ins = sp
        # SPDI: 0-based position, deleted, inserted -> VCF-like 1-based with possibly empty alleles
        key = normalise(acc, p0 + 1, d, ins, R)
        hits = cv.get(key, [])
        if hits:
            vids = sorted({x[0] for x in hits})
            sigs = sorted({x[1] for x in hits})
            revs = sorted({x[2] for x in hits})
            st = status(" | ".join(sigs))
        else:
            vids, sigs, revs, st = [], [], [], "Not in ClinVar"
        out.append(dict(gene=g, hgvs=h, chrom=acc, pos=key[1], ref=key[2], alt=key[3], vcv=";".join(vids),
                        significance=" | ".join(sigs), review_status=" | ".join(revs), status=st))
    o = pd.DataFrame(out)
    o.to_csv(OUT / "allele_map.csv", index=False)
    print(o.status.value_counts().to_string())


if __name__ == "__main__":
    main()
