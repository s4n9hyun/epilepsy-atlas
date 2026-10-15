"""Reference consistency of every atlas variant on the MANE Select transcript, with estimated corrections for source typos.

Clinician request 2026-10-12 (runs/1012/mail.md): when the reference amino acid (or base) printed in the article does not match the
reference sequence, estimate where the typo is, correct to a reference-consistent variant and mark the value as estimated.

Checks per canonical row (MANE transcript only; exons absent from MANE and CNVs skipped):
  - cDNA substitution: reference base at the coding position; predicted protein change.
  - protein change (missense/nonsense/synonymous): reference amino acid at the position.
Estimation (single candidate required, otherwise left unresolved for manual review):
  1. cDNA consistent, printed protein differs  -> protein typo; protein taken from the cDNA.
  2. protein consistent, cDNA reference base wrong -> cDNA candidates at positions one digit edit away (substitution, transposition,
     +/-1 codon base) that give the printed protein change.
  3. protein reference residue wrong, no usable cDNA -> positions one digit edit away (or legacy SCN1A NM_006920 offset +11) whose
     residue is the printed reference residue.
  4. cDNA and protein each match the reference but disagree -> single-base cDNA edits within the printed codon that give the printed
     protein; if exactly one, the cDNA is the typo.
Writes runs/1012/reference_check/{all_rows,flagged,proposals}.csv
"""
import itertools
import re
import sys
from pathlib import Path

import pandas as pd

W = Path(__import__("os").environ.get("EPILEPSY_WORKSPACE", Path(__file__).resolve().parents[2]))
sys.path.insert(0, str(W / "paper/analysis"))
from common import GENES, MANE, fetch_cds, load_rows  # noqa: E402

OUT = W / "runs/1012/reference_check"
CODON = {a + b + c: aa for (a, b, c), aa in zip(itertools.product("TCAG", repeat=3),
                                                 "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG")}
AA3 = {"Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C", "Gln": "Q", "Glu": "E", "Gly": "G", "His": "H", "Ile": "I",
       "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P", "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V",
       "Ter": "*", "Stop": "*", "X": "*"}
ONE3 = {v: k for k, v in AA3.items() if k not in ("Stop", "X")}
COMP = {"A": "T", "C": "G", "G": "C", "T": "A"}


def parse_p(p):
    s = re.sub(r"\s", "", str(p))
    s = re.sub(r"^p\.?\(?", "", s).rstrip(")")
    m = re.match(r"^([A-Z][a-z]{2}|[A-Z*])(\d+)([A-Z][a-z]{2}|[A-Z*]|=)$", s)
    if not m:
        return None
    conv = lambda x: AA3.get(x, x if len(x) == 1 else None)
    ref, alt = conv(m.group(1)), (m.group(3) if m.group(3) == "=" else conv(m.group(3)))
    if not ref or not alt:
        return None
    return ref, int(m.group(2)), (ref if alt == "=" else alt)


def parse_c(c):
    s = re.sub(r"\s", "", str(c)).split(":")[-1]
    m = re.match(r"^c\.(\d+)([ACGT])>([ACGT])$", s)
    return (int(m.group(1)), m.group(2), m.group(3)) if m else None


def fmt_p(ref, pos, alt):
    return f"p.{ONE3[ref]}{pos}{'=' if alt == ref else ONE3[alt]}"


def predict(cds, pos, alt):
    i = (pos - 1) // 3
    codon = cds[i * 3:i * 3 + 3]
    new = codon[:(pos - 1) % 3] + alt + codon[(pos - 1) % 3 + 1:]
    return CODON.get(codon), i + 1, CODON.get(new)


def digit_neighbours(n):
    s = str(n)
    out = set()
    for i in range(len(s)):
        for d in "0123456789":
            if d != s[i]:
                out.add(s[:i] + d + s[i + 1:])
        if i + 1 < len(s):
            out.add(s[:i] + s[i + 1] + s[i] + s[i + 2:])
        out.add(s[:i] + s[i + 1:])
    for d in "0123456789":
        for i in range(len(s) + 1):
            out.add(s[:i] + d + s[i:])
    return {int(x) for x in out if x and not x.startswith("0")} - {n}


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    rows = load_rows()
    cds = {g: fetch_cds(MANE[g]) for g in GENES}
    prot = {g: "".join(CODON.get(s[i:i + 3], "?") for i in range(0, len(s) - 2, 3)) for g, s in cds.items()}
    recs = []
    for r in rows.itertuples():
        g = r.gene
        on_mane = str(r.gnomad_transcript).split(".")[0] == MANE[g].split(".")[0]
        if r.consequence == "CNV / exon-level" or (str(r.gnomad_transcript) not in ("", "nan") and not on_mane):
            continue
        c = parse_c(r.variant_c) if on_mane else None
        p = parse_p(r.variant_p) or parse_p(r.protein)
        rec = dict(gene=g, canonical_patient_id=r.canonical_patient_id, pmid=r.pmid, patient_id=r.patient_id,
                   supporting_mention_ids=r.supporting_mention_ids, variant_c=r.variant_c, variant_p=r.variant_p,
                   status=str(r.gnomad_conversion_status)[:80], c_ok="", p_ok="", agree="", issue="", proposal_c="", proposal_p="",
                   method="", candidates="")
        S, P = cds[g], prot[g]
        if c:
            pos, ref, alt = c
            rec["c_ok"] = pos <= len(S) and S[pos - 1] == ref
        if p:
            pref, ppos, palt = p
            rec["p_ok"] = ppos <= len(P) and P[ppos - 1] == pref
        if c and rec["c_ok"]:
            r0, ppos_c, a0 = predict(S, c[0], c[2])
            if p:
                rec["agree"] = (p[1] == ppos_c and p[2] == a0)
        if not c and not p:
            continue
        ok = (not c or rec["c_ok"]) and (not p or rec["p_ok"]) and (rec["agree"] in ("", True))
        if ok:
            recs.append(rec)
            continue
        # ---- estimation
        cands = []
        if c and rec["c_ok"] and p and rec["agree"] is False:
            pos, ref, alt = c
            r0, ppos_c, a0 = predict(S, pos, alt)
            if p[1] == ppos_c and rec["p_ok"]:
                # same codon: which single-base edit in that codon gives the printed protein?
                ed = []
                for k in range(3):
                    cp = (ppos_c - 1) * 3 + k + 1
                    for b in "ACGT":
                        if b != S[cp - 1] and predict(S, cp, b)[2] == p[2]:
                            ed.append(f"c.{cp}{S[cp - 1]}>{b}")
                rec["issue"] = "cDNA and protein disagree (same codon)"
                cands = [("c", e, fmt_p(*p)) for e in ed]
                if not ed:
                    rec["issue"] += "; no single-base cDNA gives the printed protein"
            else:
                rec["issue"] = "cDNA and protein disagree (different codon)"
                cands = [("p", f"c.{pos}{ref}>{alt}", fmt_p(r0, ppos_c, a0))]
        elif c and not rec["c_ok"]:
            pos, ref, alt = c
            rec["issue"] = "cDNA reference base mismatch"
            for q in sorted(digit_neighbours(pos)):
                if q <= len(S) and S[q - 1] == ref:
                    r0, pp, a0 = predict(S, q, alt)
                    if p and rec["p_ok"] and pp == p[1] and a0 == p[2]:
                        cands.append(("c", f"c.{q}{ref}>{alt}", fmt_p(*p)))
                    elif p and not rec["p_ok"]:
                        pass
            if not cands and p and rec["p_ok"]:
                # strand/allele swap: printed alt is the reference base
                if pos <= len(S) and S[pos - 1] == alt:
                    r0, pp, a0 = predict(S, pos, ref)
                    if pp == p[1] and a0 == p[2]:
                        cands.append(("c", f"c.{pos}{alt}>{ref}", fmt_p(*p)))
        elif p and not rec["p_ok"]:
            pref, ppos, palt = p
            rec["issue"] = "protein reference residue mismatch"
            pool = set(digit_neighbours(ppos))
            if g == "SCN1A" and ppos >= 671 - 11:
                pool.add(ppos + 11)
            for q in sorted(pool):
                if q <= len(P) and P[q - 1] == pref:
                    if c and rec["c_ok"]:
                        continue
                    cands.append(("p", "", fmt_p(pref, q, palt)))
        rec["candidates"] = " | ".join(f"{a} {b}".strip() for _, a, b in cands)
        if len(cands) == 1:
            kind, cc, pp = cands[0]
            rec["proposal_c"] = cc
            rec["proposal_p"] = pp
            rec["method"] = {"c": "cDNA typo (single candidate consistent with printed protein)",
                             "p": "protein typo (single reference-consistent candidate)"}[kind]
        recs.append(rec)
    d = pd.DataFrame(recs)
    d.to_csv(OUT / "all_rows.csv", index=False)
    f = d[d.issue != ""]
    f.to_csv(OUT / "flagged.csv", index=False)
    f[f.proposal_c.ne("") | f.proposal_p.ne("")].to_csv(OUT / "proposals.csv", index=False)
    print("checked", len(d), "flagged", len(f), "with single proposal", int((f.proposal_p != "").sum()))
    print(f.groupby(["gene", "issue"]).size().to_string())


if __name__ == "__main__":
    main()
