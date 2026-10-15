"""Sequence-based transcript QC of atlas variants (reviewer comments R3-1, R4-M1, R5-2).

For every canonical patient-variant row with a cDNA variant, the reference base and the encoded protein change are checked
on the MANE Select transcript and on the alternative transcript used by the source article. Rows are then:
- relabelled to MANE when the variant lies in a region where the alternative transcript is identical to MANE (or fits MANE
  only), so legacy-labelled variants enter MANE-based analyses;
- relabelled to the alternative transcript when the variant only fits an exon absent from MANE (SCN8A exon 5N, the 3' exons
  of CDKL5 NM_003159);
- given a "c." prefix when it was missing;
- flagged for source review when the cDNA and the stated protein change are inconsistent on every transcript;
- excluded when the variant lies outside the gene's coding sequence on every transcript (wrong gene assignment).
Writes runs/review_fixes/variant_qc/{rules.csv, report.csv, source_check_needed.csv}.
"""
import re
import sqlite3
import sys
from pathlib import Path

import pandas as pd

W = Path(__import__("os").environ.get("EPILEPSY_WORKSPACE", Path(__file__).resolve().parents[2]))
sys.path.insert(0, str(W / "scripts/core"))
sys.path.insert(0, str(W / "paper/analysis"))
from audit_variant_consistency import fetch_cds  # noqa: E402
from common import GENES, MANE, load_rows  # noqa: E402

OUT = W / "runs/review_fixes/variant_qc"
ALT = {"GABRA1": "NM_000806.5", "SCN8A": "NM_014191.4", "STXBP1": "NM_003165.6", "PCDH19": "NM_001105243.2",
       "CDKL5": "NM_003159.3", "SCN1A": "NM_006920.6", "SCN2A": "NM_021007.3"}
# coding positions (inclusive) where the alternative transcript is identical to MANE (checked against the sequences below)
CODON = {"TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L", "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L", "ATT": "I", "ATC": "I",
         "ATA": "I", "ATG": "M", "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V", "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
         "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P", "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T", "GCT": "A", "GCC": "A",
         "GCA": "A", "GCG": "A", "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*", "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
         "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K", "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E", "TGT": "C", "TGC": "C",
         "TGA": "*", "TGG": "W", "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R", "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
         "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G"}
A3 = {"Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C", "Gln": "Q", "Glu": "E", "Gly": "G", "His": "H", "Ile": "I",
      "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P", "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V",
      "Ter": "*", "Xaa": "X"}
SUB = re.compile(r"^c\.(\d+)([+-]\d+)?([ACGT])>([ACGT])$")
POS = re.compile(r"^c\.([-*]?\d+)([+-]\d+)?(?:_([-*]?\d+)([+-]\d+)?)?")
PROT = re.compile(r"p\.\(?([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2}|=|Ter)")


def identical_until(a, b):
    """Return the number of leading coding nucleotides that are identical in the two CDS."""
    n = 0
    for x, y in zip(a, b):
        if x != y:
            break
        n += 1
    return n


def identical_spans(a, b):
    """Coding positions (1-based) at which both CDS carry the same base in the same frame (simple prefix/suffix model)."""
    pre = identical_until(a, b)
    suf = identical_until(a[::-1], b[::-1])
    return pre, suf


def norm_c(v):
    v = re.sub(r"\s", "", str(v)).strip()
    v = re.sub(r"^NM_\d+(\.\d+)?:", "", v)
    if re.match(r"^\d+([+-]\d+)?[ACGT]>[ACGT]$", v) or re.match(r"^\d+([+-]\d+)?(_\d+([+-]\d+)?)?(del|dup|ins)", v):
        v = "c." + v
    return v


def in_identical(g, pos, seq, alt, spans):
    """True if coding position pos has the same number and base context on MANE and the alternative transcript."""
    pre, suf = spans[g]
    if pos <= pre:
        return True
    return len(seq[g]) == len(alt[g]) and pos > len(seq[g]) - suf


def protein_fits(cds, pos, ref, alt, prot):
    """True/False whether the substitution at coding position pos encodes the stated protein change; None if not checkable."""
    if pos < 1 or pos > len(cds) or cds[pos - 1] != ref:
        return False
    m = PROT.search(str(prot))
    if not m:
        return None
    i = (pos - 1) // 3
    codon = cds[3 * i:3 * i + 3]
    if len(codon) < 3:
        return False
    mut = list(codon)
    mut[(pos - 1) % 3] = alt
    a0, a1 = CODON.get(codon, "?"), CODON.get("".join(mut), "?")
    want0 = A3.get(m.group(1), "?")
    want1 = a0 if m.group(3) == "=" else A3.get(m.group(3), "?")
    return int(m.group(2)) == i + 1 and want0 == a0 and want1 == a1


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    rows = load_rows()
    seq = {g: fetch_cds(MANE[g]) for g in GENES}
    alt = {g: fetch_cds(t) for g, t in ALT.items()}
    spans = {g: identical_spans(seq[g], alt[g]) for g in ALT}
    rules, report, check = [], [], []
    for g in GENES:
        con = sqlite3.connect(W / f"databases/{g}.sqlite")
        db = {str(m): (vc, vg, tr) for m, vc, vg, tr in con.execute(
            "select mention_id, variant_c_original, variant_c_gnomad, gnomad_transcript from patient_mentions")}
        printed_tx = {str(m): f"{vc} {nm}" for m, vc, nm in con.execute(
            "select mention_id, variant_c_original, nm_transcript from patient_mentions")}
        for x in rows[rows.gene == g].to_dict("records"):
            mid = str(x["mention_id"])
            vc_orig, vg, tr = db.get(mid, ("", "", ""))
            c = norm_c(vg if "c." in str(vg) or re.match(r"^NM_", str(vg)) else x["variant_c"])
            if not c.startswith("c."):
                c2 = norm_c(x["variant_c"])
                c = c2 if c2.startswith("c.") else c
            if not c.startswith("c."):
                continue
            tr = str(tr or x["gnomad_transcript"])
            on_alt = g in ALT and tr.split(".")[0] == ALT[g].split(".")[0]
            on_mane = tr.split(".")[0] == MANE[g].split(".")[0]
            m = SUB.match(c)
            pm = POS.match(c)
            target, why = None, ""
            if m and not m.group(2):
                pos, ref, a = int(m.group(1)), m.group(3), m.group(4)
                fit_mane = protein_fits(seq[g], pos, ref, a, x["protein"])
                fit_alt = protein_fits(alt[g], pos, ref, a, x["protein"]) if g in ALT else False
                ref_mane = 0 < pos <= len(seq[g]) and seq[g][pos - 1] == ref
                ref_alt = g in ALT and 0 < pos <= len(alt[g]) and alt[g][pos - 1] == ref
                if g == "SCN8A" and 615 <= pos <= 702:
                    printed = printed_tx.get(mid, "")
                    if "NM_014191" in printed or (fit_alt and fit_mane is False):
                        target, why = "ALT", "SCN8A exon 5N (NM_014191) as printed in the article or only fitting exon 5N"
                    elif re.search(r"NM_001330260|CCDS81692", printed) or (fit_mane and fit_alt is False):
                        target, why = "MANE", "SCN8A exon 5A (MANE) as printed in the article or only fitting exon 5A"
                    else:
                        target, why = "MANE", "SCN8A exon 5 isoform not stated; kept on MANE numbering (flagged)"
                        check.append(dict(gene=g, canonical_patient_id=x["canonical_patient_id"], mention_id=mid, pmid=x["pmid"],
                                          patient_id=x["patient_id"], variant_as_printed=vc_orig, variant_c=c, transcript=tr,
                                          variant_p=x["protein"], ref_on_mane=ref_mane, problem="SCN8A exon 5A/5N not stated (flag only)"))
                elif pos > max(len(seq[g]), len(alt.get(g, ""))):
                    target, why = "EXCLUDE", f"position beyond the {g} coding sequence on MANE and {ALT.get(g, '-')}: wrong gene"
                elif on_alt and fit_mane is not False and ref_mane and (fit_alt is False or in_identical(g, pos, seq, alt, spans)):
                    target, why = "MANE", "fits MANE (identical region or MANE-only fit)"
                elif on_mane and fit_mane is False and fit_alt and g == "SCN1A" and pos > len(alt[g]) - spans[g][1]:
                    shift = len(seq[g]) - len(alt[g])
                    if protein_fits(seq[g], pos + shift, ref, a, re.sub(r"(\d+)", lambda mm: str(int(mm.group(1)) + shift // 3), str(x["protein"]), count=1)):
                        c = f"c.{pos + shift}{ref}>{a}"
                        target, why = "MANE", f"legacy NM_006920 numbering converted to MANE (+{shift} nt)"
                        rules.append(dict(gene=g, mention_id=mid, pmid=x["pmid"], patient_id=x["patient_id"], field="variant_p",
                                          value=re.sub(r"(\d+)", lambda mm: str(int(mm.group(1)) + shift // 3), str(x["protein"]), count=1),
                                          previous_value=x["protein"], evidence=f"Variant QC: {why}"))
                        on_mane = False  # force relabel below with the new numbering
                elif on_mane and fit_mane is False and fit_alt:
                    target, why = "ALT", f"fits only {ALT[g]} (exon absent from MANE or alternative 3' region)"
                elif on_mane and not ref_mane and ref_alt and pos > len(seq[g]):
                    target, why = "ALT", f"position beyond MANE CDS, valid on {ALT[g]}"
                elif fit_mane is False and not fit_alt:
                    check.append(dict(gene=g, canonical_patient_id=x["canonical_patient_id"], mention_id=mid, pmid=x["pmid"],
                                      patient_id=x["patient_id"], variant_as_printed=vc_orig, variant_c=c, transcript=tr,
                                      variant_p=x["protein"], ref_on_mane=ref_mane,
                                      problem="cDNA and stated protein inconsistent on MANE" + (f" and {ALT[g]}" if g in ALT else "")))
            elif pm and on_alt:
                # intronic / indel: relabel when the anchor positions lie in the region identical to MANE
                p1 = pm.group(1)
                if p1 and not p1.startswith(("-", "*")):
                    p1 = int(p1)
                    p2 = int(pm.group(3)) if pm.group(3) and not pm.group(3).startswith(("-", "*")) else p1
                    if in_identical(g, p1, seq, alt, spans) and in_identical(g, p2, seq, alt, spans):
                        target, why = "MANE", "anchor positions in the coding region identical to MANE"
            vg_new, tr_new = vg, tr
            if target == "MANE":
                vg_new, tr_new = f"{MANE[g]}:{c}", MANE[g]
            elif target == "ALT":
                vg_new, tr_new = f"{ALT[g]}:{c}", ALT[g]
            elif on_mane and not re.match(r"^(NM_\d+\.\d+:)?c\.", str(vg).strip()) and c.startswith("c."):
                vg_new, why = f"{MANE[g]}:{c}", "c. prefix or transcript prefix restored"
            elif str(vg).startswith(MANE[g].split(".")[0]) and not on_mane and target is None and not (m and not m.group(2)):
                # intronic/indel stored on MANE numbering but labelled with another transcript
                tr_new, why = MANE[g], "transcript label aligned with the MANE-numbered variant"
            if target == "EXCLUDE":
                rules.append(dict(gene=g, mention_id=mid, pmid=x["pmid"], patient_id=x["patient_id"], field="mention_source_type",
                                  value="manual_exclude", previous_value="original_patient", evidence=f"Variant QC: {c} {why}"))
            else:
                if vg_new != vg:
                    rules.append(dict(gene=g, mention_id=mid, pmid=x["pmid"], patient_id=x["patient_id"], field="variant_c_gnomad",
                                      value=vg_new, previous_value=vg, evidence=f"Variant QC: {why}"))
                if tr_new != tr:
                    rules.append(dict(gene=g, mention_id=mid, pmid=x["pmid"], patient_id=x["patient_id"], field="gnomad_transcript",
                                      value=tr_new, previous_value=tr, evidence=f"Variant QC: {why}"))
            report.append(dict(gene=g, mention_id=mid, variant=c, transcript=tr, decision=target or "", why=why))
    r = pd.DataFrame(rules)
    r["source_run"] = "runs/review_fixes/variant_qc"
    r["status"] = "active"
    r.to_csv(OUT / "rules.csv", index=False)
    pd.DataFrame(report).to_csv(OUT / "report.csv", index=False)
    pd.DataFrame(check).to_csv(OUT / "source_check_needed.csv", index=False)
    print("identical prefix/suffix:", {g: spans[g] for g in spans}, "CDS lengths:", {g: (len(seq[g]), len(alt[g])) for g in ALT})
    if len(r):
        print(r.groupby(["gene", "field"]).size().to_string())
        print(r.evidence.value_counts().head(10).to_string())
    print("source check needed:", len(check))


if __name__ == "__main__":
    main()
