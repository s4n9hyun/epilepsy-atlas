"""Shared loaders and harmonization for the DEE atlas manuscript.

Reads the frozen canonical patient-variant CSVs from the curation workspace (read-only).
"""
import glob
import re
import sqlite3
from pathlib import Path

import pandas as pd

WORKSPACE = Path(__import__("os").environ.get("EPILEPSY_WORKSPACE", Path(__file__).resolve().parents[2]))
PROJECT = Path(__file__).resolve().parents[1]
FIG = PROJECT / "figures"
TAB = PROJECT / "tables"
GENES = ["SCN1A", "SCN2A", "SCN8A", "KCNQ2", "KCNT1", "CDKL5", "PCDH19", "STXBP1", "SYNGAP1", "GABRA1", "FOXG1", "SMC1A", "EEF1A2"]
GENES6 = ["SCN1A", "SCN2A", "KCNQ2", "CDKL5", "SYNGAP1", "EEF1A2"]  # first release (curated with clinician CRF review)
NS = {"", "not specified", "unknown", "not available", "not applicable", "nan", "none", "na", "n/a", "nr"}


def is_missing(v):
    return str(v).strip().lower() in NS


def load_rows(genes=GENES):
    frames = []
    for g in genes:
        f = sorted(glob.glob(str(WORKSPACE / f"outputs/current/canonical_patients/{g}_canonical_patients/*.csv")))[-1]
        df = pd.read_csv(f, dtype=str, keep_default_na=False, encoding="utf-8-sig")
        df["gene"] = g
        frames.append(df)
    df = pd.concat(frames, ignore_index=True)
    df["sex_h"] = df["sex"].map(lambda s: {"male": "Male", "m": "Male", "female": "Female", "f": "Female"}.get(s.strip().lower(), "NS"))
    df["onset_months"] = df["age_of_onset"].map(onset_months)
    df["phenotype_h"] = [harmonize_phenotype(g, p, s) for g, p, s in zip(df.gene, df.phenotype_review_category, df.seizure_presence)]
    df["protein"] = [fill_protein(g, c, p) for g, c, p in zip(df.gene, df.variant_c, df.variant_p)]
    df["consequence"] = [consequence(c, p, t) for c, p, t in zip(df.variant_c, df.protein, df.variant_type)]
    df["cnv_type"] = [cnv_type(c, k) for c, k in zip(df.variant_c, df.consequence)]
    df["aa_pos"] = df["protein"].map(aa_position)
    return df


def load_patients(rows=None):
    """One row per canonical patient; first variant row carries patient-level fields."""
    rows = load_rows() if rows is None else rows
    first = rows.drop_duplicates(["gene", "canonical_patient_id"]).copy()
    nvar = rows.groupby(["gene", "canonical_patient_id"]).size().rename("n_variant_rows")
    return first.merge(nvar, on=["gene", "canonical_patient_id"])


def load_mentions(genes=GENES):
    frames = []
    for g in genes:
        con = sqlite3.connect(WORKSPACE / "databases" / f"{g}.sqlite")
        m = pd.read_sql("select mention_id, pmid, patient_id, mention_source_type, cited_source, sex, "
                        "phenotype_review_category, seizure_presence, age_of_onset, variant_c_original, variant_c_gnomad, variant_p "
                        "from patient_mentions", con, dtype=str)
        m["gene"] = g
        frames.append(m)
    return pd.concat(frames, ignore_index=True).fillna("")


# ---------------------------------------------------------------- source-verified treatment responses
VERIFIED = [("scb_source_verified.csv", "Sodium channel blocker"),          # SCN1A, SCN2A (runs/paper/source_qc)
            ("treatment_verified_7genes.csv", None)]                        # SCN8A SCB, KCNT1 quinidine (runs/paper/treatment_verify_7genes)


def verified_labels():
    frames = []
    for f, tr in VERIFIED:
        path = PROJECT / "analysis" / f
        if path.exists():
            d = pd.read_csv(path, dtype=str, keep_default_na=False)
            if tr:
                d["treatment"] = tr
            frames.append(d[["gene", "canonical_patient_id", "pmid", "treatment", "source_label", "db_text_correct"]])
    return pd.concat(frames, ignore_index=True)


# ---------------------------------------------------------------- age of onset
UNIT = {"day": 1 / 30.44, "d": 1 / 30.44, "week": 7 / 30.44, "wk": 7 / 30.44, "w": 7 / 30.44,
        "month": 1.0, "mo": 1.0, "m": 1.0, "year": 12.0, "yr": 12.0, "y": 12.0, "hour": 1 / 730.5, "h": 1 / 730.5}
WORDS = [(r"\b(neonat|newborn|first (few )?(days|week)|1st week|birth|day of life|days? of (postnatal )?life|perinatal|DOL|few days|several days)", 0.1),
         (r"first month", 0.5),
         (r"\bin utero|prenatal|fetal\b", 0.0)]


def onset_months(text):
    t = str(text).strip().lower()
    if is_missing(t):
        return None
    t = t.replace(",", ".").replace("–", "-").replace("~", "-")
    m = re.search(r"(\d+(?:\.\d+)?)\s*(?:-|to)\s*(\d+(?:\.\d+)?)\s*(hours?|h|days?|d|weeks?|wks?|w|months?|mos?|m|years?|yrs?|y)\b", t)
    if m:
        u = re.sub(r"s$", "", m.group(3))
        return (float(m.group(1)) + float(m.group(2))) / 2 * UNIT.get(u, UNIT.get(u[0], 1))
    total, hit = 0.0, False
    for num, u in re.findall(r"(\d+(?:\.\d+)?)\s*(hours?|h|days?|d|weeks?|wks?|w|months?|mos?|m|years?|yrs?|y)\b", t):
        u = re.sub(r"s$", "", u)
        total += float(num) * UNIT.get(u, UNIT.get(u[0], 1))
        hit = True
    if hit:
        return total
    for pat, val in WORDS:
        if re.search(pat, t, re.I):
            return val
    m = re.fullmatch(r"(\d+(?:\.\d+)?)", t)
    return None if not m else None  # bare numbers have no unit: treat as missing


# ---------------------------------------------------------------- phenotype
DEE = "DEE"
PHENO_ORDER = ["DEE", "Self-limited (familial) neonatal/infantile epilepsy", "GEFS+ / febrile spectrum",
               "Generalized epilepsy", "Other epilepsy", "No seizures", "Not specified"]


def harmonize_phenotype(gene, cat, seizure):
    c = str(cat).lower()
    if is_missing(c):
        return "No seizures" if str(seizure).strip().lower() == "no" else "Not specified"
    if "non-seizure" in c or "never reported" in c or "no epilepsy" in c:
        return "No seizures"
    if any(k in c for k in ["dravet", "encephalopathy", "ohtahara", "west", "spasms", "lennox", "migrating focal"]):
        return DEE
    if gene == "CDKL5" and c == "epilepsy":
        return DEE  # CDKL5 deficiency disorder is a DEE by definition
    if "benign" in c or "self-limited" in c:
        return "Self-limited (familial) neonatal/infantile epilepsy"
    if "gefs" in c or "febrile" in c:
        return "GEFS+ / febrile spectrum"
    if gene == "SYNGAP1" and re.search(r"myoclonic[- ]a(tonic|static)|doose", c):
        return DEE  # epilepsy with myoclonic-atonic seizures is a DEE in ILAE 2022 (Specchio 2022)
    if gene == "SYNGAP1" and any(k in c for k in ["generalized", "myoclonic", "absence"]):
        return "Generalized epilepsy"
    return "Other epilepsy"


# ---------------------------------------------------------------- variant consequence
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts" / "core"))
from audit_variant_consistency import check_substitution, fetch_cds  # noqa: E402

MANE = {"SCN1A": "NM_001165963.4", "SCN2A": "NM_001040142.2", "CDKL5": "NM_001323289.2", "EEF1A2": "NM_001958.5",
        "KCNQ2": "NM_172107.4", "SYNGAP1": "NM_006772.3", "SCN8A": "NM_001330260.2", "KCNT1": "NM_020822.3",
        "PCDH19": "NM_001184880.2", "STXBP1": "NM_001032221.6", "GABRA1": "NM_001127644.2", "FOXG1": "NM_005249.5",
        "SMC1A": "NM_006306.4"}
_CDS = {}
CNV_RE = re.compile(r"\b(Mb|kb)\b|exon|whole[- ]gene|gene deletion|microdeletion|microduplication|\barr\b|"
                    r"del(etion)?\s+(of\s+)?(SCN|KCNQ|KCNT|CDKL|SYNGAP|EEF|STXBP|PCDH|FOXG|SMC|GABR)|\b\d+q\d|\bXp\d|\bXq\d|^del|^dup|^deletion|duplication|外显子", re.I)
P_MISSENSE = re.compile(r"p\.\(?[A-Z][a-z]{2}\d+(?!Ter)[A-Z][a-z]{2}\)?$")


def fill_protein(gene, c, p):
    """Translate MANE cDNA substitutions when the protein change was not reported."""
    if not is_missing(p) or gene not in MANE:
        return p
    m = re.fullmatch(r"c\.(\d+)([ACGT])>([ACGT])", str(c).strip())
    if not m:
        return p
    cds = _CDS.setdefault(gene, fetch_cds(MANE[gene]))
    try:
        ref_ok, _, expected = check_substitution(cds, int(m.group(1)), m.group(2), m.group(3), "")
    except Exception:
        return p
    return expected if ref_ok and expected else p


def consequence(c, p, vtype):
    """Predicted consequence from the MANE cDNA first, then the protein change (reviewer comments R3-5, R4 minor 18).
    Truncating = nonsense, frameshift, canonical +/-1/2 splice or start loss; in-frame delins/del/dup = in-frame indel;
    intronic beyond +/-2 and synonymous = other splice/intronic/synonymous."""
    c, p, vt = str(c).strip(), str(p).strip(), str(vtype).lower()
    if "cnv" in vt or (CNV_RE.search(c) and not c.startswith("c.")):
        return "CNV / exon-level"
    if re.search(r"\d{2,3},?\d{3},?\d{3}\s*[-–_]", c) or re.search(r"large deletion|heterozygous deletion|complete .*deletion", c, re.I):
        return "CNV / exon-level"
    if c.startswith("c."):
        if re.search(r"c\.-?\d+[+-][12](?!\d)[ACGT_d]|c\.-?\d+[+-][12]$|_-?\d+[+-][12](?!\d)", c):
            return "Truncating"  # canonical splice site
        if re.search(r"c\.-?\d+[+-]\d+|c\.\*|c\.-\d", c):
            return "Other splice/intronic/synonymous"  # non-canonical intronic or UTR
        if re.match(r"c\.[123][ACGT]>[ACGT]$", c) and re.match(r"c\.(1A>[CGT]|2T>[ACG]|3G>[ACT])$", c):
            return "Truncating"  # start loss
        fs = _frameshift_c(c)
        if fs is True:
            return "Truncating"
        if fs is False and re.search(r"del|dup|ins", c):
            return "In-frame indel"
    if re.search(r"=\)?$", p):
        return "Other splice/intronic/synonymous"
    if re.search(r"fs|trunc", p, re.I) or re.search(r"(Ter|\*|X)\)?$", p) or re.search(r"fsX|X\d+$|^Truncat", c):
        return "Truncating"
    if re.search(r"^[A-Z]?\d+[+-][12][ACGT]>|(?:IVS)-?\d+[+-][12](?!\d)", c):
        return "Truncating"
    if P_MISSENSE.search(p):
        return "Missense"
    if re.search(r"p\.\(?[A-Z][a-z]{2}\d+(_[A-Z][a-z]{2}\d+)?(del|ins|dup|delins)", p):
        return "In-frame indel"
    if re.search(r"(?:c\.|IVS)-?\d+[+-]\d+", c):
        return "Other splice/intronic/synonymous"
    return "Unclassified"


def cnv_type(c, cons):
    """gain / loss / unknown for copy-number rows (FOXG1 14q12 gains and losses have different phenotypes)."""
    if cons != "CNV / exon-level":
        return ""
    t = str(c).lower()
    if re.search(r"dup|gain|triplicat|amplif|trisomy|x4|x3\b|supernumerary", t):
        return "gain"
    if re.search(r"del|loss|monosomy|x1\b", t):
        return "loss"
    return "unknown"


def _frameshift_c(c):
    """True/False for frameshift vs in-frame indels on the cDNA; None if not an indel or length unknown."""
    m = re.search(r"c\.(\d+)(?:_(\d+))?delins([ACGT]+)$", c)
    if m:
        dl = (int(m.group(2)) - int(m.group(1)) + 1) if m.group(2) else 1
        return (len(m.group(3)) - dl) % 3 != 0
    m = re.search(r"c\.(\d+)(?:_(\d+))?(del|dup)([ACGT]*)$", c)
    if m:
        n = len(m.group(4)) if m.group(4) else (int(m.group(2)) - int(m.group(1)) + 1 if m.group(2) else 1)
        return n % 3 != 0
    m = re.search(r"ins([ACGT]+)$", c)
    if m:
        return len(m.group(1)) % 3 != 0
    return None


# ---------------------------------------------------------------- family grouping (reviewer comments R2-M2, R3-2/6)
FAM_RE = re.compile(r"\b(?:family|fam\.?|kindred|pedigree|ped\.?)\s*#?\s*([A-Za-z0-9]+(?:[-_][A-Za-z0-9]+)*)", re.I)
FAM_CODE = re.compile(r"^([A-Z]{0,3}\d*[A-Za-z]?)\s*[:\s]\s*(?=[IVX]{1,4}\s*[-:.]\s*\d)")
REL_RE = re.compile(r"\b(?:mother|father|parent|sister|brother|sibling|sib|son|daughter|twin|cousin|aunt|uncle|grand\w*|half-\w+)\b"
                    r"\s+of\s+(?:the\s+)?(?:patient|proband|case|index|subject|individual|pt\.?|p)?\s*#?\s*([A-Za-z]?\d+[A-Za-z]?)", re.I)
GEN_RE = re.compile(r"\b(?:[IVX]{1,4})\s*[-:.]\s*\d+\b")
PROBAND_RE = re.compile(r"proband|index|propositus|\bpatient\b(?!.*\bof\b)", re.I)


def family_ids(df):
    """Heuristic family identifier within each article: explicit family/kindred/pedigree tokens, 'mother/sister of patient N'
    labels, and pedigree generation labels (e.g. II-3) of inherited carriers of the same variant. Relatives first reported in
    different articles are not linked. Returns (family_id, is_proband) aligned with df."""
    fam, prob = [], []
    for g, pmid, lab, c, inh in zip(df.gene, df.pmid, df.patient_id.astype(str), df.variant_c.astype(str),
                                    df.get("inheritance_h", pd.Series([""] * len(df))).astype(str)):
        m = FAM_RE.search(lab) or FAM_CODE.search(lab.strip())
        r = REL_RE.search(lab)
        if m:
            fid = f"{g}|{pmid}|F{m.group(1).upper()}"
        elif r:
            fid = f"{g}|{pmid}|P{r.group(1).upper()}"
        elif GEN_RE.search(lab):
            fid = f"{g}|{pmid}|V{re.sub(r'\s', '', c)}"
        else:
            fid = f"{g}|{pmid}|{lab}"
        fam.append(fid)
    df = df.assign(_fid=fam)
    # proband: an explicitly labelled proband/index, else the first row of the family
    order = df.assign(_p=[0 if re.search(r"proband|index|propositus", l, re.I) else 1 for l in df.patient_id.astype(str)])
    first = order.sort_values(["_fid", "_p", "mention_id"]).drop_duplicates("_fid").index
    prob = df.index.isin(first)
    return pd.Series(fam, index=df.index), pd.Series(prob, index=df.index)


def aa_position(p):
    m = re.search(r"p\.\(?[A-Z][a-z]{2}(\d+)", str(p))
    return int(m.group(1)) if m else None
