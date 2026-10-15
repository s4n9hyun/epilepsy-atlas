"""Rule-based coding of reported treatment response per patient and treatment class.

Evidence text comes from targeted_treatment_response and response_to_asm (narrative), responsive_asm and asm_used.
Text is split into segments (; newline sentence) and then clauses (, and with while). A statement is coded for a
class only when a specific drug of the class is named in the same clause as the outcome words, or in a segment that
names no drug of another class. Extraction-template prefixes ("Sodium channel blocker benefit:") are removed first so
that the generic class name never counts as evidence.

Per (patient, class): Worsening if any statement reports worsening, else Benefit, else No benefit, else
Used, response not reported. v2 rules after validation round 1 (see analysis/treatment_validation_reviewed.csv).
"""
import re

import pandas as pd

from common import is_missing

CLASSES = {
    "Sodium channel blocker": r"carbamazepine|\bCBZ\b|oxcarbazepine|\bOXC\b|phenytoin|fosphenytoin|\bPHT\b|\bDPH\b|lacosamide|"
                              r"\bLCM\b|lamotrigine|\bLTG\b|eslicarbazepine|\bESL\b|mexiletine|lidocaine|"
                              r"(?<!\w)(?:sodium|Na\+?) ?channel[- ]block(?:er|ers|ing|ade)?\b(?! benefit| worsening| data| treatment response)",
    "Ketogenic diet": r"ketogenic|\bKD\b|modified atkins|\bMAD\b|low glycemic",
    "Cannabidiol": r"cannabidiol|\bCBD\b|epidiolex",
    "Fenfluramine": r"fenfluramine|\bFFA\b|\bFFL\b",
    "Stiripentol": r"stiripentol|\bSTP\b",
    "Retigabine": r"retigabine|ezogabine",
    "Vigabatrin/ACTH/steroids": r"vigabatrin|\bVGB\b|\bACTH\b|prednis|corticotropin|steroid",
    "Quinidine": r"quinidine|\bQND\b",
}
OTHER_ASM = r"levetiracetam|\bLEV\b|valpro|\bVPA\b|topiramate|\bTPM\b|clobazam|\bCLB\b|clonazepam|\bCZP\b|phenobarb|\bPB\b|\bPHB\b|" \
            r"zonisamide|\bZNS\b|perampanel|\bPER\b|brivaracetam|rufinamide|ethosuximide|sulthiame|\bSTM\b|pyridox|\bPLP\b|\bB6\b|" \
            r"nitrazepam|diazepam|midazolam|lorazepam|bromide|acetazolamide|memantine|gabapentin|pregabalin|felbamate|" \
            r"cenobamate|everolimus|IVIG|immunoglobulin"
TEMPLATE = re.compile(r"^\s*(?:(?:no |no clear )?(?:sodium channel blocker|SCB)s?|benefit from [\w\s/-]+?|[\w\s/-]+? (?:benefit|worsening))"
                      r"\s*(?:benefit|worsening|benefit or worsening)?(?:/\w+)?\s*:\s*", re.I)
GENERAL = r"may (aggravat|worsen|exacerbat)|might (aggravat|worsen|exacerbat)|can (aggravat|worsen)|known to|contraindicat|" \
          r"in general|generally|in (dravet|the) (syndrome|literature)|cohort|aggregate|all such patients|group|some patients|" \
          r"not linked|\b\d+ (of|/) ?\d+\b|patients with \w+ variants"
NOT_GIVEN = r"not used|were not (given|tried|used)|avoid|no (sodium channel blocker|SCB) data|never (received|given|treated)"
UNSPEC = r"not specified|not explicitly|no explicit|not stated|not reported|not described|unclear|unknown|timeline|" \
         r"not isolated|cannot be (assessed|determined)|not clearly|cohort[- ]level|not linked|dose (was )?reduced|used before|" \
         r"before (accurate |genetic )?diagnosis|no \w*\s?(benefit|response|worsening|effect) (was )?(reported|specified|described)|" \
         r"not (evaluable|assessed)|individual response|" \
         r"no (clear )?(\w+ )?(benefit|response|worsening|effect)( from [\w\s,/+-]+?)? (was )?(reported|specified|described)"
IMPROVE_OFF = r"(improv|seizure[- ]free|reduc|decreas|resolv|ceased|stopped)\w*\s.{0,40}(after|following|upon|on) " \
              r"(stopping|withdraw|discontinu|taper|cessation|removal|switch)|(withdraw|discontinu|taper|stopp)\w*\s.{0,70}(improv|seizure[- ]free|reduc|decreas)"
WORSE_OFF = r"(worse|recur|increas|relaps|exacerbat)\w*\s.{0,40}(after|following|upon|on|when|during) ([\w-]+ ){0,3}(was |were )?(reduc|withdraw|discontinu|taper|wean|stopp)|" \
            r"(reduc|withdraw|discontinu|taper|wean)\w*\s.{0,40}(associated with|led to|followed by|caused) (increas|worse|recur)"
WORSE = r"worsen|worse\b|aggravat|exacerbat|deteriorat|increase[sd]? (in )?(seizure|frequency)|increased seizures|provok|paradoxical|" \
        r"induced (seizure|status)|new seizure type"
NOT_WORSE = r"no (clear |obvious |significant )?(seizure )?(worsen|aggravat|exacerbat)|without (worsen|aggravat|exacerbat)|" \
            r"not (worsen|aggravat)|not found to (exacerbat|worsen|aggravat)|attenuat|kainate|KA-induced|" \
            r"worsening (of|in) (irritab|sleep|behav|develop)|(behavio|sleep)[\w\s]{0,30}deteriorat"
NO_BENEFIT = r"no (clear |significant |obvious |sustained |)(benefit|effect|response|improvement|change)|ineffective|not effective|" \
             r"without (significant |clear |any )?(benefit|effect|improvement)|not found to be effective|persist|N/E\b|" \
             r"no seizure freedom|unresponsive|refractory to|failed|no seizure reduction|did not (respond|improve|help)|" \
             r"not (helpful|successful|beneficial|valid|observed|preventing|useful)|resistant to|<\s?50%|less than 50%|non-?responder|" \
             r"lack of (efficacy|response)|response 0%|0% (response|reduction)|good efficacy\s*=\s*no|no effect"
BENEFIT = r"benefit|effective|efficac|seizure[- ]free|controlled|control of|responded|respond(ed|er|s)? (well|to)|improv|reduc|" \
          r"cessation|stopp?ed (the )?seizures|≥\s?50%|>\s?50%|50[-–]99%|100%|decrease|remission|dramatic|good response|resolved|" \
          r"helpful|partial response|responsive"
DESPITE = r"despite|although on|while on|in spite of"


def _segments(text):
    return [s.strip() for s in re.split(r"[;\n]|(?<=[a-z\)])\. ", str(text)) if s.strip()]


def _clauses(seg):
    return [c.strip() for c in re.split(r",\s|\s(?:while|whereas|but)\s", seg) if c.strip()]


def _code(text):
    """Outcome of one statement already attributed to the class."""
    t = text
    if re.search(GENERAL, t, re.I):
        return "Used, response not reported"
    if re.search(IMPROVE_OFF, t, re.I):
        return "Worsening"  # seizures improved once the drug was stopped
    if re.search(WORSE_OFF, t, re.I):
        return "Benefit"  # seizures worsened when the drug was reduced or stopped
    if re.search(UNSPEC, t, re.I) and not re.search(r"seizure[- ]free|\d+\s?%|improv|reduc", t, re.I):
        return "Used, response not reported"
    if re.search(WORSE, t, re.I) and not re.search(NOT_WORSE, t, re.I):
        return "Worsening"
    if re.search(NO_BENEFIT, t, re.I):
        return "No benefit"
    if re.search(DESPITE, t, re.I):
        return "Used, response not reported"
    if re.search(BENEFIT, t, re.I):
        return "Benefit"
    return "Used, response not reported"


def classify(row, pattern):
    """Return (outcome, evidence) for one treatment class, or (None, None) if the class is never named.

    evidence: 'narrative' (coded from a free-text statement), 'responsive_field' (drug listed as responsive ASM only),
    'used_only' (named, no response).
    """
    drug = re.compile(pattern, re.I)
    others = re.compile("|".join([OTHER_ASM] + [p for p in CLASSES.values() if p != pattern]), re.I)
    order = {"Worsening": 3, "Benefit": 2, "No benefit": 1, "Used, response not reported": 0}
    found = []
    for field in ("targeted_treatment_response", "response_to_asm"):
        for seg in _segments(row.get(field, "")):
            seg = TEMPLATE.sub("", seg)
            if not drug.search(seg) or re.search(NOT_GIVEN, seg, re.I):
                continue
            if re.search(GENERAL, seg, re.I) or re.search(IMPROVE_OFF, seg, re.I) or re.search(WORSE_OFF, seg, re.I):
                found.append((_code(seg), "narrative" if not re.search(GENERAL, seg, re.I) else "used_only"))
                continue
            outcome_words = f"{WORSE}|{NO_BENEFIT}|{BENEFIT}|{UNSPEC}"
            hits = [c for c in _clauses(seg) if drug.search(c) and re.search(outcome_words, c, re.I)]
            if hits:
                found += [(_code(c), "narrative") for c in hits]
            elif not others.search(drug.sub("", seg)):
                found.append((_code(seg), "narrative"))
            else:
                found.append(("Used, response not reported", "used_only"))
    ra = str(row.get("responsive_asm", ""))
    if not is_missing(ra) and drug.search(ra):
        found.append(("Benefit", "responsive_field"))
    if not found and drug.search(str(row.get("asm_used", ""))):
        found.append(("Used, response not reported", "used_only"))
    if not found:
        return None, None
    best = max(found, key=lambda x: (order[x[0]], x[1] == "narrative"))
    return best[0], ("used_only" if best[0] == "Used, response not reported" else best[1])


def code_treatments(patients):
    out = []
    for _, r in patients.iterrows():
        rec = r.to_dict()
        for cls, pat in CLASSES.items():
            o, ev = classify(rec, pat)
            if o:
                out.append({"gene": r.gene, "canonical_patient_id": r.canonical_patient_id, "pmid": r.pmid,
                            "treatment": cls, "outcome": o, "evidence": ev})
    return pd.DataFrame(out)
