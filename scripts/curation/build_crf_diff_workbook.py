"""0904 통합본(교수님 CRF 통합본)을 그대로 복사해 현재 DB와의 차이를 git diff처럼 표시한 대조본을 만든다.

- 시트·열·헤더·드롭다운은 0904 통합본과 같고, 맨 오른쪽에 '비고' 열 하나만 추가한다.
  · 변경 없음: 통합본 행 그대로
  · 수정: 한 행(연노랑). 바뀐 칸만 진노랑+굵게 현재 값, 0904 값은 칸 메모와 비고
  · 추가(현재 DB에만): 연초록 행 — 같은 PMID 행들 뒤, 없으면 시트 끝
  · 삭제(현재 DB에 없음): 연빨강+취소선 행
  · 비고: 임상의가 읽는 문장으로, 무엇을 왜 바꿨는지와 원문 근거
- 정규화 후 같은 값(대소문자, N/A·Not specified, 단위 표기, 라벨 표기, chr 접두어)은 차이로 보지 않고 통합본 값을 유지.
- 행 짝짓기: scripts/curation/compare_crf_vs_canonical.py, 현재 행 변환: scripts/curation/export_crf_workbook.py,
  사유: outputs/current/curation_change_log.csv, 없으면 DB 언급 분류·제외 사유, 통합본 바탕 DB(6/16)와 비교.
"""
import glob
import importlib.util
import re
import shutil
from collections import defaultdict
from pathlib import Path

import pandas as pd
from openpyxl import Workbook
from openpyxl.styles import Alignment, Border, Font, PatternFill, Side
from openpyxl.utils import get_column_letter
from openpyxl.worksheet.datavalidation import DataValidation

ROOT = Path(__file__).resolve().parents[2]
DECISIONS = ROOT / "config/crf_0904"  # 0904 통합본 대조 판정 (임상의 입력값 미채택·충돌, 짝 없는 행 판정)
spec = importlib.util.spec_from_file_location("cmp", ROOT / "scripts/curation/compare_crf_vs_canonical.py")
cmp = importlib.util.module_from_spec(spec)
spec.loader.exec_module(cmp)

GENES = cmp.GENES
CORE = ["variant_c", "variant_p", "sex", "age_of_onset", "last_followup_age", "phenotype"]
SHOW = [("variant_c", "cDNA 변이"), ("variant_p", "단백질 변이"), ("sex", "성별"), ("age_of_onset", "발병 나이"),
        ("last_followup_age", "마지막 추적 나이"), ("phenotype", "표현형"), ("targeted_response", "표적 치료 반응"),
        ("clinvar_pathogenicity", "ClinVar 분류"), ("clinvar_id", "ClinVar ID")]
FIELD_KO = dict(SHOW)
LOG_FIELD = {"variant_c": "Variant cDNA change", "variant_p": "Variant protein change", "sex": "Sex",
             "age_of_onset": "Age at seizure onset", "last_followup_age": "Age at last follow-up", "phenotype": "Phenotype",
             "targeted_response": "Response to targeted ASM", "clinvar_pathogenicity": "ClinVar pathogenicity",
             "clinvar_id": "ClinVar ID"}
KIND_KO = {"modified": "수정", "added": "추가 (현재 DB에만)", "crf_only": "CRF에만 있음"}

RED, GREEN, GREY = PatternFill("solid", fgColor="F8D7DA"), PatternFill("solid", fgColor="D4EDDA"), PatternFill("solid", fgColor="F2F2F2")
RED_CELL, GREEN_CELL = PatternFill("solid", fgColor="F1A7AE"), PatternFill("solid", fgColor="8FD19E")
HEAD, PAPER = PatternFill("solid", fgColor="1F4E79"), PatternFill("solid", fgColor="D9E1F2")
THIN = Border(bottom=Side(style="thin", color="BFBFBF"))
_GS = Side(style="thin", color="FFBFBFBF")
GRID = Border(left=_GS, right=_GS, top=_GS, bottom=_GS)

# ---------------------------------------------------------------- reference data
LOG = pd.read_csv(ROOT / "outputs/current/curation_change_log.csv", dtype=str, keep_default_na=False, encoding="utf-8-sig")
LABEL_MAP = pd.read_csv(ROOT / "config/crf_phenotype_label_map.csv", dtype=str, keep_default_na=False)
LABEL = {(g, c.strip().lower()): p.strip().lower() for g, c, p in zip(LABEL_MAP.gene, LABEL_MAP.canonical_value, LABEL_MAP.crf_phenotype)}
TITLES = {}
for g in GENES:
    import sqlite3
    con = sqlite3.connect(ROOT / "databases" / f"{g}.sqlite")
    TITLES.update(dict(con.execute("select pmid, title from pubmed_articles").fetchall()))
MENTIONS = {}
for g in GENES:
    frames = [pd.read_csv(sorted(glob.glob(str(ROOT / f"outputs/auxiliary/{g}/{g}_{k}_*.csv")))[-1], dtype=str,
                          keep_default_na=False, encoding="utf-8-sig") for k in ("patient_mentions", "excluded_mentions")]
    MENTIONS[g] = pd.concat(frames).drop_duplicates("mention_id")
SOURCE_TYPE_KO = {"cited_patient": "이전 논문에 보고된 환자를 다시 인용한 행이라, 원래 논문의 환자로만 남겼습니다.",
                  "previously_reported_patient": "이미 발표된 환자를 다시 보고한 행이라, 원래 논문의 환자로만 남겼습니다.",
                  "reanalyzed_previous_patient": "이전에 보고된 환자를 재분석한 행이라, 원래 논문의 환자로만 남겼습니다.",
                  "aggregate_only": "코호트 요약만 있고 개별 환자 정보가 없어 제외했습니다.",
                  "variant_only": "변이만 있고 환자별 임상 정보가 없어 제외했습니다.",
                  "unclear": "원래 보고된 환자인지 원문으로 확인되지 않아 제외했습니다.",
                  "manual_exclude": "원문을 다시 확인해 제외했습니다."}
EXCLUSION_KO = {"no_patient_level_clinical_field": "환자별 임상 정보가 없어 제외했습니다.",
                "variant_not_recorded": "원문에 이 환자의 변이가 기재되어 있지 않아 제외했습니다.",
                "benign_or_polymorphism_in_paper": "원문에서 양성 변이 또는 다형성으로 분류되어 제외했습니다.",
                "cited_from_prior_literature": SOURCE_TYPE_KO["cited_patient"],
                "aggregate_only_no_individual_patient": SOURCE_TYPE_KO["aggregate_only"],
                "variant_only_no_patient_level_clinical_information": SOURCE_TYPE_KO["variant_only"],
                "previously_reported_or_reanalyzed_patient": SOURCE_TYPE_KO["previously_reported_patient"],
                "manual_excluded_after_source_review": SOURCE_TYPE_KO["manual_exclude"],
                "unclear_originality": SOURCE_TYPE_KO["unclear"],
                "non_causal_secondary_variant_manual_exclude": "원인 변이가 아닌 부수 변이라 제외했습니다."}
REASON_KO = [("variant_prefix_fix", "이전 파일에서 잘려 있던 변이 표기(염색체 좌표 부분)를 원문 표기대로 복원했습니다."),
             ("crf_clinician_values", "0904 통합본에 입력하신 값을 반영했습니다."),
             ("clinvar", "ClinVar를 다시 조회해 최신 정보로 바꿨습니다."),
             ("reextract", "보충자료(supplementary)를 다시 추출하면서 원문을 확인해 바로잡았습니다."),
             ("professor_feedback", "교수님 검토 의견을 원문에서 확인해 반영했습니다."),
             ("cohort_value", "코호트 전체 요약값이 개별 환자 값으로 잘못 들어가 있어, 원문을 확인해 환자별 값으로 바로잡았습니다."),
             ("duplicate", "같은 환자가 여러 논문에 중복 보고되어, 먼저 보고한 논문 쪽으로 합쳤습니다."),
             ("rule8", "처음으로 유전형을 보고한 논문을 원래 보고로 정리했습니다."),
             ("_add", "검색식(유전자 AND epilepsy)에 걸리지 않던 논문을 추가로 검토해 반영했습니다."),
             ("regroup", "같은 변이의 행을 원문 기준으로 다시 정리했습니다."),
             ("_exclude", "원문을 다시 확인해 제외했습니다."),
             ("crf_adopted_source_check", "0904 통합본 입력값을 원문과 다시 대조해 원문 기준으로 바로잡았습니다."),
             ("full_source_check", "원문을 전수 대조해 바로잡았습니다."),
             ("", "원문을 다시 확인해 바로잡았습니다.")]  # paper_qc, data_cleanup, crf_review, variant_review 등
DELETE_KO = [("full_source_check_followup", "원문을 전수 확인한 결과 증상이 없는 보인자 가족이나 유전자 검사를 받지 않은 가족이거나, 먼저 보고된 논문의 같은 환자(재보고), 또는 개별 환자가 아닌 행이라 뺐습니다."),
             ("full_source_check", "원문을 전수 대조한 결과 먼저 보고된 논문의 같은 환자(재보고)이거나, 개별 환자가 아닌 행(가족·코호트 요약, 세포·동물 연구, 논문에 없는 행)이라 뺐습니다."),
             ("cell_line_rows_check", "원문을 확인한 결과 먼저 보고된 논문의 같은 환자이거나(중복), 이 논문에는 환자 개인의 임상 정보가 없어(환자 유래 세포주 연구) 뺐습니다."),
             ("crf_adopted_source_check", "원문을 확인한 결과 사람 환자가 아니라 동물 모델(쥐·생쥐)이어서 제외했습니다."),
             ("duplicate", "같은 환자가 여러 논문에 중복 보고되어, 먼저 보고한 논문의 환자로 합쳤습니다."),
             ("regroup", "같은 환자의 같은 변이 행이 중복되어 하나로 합쳤습니다."),
             ("reextract", "보충자료(supplementary)를 다시 추출하면서 원문을 확인해 제외했습니다."),
             ("professor_feedback", "교수님 검토 의견을 원문에서 확인해 제외했습니다."),
             ("", "원문을 다시 확인해 제외했습니다.")]
ADD_KO = [("first_report_fix", "논문 원문을 전수 확인하면서 찾은, 이 논문에서 처음 보고된 환자입니다(이전 자료에서 빠져 있었습니다)."),
          ("full_source_check", "논문 원문을 전수 확인하면서 찾은, 이 논문에서 처음 보고된 환자입니다(이전 자료에서 빠져 있었습니다)."),
          ("reextract", "보충자료(supplementary)를 다시 추출해 찾은 환자입니다."),
          ("crf_added_papers", "0904 통합본에 임상의가 추가하신 환자로, 원문을 확보해 확인했습니다."),
          ("_add", "여러 유전자를 함께 다룬 코호트 논문으로, 다른 유전자 목록에는 이미 있었으나 이 유전자의 검색식(유전자 AND epilepsy)에는 "
                   "걸리지 않았습니다. 보충자료 표에 이 유전자 환자가 환자별로 실려 있어 추가했습니다."),
          ("rule8", "처음으로 유전형을 보고한 논문을 원래 보고로 정리하면서 추가된 환자입니다."),
          ("", "원문을 다시 확인해 찾은 환자입니다.")]


def batch_ko(batch, table=None):
    for k, v in table or REASON_KO:
        if k in batch:
            return v
    return batch


def log_for(gene, cp, crf_id, mid, field=None, types=None):
    d = LOG[LOG.gene == gene]
    d = d[(d.current_patient_id == cp) | (d.crf_patient_id == crf_id) | ((d.mention_id == mid) & (mid != ""))]
    if field:
        d = d[d.field == LOG_FIELD.get(field, field)]
    if types:
        d = d[d.change_type.isin(types)]
    return d.sort_values("change_id").tail(1)


def clean_ev(text):
    """변경 기록의 근거 문장에서 임상의에게 의미 없는 내부 표현(파일명, 규칙 번호, DB 용어)을 걷어 낸다."""
    t = str(text or "")
    _lab = {"Benefit": "효과 있음", "No benefit": "효과 없음", "Worsening": "악화", "unclear": "개별 환자 반응 확인 불가",
            "Used, response not reported": "사용했으나 반응 기록 없음", "Not used": "사용 기록 없음"}
    t = re.sub(r"^0904 통합본 입력값을 원문과 대조한 결과 다름 → 원문 값 유지\. ([A-Za-z ,]+): ",
               lambda m: f"원문 판단: {_lab.get(m.group(1), m.group(1))}. ", t)
    t = re.sub(r"^0904 통합본 입력값 원문 대조: [^→]+→ 원문 기준 값\. ", "", t)
    t = re.sub(r"^원문 대조\([^)]*\): ", "", t)
    t = re.sub(r"^원문 전수 확인(\([^)]*\))?: ", "", t)
    t = re.sub(r"re-extraction \S+\.json,?\s*\(?", "", t)
    t = re.sub(r"\b[A-Za-z0-9]+(?:_[A-Za-z0-9.+>-]+){2,}:\s*", "", t)          # 규칙 ID 접두어
    t = re.sub(r"\((?:extraction_error|[a-z]+(?:_[a-z]+)+)\)", "", t)          # (extraction_error) 등
    t = re.sub(r"Manual source review [0-9-]+\s*(\([^)]*\))?:?", "Source review:", t)
    t = re.sub(r"SCB response source review \([^)]*\):", "Source check of sodium channel blocker response:", t)
    t = re.sub(r":?\s*excluded from canonical(?: output)?\.?", " excluded.", t)
    t = re.sub(r"\S+\.json\S*", "", t)
    t = t.replace("row created by earlier extraction", "row in our earlier results")
    t = t.replace("CRF 통합본", "0904 통합본")
    t = re.sub(r"\bCRF\b", "0904 통합본", t)
    t = re.sub(r"\bDB\b", "previous version", t)
    t = re.sub(r"\bmentions?\b", "records", t)
    t = re.sub(r"\bcanonical\b", "final list", t)
    t = re.sub(r"cohort_index_(\d+)", r"row \1", t)
    return re.sub(r"\s{2,}", " ", t).strip(" ;,")


def reason_text(row):
    r = row.iloc[0]
    return batch_ko(r.batch_id), clean_ev(r.evidence or r.reason)


def norm_key(c, p):
    return cmp.norm_c(c) or cmp.norm_p(p)


_DEC_PATH = DECISIONS / "unmatched_decisions.csv"
CRF_DECISIONS = {}
if _DEC_PATH.exists():
    for _, _r in pd.read_csv(_DEC_PATH, dtype=str, keep_default_na=False).iterrows():
        CRF_DECISIONS[(_r.gene, _r.crf_patient_id)] = _r
# 원문 확인 후 개별로 사유를 적은 칸 (gene, 0904 Patient ID, 열)
NOTE_OVERRIDES = {
    **{("SCN2A", cp, "Age at seizure onset"): "[발병 나이] 원문의 구체적인 값으로 적었습니다. 0904 통합본의 범위 표기와 같은 범위입니다."
       for cp in ("CP0000907", "CP0000928", "CP0000947", "CP0000968")},
    ("SCN2A", "CP0000920", "Age at seizure onset"): "[발병 나이] 발작이 없는 환자라 발병 나이 칸을 비워 두었습니다(표현형 칸에 발작 없음으로 표시).",
    ("SYNGAP1", "SYNGAP1_004", "Response to targeted ASM"): "[표적 치료 반응] 같은 뜻이며, 원문 표현(발프로산으로 발작 완전 조절)을 그대로 적었습니다.",
    ("SYNGAP1", "SYNGAP1_009", "Response to targeted ASM"): "[표적 치료 반응] 같은 뜻이며, 원문 표현(발프로산 효과 있음)을 그대로 적었습니다.",
    **{("SYNGAP1", cp, "Mutation Type"): "[Mutation Type] 원문상 작은 frameshift 변이라 Large CNV가 아니어서 Small mutation으로 적었습니다."
       for cp in ("SYNGAP1_333", "SYNGAP1_334")},
}
INCLUSIONS = pd.read_csv(ROOT / "config/manual_article_inclusions.csv", dtype=str, keep_default_na=False)


def mention_status(gene, pmid, c, p):
    """CRF에만 있는 행: DB에 같은 논문·변이 언급이 있으면 그 분류로 사유를 만든다."""
    if not pmid:
        return "PMID가 없는 미발표 코호트 행입니다. 출판 문헌을 대상으로 하는 자료라 이번 수정본에는 넣지 않았습니다.", ""
    inc = INCLUSIONS[(INCLUSIONS.gene == gene) & (INCLUSIONS.pmid == pmid) & INCLUSIONS.status.str.startswith("not_added")]
    if len(inc):
        return {"not_added_no_patient": "원문을 확보해 확인한 결과 이 논문에는 해당 변이를 가진 환자가 없습니다",
                "not_added_prior_report": "원문을 확보해 확인한 결과 이 논문은 다른 논문의 환자를 인용한 리뷰입니다"}.get(
                    inc.iloc[0].status, "원문을 확보해 확인했습니다") + f" ({inc.iloc[0].notes}). 확인 부탁드립니다.", ""
    m = MENTIONS[gene]
    m = m[m.pmid == pmid]
    if m.empty:
        in_db = pmid in TITLES
        return ("이 논문을 추출했지만 개별 환자 행이 나오지 않았습니다. 원문 확인이 필요합니다." if in_db
                else "검색식(유전자 AND epilepsy)에 걸리지 않는 논문에서 0904 통합본에 추가하신 행입니다. 이번 수정본에는 아직 넣지 않았습니다. "
                     "포함 여부를 정해 주십시오."), ""
    kc, kp = cmp.norm_c(c), cmp.norm_p(p)
    hit = m[[bool((kc and cmp.norm_c(a) == kc) or (kp and cmp.norm_p(b) == kp)) for a, b in zip(m.variant_c, m.variant_p)]]
    if hit.empty:
        return "같은 논문은 추출했지만 이 환자·변이 행을 찾지 못했습니다. 원문 확인이 필요합니다.", ""
    h = hit.iloc[0]
    if h.canonical_status == "included":
        return (f"이번 수정본에도 있는 환자입니다(논문 내 {h.patient_id}). 변이 표기가 달라 같은 행으로 짝지어지지 않았고, "
                f"같은 논문의 초록색 추가 행과 같은 환자일 수 있습니다."), h.location_in_paper
    txt = EXCLUSION_KO.get(h.exclusion_reason, SOURCE_TYPE_KO.get(h.mention_source_type, "원문을 다시 확인해 제외했습니다."))
    if h.cited_source and h.cited_source.lower() != "not specified":
        txt += f" (원 논문: {h.cited_source})"
    return txt, clean_ev(h.manual_exclusion_evidence or h.originality_evidence)


def label_only(gene, field, crf_v, cur_v):
    return field == "phenotype" and LABEL.get((gene, str(cur_v).strip().lower())) == str(crf_v).strip().lower()


# ---------------------------------------------------------------- build records
# ---------------------------------------------------------------- diff workbook (CRF format, no extra columns)
import copy  # noqa: E402
import sys  # noqa: E402

from openpyxl import load_workbook  # noqa: E402
from openpyxl.comments import Comment  # noqa: E402

sys.path.insert(0, str(ROOT / "scripts" / "curation"))
import export_crf_workbook as ex  # noqa: E402

CRF_PATH = ROOT / "feedback/20260923/유전자별_Phenotype_CRF_통합_0904.xlsx"
OUT_DIR = ROOT / "outputs/current"  # 최신본 하나만 두고, 이전본은 backups/output_history/crf_diff/로 옮긴다
OUT_PATH = OUT_DIR / f"유전자별_Phenotype_CRF_통합_0904_변경대조_{__import__('datetime').datetime.now():%Y%m%d%H%M}.xlsx"
COL_KEY = {"Variant cDNA change": "variant_c", "Variant protein change": "variant_p", "ClinVar pathogenicity": "clinvar_pathogenicity",
           "ClinVar ID": "clinvar_id", "Sex": "sex", "Age at seizure onset": "age_of_onset", "Age at last follow-up": "last_followup_age",
           "Phenotype": "phenotype", "Response to targeted ASM": "targeted_response"}
COL_KO = {"Variant cDNA change": "cDNA 변이", "Variant protein change": "단백질 변이", "Mutation Type": "Mutation Type",
          "ClinVar pathogenicity": "ClinVar 분류", "ClinVar ID": "ClinVar ID", "Source in PMID": "PMID", "Sex": "성별",
          "Age at seizure onset": "발병 나이", "Age at last follow-up": "마지막 추적 나이", "Phenotype": "표현형",
          "Phenotype_subcategory": "표현형 하위분류", "Response to targeted ASM": "표적 치료 반응", "Patient ID": "Patient ID"}
OLD_ROW, OLD_CELL = PatternFill("solid", fgColor="FDE9EA"), PatternFill("solid", fgColor="F4A6AD")
NEW_ROW, NEW_CELL = PatternFill("solid", fgColor="E8F5E9"), PatternFill("solid", fgColor="93D3A2")
MOD_ROW, MOD_CELL = PatternFill("solid", fgColor="FFF8DC"), PatternFill("solid", fgColor="FFD966")
DV = {"Mutation Type": "_lists!$H$2:$H$3", "Sex": "_lists!$A$2:$A$3"}
PHEN_LIST = {"SCN1A": "_lists!$B$2:$B$5", "SCN2A": "_lists!$C$2:$C$8", "KCNQ2": "_lists!$D$2:$D$9", "CDKL5": "_lists!$E$2:$E$3",
             "SYNGAP1": "_lists!$F$2:$F$7", "EEF1A2": "_lists!$G$2:$G$8"}


def same(col, a, b):
    if col == "Source in PMID":
        return cmp.first_pmid(a) == cmp.first_pmid(b)
    if col == "Patient ID":
        return True  # 짝지어진 행은 같은 환자 (SYNGAP1/EEF1A2는 CRF 자체 ID 유지)
    key = COL_KEY.get(col)
    if key:
        if key == "phenotype" and label_only("", key, a, b):
            return True
        return cmp.norm_field(key, a) == cmp.norm_field(key, b)
    return cmp.norm_text(a) == cmp.norm_text(b)


def note(text):
    c = Comment(text[:1800], "변경 대조")
    c.width, c.height = 420, 180
    return c


_NOT_ADOPTED = {}
CHECK_VERDICT = {}
for _f, _why in [("crf_values_inherited_not_adopted.csv", None), ("crf_values_not_adopted_allele_ids.csv", "allele"),
                 ("crf_values_conflicting_with_source_corrections.csv", "conflict")]:
    if (DECISIONS / _f).exists():
        for _, _r in pd.read_csv(DECISIONS / _f, dtype=str, keep_default_na=False).iterrows():
            _k = _r.get("field", "")
            _NOT_ADOPTED[(_r.gene, _r.mention_id, _k)] = _why or _r.get("inherited", "")
            if _why == "conflict":
                CHECK_VERDICT[(_r.gene, _r.mention_id, _k)] = _r.get("check_verdict", "")
CUR_VARIANT = {}
for _g in GENES:
    _cur = pd.read_csv(sorted(glob.glob(str(ROOT / f"outputs/current/canonical_patients/{_g}_canonical_patients/*.csv")))[-1],
                       dtype=str, keep_default_na=False, encoding="utf-8-sig")
    CUR_VARIANT.update({(_g, m): v for m, v in zip(_cur.mention_id, _cur.variant_c)})
CLINVAR_IDS = {}
CLINVAR_STATUS = {}
for _g in GENES:
    _con = sqlite3.connect(ROOT / "databases" / f"{_g}.sqlite")
    for _mid, _t, _cls in _con.execute("select mention_id, clinvar, clinvar_classification from patient_mentions"):
        _t = _t or ""
        CLINVAR_STATUS[(_g, str(_mid))] = (_t, _cls or "")
        CLINVAR_IDS[(_g, str(_mid))] = (set(re.findall(r"AlleleID=(\d+)", _t)),
                                        {str(int(v)) for v in re.findall(r"VCV0*(\d+)", _t)})
RULES = pd.read_csv(ROOT / "config/manual_value_corrections.csv", dtype=str, keep_default_na=False)
FIELD_FROM_RULE = {"variant_type": "Mutation Type", "phenotype_review_category": "Phenotype", "clinvar_id": "ClinVar ID",
                   "targeted_treatment_response": "Response to targeted ASM"}
NOT_ADOPTED_TEXT = {
    "바탕 DB 값과 같음": "(0904 통합본의 이 값은 저희가 이전에 보내드린 추출 결과의 값이 그대로 옮겨진 것입니다.)",
    "ASM 반응 칸에서 옮긴 값": "0904 통합본에는 전체 ASM 반응(예: Refractory)이 적혀 있습니다. 이번 수정본에는 원문에서 찾은 표적 치료 반응 서술을 표시했습니다.",
    "사용 약물 칸에서 옮긴 값": "0904 통합본에는 사용 약물 목록이 적혀 있습니다. 이번 수정본에는 원문에서 찾은 표적 치료 반응 서술을 표시했습니다.",
    "allele": "0904 통합본의 번호는 예전 ClinVar Allele ID입니다. 현재는 ClinVar Variation ID로 표시합니다.",
    "conflict": "0904 통합본 값과 원문을 다시 확인한 결과가 다릅니다. 원문 값을 표시했으며, 확인 부탁드립니다."}
BASE = {}
for _g in GENES:
    _f = sorted(ROOT.glob(f"archive/backups_legacy/output_archives/previous_active_canonical_20260923_*/{_g}_canonical_patients/*_20260616_161440.csv"))
    BASE[_g] = pd.read_csv(_f[-1], dtype=str, keep_default_na=False).drop_duplicates("mention_id").set_index("mention_id") if _f else pd.DataFrame()
BASE_COL = {"Variant cDNA change": "variant_c", "Variant protein change": "variant_p", "Sex": "sex",
            "Age at seizure onset": "age_of_onset", "Age at last follow-up": "last_followup_age",
            "Phenotype": "phenotype_review_category", "ClinVar ID": "ClinVar ID", "ClinVar pathogenicity": "clinvar_pathogenicity",
            "Response to targeted ASM": "targeted_treatment_response"}


def cleared_by_clinician(gene, col, mid, cur_v):
    """0904 통합본 칸이 비어 있고, 통합본의 바탕 DB(6/16)에는 같은 값이 있었으면 임상의가 비운 것."""
    b = BASE.get(gene)
    if b is None or b.empty or mid not in b.index or col not in BASE_COL:
        return False
    bv = b.at[mid, BASE_COL[col]]
    return not cmp.missing(bv) and cmp.norm_text(bv) == cmp.norm_text(cur_v)


def not_adopted(gene, col, mid):
    for k, why in _NOT_ADOPTED.items():
        if k[0] == gene and k[1] == mid and (k[2] == COL_KEY.get(col) or FIELD_FROM_RULE.get(k[2]) == col):
            return NOT_ADOPTED_TEXT.get(why, why)
    return None


def hide_clinvar(gene, mid, crf_v, cur_v):
    """ClinVar 칸 차이 중 표시하지 않을 것: 0904 이후 추가 환자, 이전 조회에서 못 찾은 기록(변이 같음), 재조회로 값 변경.
    표시하는 것: 변이 표기를 바로잡아 새로 맞춰진 기록, 기록이 여러 건이라 비운 칸."""
    if cmp.missing(cur_v):
        return False
    b = BASE.get(gene)
    if b is None or b.empty or mid not in b.index:
        return True
    same_variant = cmp.norm_c(b.at[mid, "variant_c"]) == cmp.norm_c(CUR_VARIANT.get((gene, mid), ""))
    return same_variant


def allele_to_variation(gene, mid, crf_v):
    """0904 통합본 ClinVar ID가 이 변이 기록의 Allele ID이면 True (현재는 같은 기록의 Variation ID로 표기)."""
    alle, _ = CLINVAR_IDS.get((gene, mid), (set(), set()))
    ids = {x.strip().split(".")[0] for x in str(crf_v).split(";") if x.strip()}
    return bool(ids) and ids <= alle


def field_reason(gene, col, crf_v, cur_v, cp, crf_id, mid):
    key = COL_KEY.get(col)
    label = f"[{COL_KO[col]}]"
    if col == "Phenotype_subcategory" and (gene, mid, "phenotype_review_category") in CHECK_VERDICT:
        # 하위분류는 표현형 분류에서 나온다: 표현형 입력값을 원문 대조로 바꾼 결과
        verdict = CHECK_VERDICT[(gene, mid, "phenotype_review_category")]
        r = RULES[(RULES.gene == gene) & (RULES.mention_id == mid) & (RULES.field == "phenotype_review_category") & (RULES.status == "active")]
        ev = clean_ev(r.iloc[-1].evidence) if len(r) else ""
        why = {"not_reported": "원문에 이 환자의 증후군명이 없습니다",
               "partial": "원문에 이 증후군명이 그대로 적혀 있지 않습니다"}.get(verdict, "원문에 적힌 증후군과 다릅니다")
        return (f"{label} 0904 통합본에 입력하신 하위분류('{str(crf_v).strip()[:80]}')는 {why}. 원문 기준으로 표시했습니다. 확인 부탁드립니다."
                + (f"\n  근거: {ev[:300]}" if ev else ""))
    na = not_adopted(gene, col, mid)
    lg = log_for(gene, cp, crf_id, mid, key) if key else LOG.iloc[0:0]
    if na == NOT_ADOPTED_TEXT["conflict"]:
        ev = clean_ev(lg.iloc[0].evidence or lg.iloc[0].reason) if len(lg) else ""
        if not ev:
            rule_field = {v: k for k, v in FIELD_FROM_RULE.items()}.get(col, COL_KEY.get(col, ""))
            r = RULES[(RULES.gene == gene) & (RULES.mention_id == mid) & (RULES.field == rule_field) & (RULES.status == "active")]
            ev = clean_ev(r.iloc[-1].evidence) if len(r) else ""
            if col == "Mutation Type":
                ev = "0904 통합본 헤더 기준으로 엑손 단위 이상의 결실·중복은 Large CNV입니다. " + ev
        rule_field = {v: k for k, v in FIELD_FROM_RULE.items()}.get(col, COL_KEY.get(col, ""))
        verdict = next((v for (g_, m_, f_), v in CHECK_VERDICT.items() if g_ == gene and m_ == mid and f_ in (rule_field, COL_KEY.get(col))), "")
        head = {
            "partial": "은 원문과 일부만 일치합니다. 원문에 맞게 보완한 값을 표시했습니다.",
            "not_reported": "의 근거를 원문에서 찾지 못했습니다(원문에 이 환자의 해당 정보가 없습니다). 원문 기준으로 표시했습니다.",
            "cohort_level_only": "은 원문에서 이 환자 개인이 아니라 여러 환자 전체에 대한 서술에서 나온 값입니다. 원문 기준으로 표시했습니다.",
            "patient_not_found": "의 근거가 되는 환자를 원문에서 특정할 수 없었습니다. 원문 기준으로 표시했습니다.",
        }.get(verdict, "과 원문 확인 결과가 다릅니다. 원문 근거에 따라 이번 수정본에는 원문 값을 표시했습니다.")
        return (f"{label} 0904 통합본에 입력하신 값('{str(crf_v).strip()[:80]}'){head} 확인 부탁드립니다."
                + (f"\n  근거: {ev[:300]}" if ev else ""))
    if cmp.missing(crf_v) and cleared_by_clinician(gene, col, mid, cur_v):
        return f"{label} 0904 통합본에는 비어 있던 칸입니다. 원문을 다시 확인해 찾은 값을 표시했습니다."
    if col == "Source in PMID":
        return (f"{label} 0904 통합본은 이 환자를 인용한 논문(PMID {cmp.first_pmid(crf_v)})으로 적혀 있으나, "
                f"원문을 확인해 이 환자를 처음 보고한 논문(PMID {cmp.first_pmid(cur_v)})으로 바꿨습니다.")
    if col == "ClinVar ID" and not cmp.missing(crf_v) and not cmp.missing(cur_v):
        alle, vcv = CLINVAR_IDS.get((gene, mid), (set(), set()))
        crf_ids = {x.strip().split(".")[0] for x in str(crf_v).split(";") if x.strip()}
        if crf_ids and crf_ids <= alle:
            return (f"{label} 0904 통합본의 번호는 ClinVar Allele ID(저희가 이전에 보내드린 자료의 표기 오류)이고, 이번 수정본은 "
                    f"같은 기록의 Variation ID(ClinVar 화면의 VCV 번호)로 표기했습니다. 가리키는 ClinVar 기록은 같습니다.")
    if col.startswith("ClinVar"):
        if cmp.missing(cur_v) and not cmp.missing(crf_v):
            if na == NOT_ADOPTED_TEXT["allele"]:
                return f"{label} {na}"
            what = "번호" if col == "ClinVar ID" else "분류"
            text, cls = CLINVAR_STATUS.get((gene, mid), ("", ""))
            if text.startswith("not_listed") or cls == "not in ClinVar":
                return (f"{label} 0904 통합본의 {what}는 저희가 이전 조회에서 잘못 붙인 다른 변이의 ClinVar 기록이었습니다. "
                        f"이 변이와 정확히 일치하는 ClinVar 기록은 없어 비웠습니다.")
            if cls:
                return (f"{label} 이 변이의 ClinVar 기록은 현재 '{cls}'로 Pathogenic/Likely pathogenic이 아니어서, "
                        f"0904 통합본 규칙(P/LP일 때만 표시)에 따라 비웠습니다. 0904 통합본의 {what}는 다른 변이 기록이 섞인 이전 조회 결과였습니다.")
            return f"{label} 0904 통합본의 {what}는 저희가 이전에 보내드린 추출 결과의 값입니다. ClinVar에서 이 변이 기록을 다시 확인하고 있습니다."
        if cmp.missing(crf_v):
            b = BASE.get(gene)
            cur_c = CUR_VARIANT.get((gene, mid), "")
            if b is not None and not b.empty and mid in b.index:
                old_c = b.at[mid, "variant_c"]
                if cmp.norm_c(old_c) != cmp.norm_c(cur_c):
                    return (f"{label} 0904 통합본에는 없던 정보입니다. 변이 표기를 원문대로 바로잡은 뒤({old_c} → {cur_c}) "
                            f"ClinVar 기록과 일치해 추가했습니다.")
                return f"{label} 0904 통합본에는 없던 정보입니다. 이전 조회에서 찾지 못했던 ClinVar 기록을 다시 조회해 찾아 추가했습니다."
            return f"{label} 0904 통합본에는 없던 정보입니다. ClinVar를 조회해 추가했습니다."
        b = BASE.get(gene)
        if b is not None and not b.empty and mid in b.index:
            old_c, cur_c = b.at[mid, "variant_c"], CUR_VARIANT.get((gene, mid), "")
            if cmp.norm_c(old_c) != cmp.norm_c(cur_c):
                return f"{label} 변이 표기를 원문대로 바로잡은 뒤({old_c} → {cur_c}) 해당 변이의 ClinVar 기록으로 바꿨습니다."
        return f"{label} ClinVar를 다시 조회해 최신 정보로 바꿨습니다."
    if len(lg):
        action, ev = reason_text(lg)
        ev = "" if "variant_prefix_fix" in lg.iloc[0].batch_id else ev
        head = ""
        if cmp.missing(crf_v):
            head = "0904 통합본에는 없던 정보입니다. "
            action = action.replace("바로잡았습니다", "찾아 추가했습니다").replace("반영했습니다", "찾아 추가했습니다")
        text = f"{label} {head}{action}"
        if na:
            text += f" {na}"
        return text + (f"\n  근거: {ev[:300]}" if ev else "")
    if na:
        return f"{label} {na}"
    if cmp.missing(crf_v):
        new_paper = ((LOG.gene == gene) & (LOG.mention_id == mid) & LOG.batch_id.str.contains("crf_added_papers")).any()
        how = ("이후 ClinVar를 다시 조회해 추가했습니다." if col.startswith("ClinVar") else
               "이 논문 원문을 새로 확보(기관 접속)해 확인하고 추가했습니다." if new_paper else
               "이후 원문을 다시 추출·확인하는 과정에서 찾아 추가했습니다.")
        return f"{label} 0904 통합본에는 없던 정보입니다. {how}"
    if cmp.missing(cur_v):
        return f"{label} 0904 통합본에는 값이 있으나 이번 수정본에는 비어 있습니다. 확인하겠습니다."
    return f"{label} 0904 통합본 값과 이번 수정본 값이 다릅니다. 변경 이력이 확인되지 않아 다시 확인하겠습니다."


def diff_sheet(ws, gene, label_map):
    header = [c.value for c in ws[2]]
    ncol = len(header)
    crf, can = cmp.load(gene)
    pairs, crf_only, can_only = cmp.match(gene, crf, can)
    # 원문 확인으로 같은 환자임이 확인된 행(변이 표기 차이로 자동 짝짓기 실패)을 짝지어 준다
    for i in list(crf_only):
        dec = CRF_DECISIONS.get((gene, str(crf.at[i, "Patient ID"])))
        if dec is None or dec.outcome != "same_as_existing" or not dec.our_canonical_patient_id:
            continue
        js = [j for j in can_only if can.at[j, "canonical_patient_id"] == dec.our_canonical_patient_id.split(";")[0].strip()]
        if js:
            pairs.append((i, js[0], "source_verified_pair"))
            crf_only.remove(i)
            can_only.remove(js[0])
    cur_rows = ex.gene_rows(gene, can.reset_index(drop=True), label_map)
    pos = {idx: k for k, idx in enumerate(can.index)}
    # 원본 행 (엑셀 행 번호 = crf_row)
    orig = {r: [ws.cell(r, c) for c in range(1, ncol + 1)] for r in range(3, ws.max_row + 1)}
    vals = {r: [c.value for c in cells] for r, cells in orig.items()}
    links = {r: [c.hyperlink.target if c.hyperlink else None for c in cells] for r, cells in orig.items()}
    style = [copy.copy(ws.cell(3, c)._style) for c in range(1, ncol + 1)]
    # 0904 통합본 칸 서식은 행마다 원래대로 유지하되, 임상의가 칠한 색은 지운다(아래에서 fill 초기화).
    row_style = {r: [copy.copy(c._style) for c in cells] for r, cells in orig.items()}
    plain = []
    for c in range(1, ncol + 1):
        _c = ws.cell(3, c)
        _saved = copy.copy(_c._style)
        _c.fill = PatternFill(fill_type=None)
        plain.append(copy.copy(_c._style))
        _c._style = _saved
    by_row = {}
    verified = set()
    for i, j, how in pairs:
        by_row[int(crf.at[i, "crf_row"])] = ("pair", i, j)
        if how == "source_verified_pair":
            verified.add(int(crf.at[i, "crf_row"]))
    for i in crf_only:
        by_row[int(crf.at[i, "crf_row"])] = ("crf_only", i, None)
    added = defaultdict(list)
    for j in can_only:
        added[can.at[j, "pmid_n"]].append(j)
    last_row_of_pmid = {}
    for r in sorted(vals):
        last_row_of_pmid[cmp.first_pmid(vals[r][header.index("Source in PMID")])] = r

    out = []  # (values, kind, changed_cols, comment, links)
    stats = defaultdict(int)

    def cur_values(j):
        row = cur_rows[pos[j]]
        return [row.get(h, "") for h in header]

    def added_block(pmid):
        for j in added.pop(pmid, []):
            v = cur_values(j)
            cp, mid = can.at[j, "canonical_patient_id"], str(can.at[j, "mention_id"])
            lg = log_for(gene, cp, "", mid, types=["row_added"])
            if len(lg):
                r_ = batch_ko(lg.iloc[0].batch_id, ADD_KO)
                e_ = clean_ev(lg.iloc[0].evidence or lg.iloc[0].reason)
            else:
                r_ = "0904 통합본을 만들 때 빠져 있던 환자로, 원문에서 확인해 추가했습니다."
                e_ = clean_ev(can.at[j, "originality_evidence"])
            txt = f"[추가] 0904 통합본에 없던 환자입니다. {r_}\n  근거: {e_ or can.at[j, 'location_in_paper']}\n  논문 내 환자: {can.at[j, 'patient_id']}" + (
                "")
            out.append((v, "added", set(), txt, [None] * ncol, {}, None))
            stats["추가"] += 1

    for r in sorted(vals):
        kind, i, j = by_row.get(r, ("keep", None, None))
        if kind == "pair":
            new = cur_values(j)
            changed = {h for h, a, b in zip(header, vals[r], new) if not same(h, a, b)}
            mid = str(can.at[j, "mention_id"])
            base_vals = list(vals[r])
            if ("ClinVar ID" in changed and not cmp.missing(new[header.index("ClinVar ID")])
                    and (allele_to_variation(gene, mid, vals[r][header.index("ClinVar ID")])
                         or not_adopted(gene, "ClinVar ID", mid) == NOT_ADOPTED_TEXT["allele"])):
                # Allele ID -> 같은 기록의 Variation ID: 변경으로 표시하지 않고 값만 바꿔 적는다
                changed.discard("ClinVar ID")
                base_vals[header.index("ClinVar ID")] = new[header.index("ClinVar ID")]
                stats["ClinVar 번호 표기 정정"] += 1
            for col in ("ClinVar ID", "ClinVar pathogenicity"):
                # ClinVar는 외부 DB 재조회 값이라 임상의 검토 대상이 아님: 변경으로 표시하지 않고 현재 값만 적는다
                if col in changed and col in header:
                    # 이전 조회에서 못 찾은 기록 추가·재조회로 값 변경·추가 환자: 변경으로 표시하지 않고 현재 값만 적는다
                    changed.discard(col)
                    base_vals[header.index(col)] = new[header.index(col)]
            if not changed:
                out.append((base_vals, "keep", set(), None, links[r], {}, r))
            else:
                merged = [b if h in changed else a for h, a, b in zip(header, base_vals, new)]
                cp, crf_id = can.at[j, "canonical_patient_id"], crf.at[i, "Patient ID"]
                reasons = [NOTE_OVERRIDES.get((gene, str(crf_id), h)) or field_reason(gene, h, a, b, cp, crf_id, mid)
                           for h, a, b in zip(header, vals[r], new) if h in changed]
                if r in verified:
                    reasons = [x for x in reasons if not x.startswith(("[cDNA 변이]", "[단백질 변이]"))]
                    reasons.insert(0, "[변이 표기] 원문 확인 결과 이번 수정본의 같은 환자이며, 변이 표기만 다릅니다(원문 기준 표기로 적었습니다).")
                cid = next((x for x in reasons if x.startswith("[ClinVar ID]")), None)
                ccl = next((x for x in reasons if x.startswith("[ClinVar 분류]")), None)
                if cid and ccl and cid.split("] ", 1)[1].replace("번호", "값") == ccl.split("] ", 1)[1].replace("분류", "값"):
                    reasons = [x for x in reasons if x is not ccl]
                    reasons[reasons.index(cid)] = "[ClinVar ID·분류] " + cid.split("] ", 1)[1].replace("번호", "번호·분류")
                # 한 줄로 표시: 수정한 칸만 새 값(색 칸), 0904 통합본 값은 칸 메모와 비고에 남긴 채
                olds = {h: a for h, a in zip(header, vals[r]) if h in changed}
                for h, a in olds.items():
                    if not cmp.missing(a):
                        lab = f"[{COL_KO.get(h, h)}]"
                        k = next((n for n, x in enumerate(reasons) if x.startswith(lab)), None)
                        line = f"  0904 통합본 값: {str(a).strip()[:200]}"
                        if k is None:
                            reasons.append(f"{lab}\n{line}")
                        else:
                            head_, _, tail_ = reasons[k].partition("\n")
                            reasons[k] = head_ + "\n" + line + ("\n" + tail_ if tail_ else "")
                lk = [links[r][c] if h not in changed else None for c, h in enumerate(header)]
                out.append((merged, "mod", changed, "[수정]\n" + "\n".join(reasons), lk, olds, r))
                stats["수정"] += 1
        elif kind == "crf_only":
            crf_id = crf.at[i, "Patient ID"]
            lg = log_for(gene, crf_id, crf_id, "", types=["row_excluded", "merged_duplicate", "same_variant_row_collapsed"])
            if len(lg):
                r_ = batch_ko(lg.iloc[0].batch_id, DELETE_KO)
                e_ = clean_ev(lg.iloc[0].evidence or lg.iloc[0].reason)
            elif (gene, str(crf_id)) in CRF_DECISIONS and CRF_DECISIONS[(gene, str(crf_id))].outcome == "same_as_existing":
                same_cp = CRF_DECISIONS[(gene, str(crf_id))].our_canonical_patient_id.split(";")[0].strip()
                r_ = (f"원문 확인 결과 이 행의 환자는 {same_cp} 환자와 같은 사람입니다. 0904 통합본에 같은 환자가 "
                      f"{same_cp} 행으로도 들어 있어 한 행으로 정리했습니다.")
                e_ = clean_ev(CRF_DECISIONS[(gene, str(crf_id))].evidence)
            elif (gene, str(crf_id)) in CRF_DECISIONS and CRF_DECISIONS[(gene, str(crf_id))].outcome != "same_as_existing":
                dec = CRF_DECISIONS[(gene, str(crf_id))]
                src = crf.at[i, "pmid_n"]
                if dec.outcome == "previously_reported":
                    early = dec.earlier_pmid if re.search(r"\d{7,}", dec.earlier_pmid or "") else ""
                    r_ = (f"통합본에 적힌 논문(PMID {src})은 이 환자를 다른 논문" + (f"(PMID {early})" if early else "") +
                          "에서 인용한 것이라, 이 논문의 환자로는 넣지 않고 원래 보고한 논문 기준으로 정리했습니다.")
                elif dec.outcome == "not_individual":
                    r_ = "원문을 확인한 결과 환자가 아니라 실험실에서 만든 변이이거나 개별 환자 정보가 없는 행입니다."
                elif dec.outcome == "unclear":
                    r_ = (f"통합본에 적힌 논문(PMID {src})에서 이 변이를 찾지 못했습니다. 출처 논문을 알려주시면 반영하겠습니다.")
                else:
                    r_ = "원문을 다시 확인했습니다."
                e_ = clean_ev(dec.evidence)
            else:
                r_, e_ = mention_status(gene, crf.at[i, "pmid_n"], crf.at[i, "Variant cDNA change"], crf.at[i, "Variant protein change"])
            out.append((vals[r], "removed", set(header), f"[삭제] 이번 수정본에서 뺀 행입니다. {r_}" + (f"\n  근거: {e_}" if e_ else ""), links[r], {}, r))
            stats["CRF에만"] += 1
        else:
            out.append((vals[r], "keep", set(), None, links[r], {}, r))
        pm = cmp.first_pmid(vals[r][header.index("Source in PMID")])
        if last_row_of_pmid.get(pm) == r:
            added_block(pm)
    for pm in sorted(added, key=lambda x: int(x) if str(x).isdigit() else 10 ** 9):
        added_block(pm)

    # rewrite data rows
    ws.delete_rows(3, ws.max_row)
    for k, (v, kind, changed, text, lk, olds, src) in enumerate(out):
        r = 3 + k
        for c, (h, val) in enumerate(zip(header, v), 1):
            cell = ws.cell(r, c, value=val if val is not None and str(val) != "nan" else None)
            cell._style = copy.copy(row_style[src][c - 1] if src else plain[c - 1])
            cell.fill = PatternFill(fill_type=None)  # 임상의 색 표시는 대조본에서 지움 (변경 표시 색만 남김)
            cell.border = GRID  # 0904 통합본 EEF1A2 시트는 테두리가 빠진 칸이 많아 모든 시트에 같은 테두리를 적용
            if lk[c - 1]:
                cell.hyperlink = lk[c - 1]
            elif kind in ("mod", "added") and val and h == "ClinVar ID" and str(val).isdigit() and h in changed | ({h} if kind == "added" else set()):
                cell.hyperlink = f"https://www.ncbi.nlm.nih.gov/clinvar/variation/{val}/"
            elif kind == "added" and h == "Source in PMID" and str(val).isdigit():
                cell.hyperlink = f"https://pubmed.ncbi.nlm.nih.gov/{val}/"
            f = copy.copy(cell.font)
            if kind in ("old", "removed"):
                cell.fill = OLD_CELL if (h in changed and kind == "old") else OLD_ROW
                if h in changed:
                    cell.font = Font(name=f.name, sz=f.sz, color=f.color, bold=f.bold, strike=True)
            elif kind == "added":
                cell.fill = NEW_ROW
            elif kind == "mod":
                cell.fill = MOD_CELL if h in changed else MOD_ROW
                if h in changed:
                    cell.font = Font(name=f.name, sz=f.sz, color=f.color, bold=True)
                    a = olds.get(h)
                    cell.comment = note("0904 통합본 값: " + ("(빈칸)" if cmp.missing(a) else str(a).strip()))
        memo = ws.cell(r, ncol + 1, value=text)
        memo._style = copy.copy(style[0])
        memo.alignment = Alignment(wrap_text=True, vertical="top")
        memo.border = GRID
        if kind in ("old", "removed"):
            memo.fill = OLD_ROW
        elif kind == "added":
            memo.fill = NEW_ROW
        elif kind == "mod":
            memo.fill = MOD_ROW
    last = 2 + len(out)
    head = ws.cell(2, ncol + 1, value="비고")
    head._style = copy.copy(ws.cell(2, ncol)._style)
    ws.column_dimensions[get_column_letter(ncol + 1)].width = 70
    from openpyxl.worksheet.datavalidation import DataValidation as _DV
    ws.data_validations.dataValidation = []
    for h, formula in {**DV, "Phenotype": PHEN_LIST[gene]}.items():
        if h in header:
            col = get_column_letter(header.index(h) + 1)
            dv = _DV(type="list", formula1=formula, allow_blank=True)
            dv.add(f"{col}3:{col}{last}")
            ws.add_data_validation(dv)
    return stats


def main():
    wb = load_workbook(CRF_PATH)
    label_map = ex.load_label_map()
    summary = {}
    for gene in GENES:
        summary[gene] = diff_sheet(wb[gene], gene, label_map)
        legend = ("변경 대조 (0904 통합본 → 이번 수정본)\n"
                  "· 연노랑 행: 수정한 행. 진노랑·굵은 칸이 바뀐 칸(이번 수정본 값)이며, 0904 통합본 값은 칸 메모와 '비고'에 있습니다\n"
                  "· 연초록 단독 행: 이번 수정본에 새로 추가한 행\n"
                  "· 연빨강·취소선 단독 행: 이번 수정본에서 뺀 행\n"
                  "· 사유·근거: 맨 오른쪽 '비고' 열\n"
                  "· ClinVar ID·분류: ClinVar를 다시 조회한 현재 값(Variation ID 기준)으로 적었고, 변경으로 표시하지 않았습니다.")
        wb[gene]["A1"].comment = note(legend)
    for old in OUT_DIR.glob("유전자별_Phenotype_CRF_통합_0904_변경대조_*.xlsx"):
        shutil.move(str(old), str(ROOT / "backups/output_history/crf_diff" / old.name))
    wb.save(OUT_PATH)
    print(OUT_PATH)
    for g, s in summary.items():
        print(g, dict(s))


if __name__ == "__main__":
    main()
