# 큐레이션 정책

유전자별 phenotype 기준, 중복 제거, 수작업 정정·변이 분류·ClinVar 기준을 정리한 문서입니다.

## 유전자별 phenotype 기준

기본 DB column은 모든 유전자에서 같습니다. 즉, `variant_c_original`, `variant_c_gnomad`, `variant_p`, `patient_id`, `seizure_types`, `response_to_asm` 같은 컬럼 구조는 유전자별로 달라지지 않습니다. 검토용 CSV에서는 이 중 cDNA 관련 컬럼을 `variant_c_source`, `variant_c`로 합쳐 표시합니다.

다만 **phenotype 해석 기준과 prompt는 유전자별로 다릅니다.** 유전자마다 중요한 syndrome category와 약물 반응 포인트가 다르기 때문입니다.

유전자별 기준은 `config/gene_phenotype_review_guidance.csv`가 source of truth입니다. `scripts/core/batch_extract_patient_mentions.py`는 이 CSV를 읽어 prompt를 생성합니다.

| gene | phenotype category 예시 | 중점적으로 보는 치료 반응 |
| --- | --- | --- |
| `SCN1A` | Dravet syndrome, GEFS+, Other epilepsy, Non-seizure | sodium channel blocker 악화/효과, stiripentol, cannabidiol, fenfluramine |
| `SCN2A` | early infantile epileptic encephalopathy/Ohtahara, benign/self-limited epilepsy, Other epilepsy | sodium channel blocker 효과 |
| `KCNQ2` | developmental and epileptic encephalopathy, West syndrome, benign/self-limited epilepsy | sodium channel blocker 효과 |
| `CDKL5` | Epilepsy, Non-seizure types | ketogenic diet, cannabidiol, fenfluramine |
| `EEF1A2` | DEE, Lennox-Gastaut syndrome, West syndrome, Dravet syndrome, epileptic encephalopathy | valproic acid 효과 |
| `SYNGAP1` | Generalized epilepsy, myoclonic astatic/atonic epilepsy, early onset absence epilepsy | valproic acid 효과 |
| `PCDH19` | Dravet syndrome, Other epilepsy, Non-seizure | 효과 있는 ASM |
| `SCN8A` | benign/self-limited epilepsy, epileptic encephalopathy/DEE, Non-seizure | sodium channel blocker 효과 |
| `GABRA1` | Dravet syndrome, DEE/EE, Other epilepsy, Non-seizure | 효과 있는 ASM |
| `STXBP1` | EIEE/Ohtahara, West/infantile spasms, Other epilepsy, Non-seizure | 효과 있는 ASM |
| `SMC1A` | DEE/EE, Other epilepsy, Non-seizure | 효과 있는 ASM |
| `FOXG1` | DEE/EE, Other epilepsy, Non-seizure | 효과 있는 ASM |
| `KCNT1` | epilepsy of infancy with migrating focal seizures, ADNFLE, Other epilepsy, Non-seizure | 효과 있는 ASM |

예를 들어 같은 `response_to_asm` 컬럼이라도, `SCN1A`에서는 sodium channel blocker가 악화시켰는지 여부가 중요하고, `SYNGAP1`에서는 valproic acid benefit이 명시되어 있는지가 중요합니다.

2026-09-21 보완 요청에 따라 새 7개 유전자(`PCDH19`, `SCN8A`, `GABRA1`, `STXBP1`, `SMC1A`, `FOXG1`, `KCNT1`) 기준도 같은 config에 추가되어 있습니다.


## 중복 제거 정책

추출 단계에서는 deduplication을 하지 않습니다.

이유:

1. 같은 환자가 여러 논문에 반복 등장할 수 있습니다.
2. review 논문이나 cohort summary에서 기존 환자를 다시 표로 정리하는 경우가 있습니다.
3. 추출 시점에 자동으로 합치면 원문 evidence를 잃거나 서로 다른 환자를 잘못 합칠 수 있습니다.

따라서 `patient_mentions`에는 모든 mention을 보존합니다. 최종 unique patient DB는 별도 후처리로 만듭니다.

Deduplication 원칙:

- `original_patient`를 중심으로 canonical patient 후보를 만듭니다.
- `cited_patient`는 중복 확인 evidence로 사용합니다.
- 같은 variant라도 patient ID, 논문, phenotype, cited source를 함께 보고 판단합니다.
- 자동 병합은 보수적으로 수행합니다.
- cross-paper original patient 병합은 자동으로 넓히지 않습니다. 같은 환자임이 문헌 근거로 확인된 경우에만 `config/manual_merge_rules.csv`에 `mention_id` 묶음을 기록하고, `scripts/core/deduplicate_patients.py`가 해당 규칙만 수동 병합합니다.
- 한 환자에게 여러 변이가 보고되었더라도 dominant gene에서는 causal/P/LP variant가 보통 1개라는 전제를 적용합니다. causal variant가 아닌 secondary variant는 원문 근거가 있을 때 `config/manual_variant_exclusions.csv`로 canonical에서 제외합니다.
- `previously reported`, `previously described`, `re-analyzed/re-analysed` 등 기존 보고 또는 재분석 환자를 시사하는 행은 strict source filter에서 canonical 후보에서 제외합니다.
- Methods/Supplement/Table note에 기존 cohort 재사용이 명시된 경우, 새 환자로 명확히 구분되는 subset만 `original_patient`로 유지합니다.
- 환자별 임상정보 없이 variant list만 제공된 row도 canonical 후보에서 제외합니다.
- Individual genotype은 있지만 phenotype이 cohort-level summary로만 제공되어 직접 매칭이 안 되는 row는 `variant_only` 또는 `unclear`로 분리하고 canonical 후보에서 제외합니다.
- 단순히 같은 유전자, 같은 variant, 같은 sex, 같은 onset이라는 이유만으로는 cross-paper 병합하지 않습니다. recurrent variant 때문에 false-positive merge가 생길 수 있기 때문입니다.
- 같은 환자임을 뒷받침하는 근거는 원문 문구, cited source, PMID mapping, 같은 PDF/source file, 같은 family/patient identifier 등 source-supported evidence여야 합니다.
- `{GENE}_canonical_id_locks.csv`로 기존 CP 번호를 유지합니다.
- 병합된 한 환자 안에서 같은 변이가 논문마다 다르게 표기된 경우(cDNA 누락, `c.` 접두어 유무, frameshift 표기 차이)는 cDNA 또는 단백질 핵심 표기가 같으면 한 행으로 묶습니다.
- 변이 분류 기준 (2026-09-23 결정, 교수님 통합본의 실제 기준과 일치): VUS(의미 불명 변이)만 있는 환자는 canonical에 포함하되 VUS임을 표시하고, 논문이 benign/polymorphism으로 분류한 변이만 있는 환자는 제외합니다. 원인 변이가 따로 있는 환자의 benign/다형성 변이 행은 `manual_variant_exclusions.csv`로 제외합니다. 기존 6개 유전자에서 통합본에 남은 VUS 환자(예: 29141279 pt16–23)는 유지합니다.
- `manual_merge_rules.csv`의 `preferred_mention_id`에 mention을 지정하면 병합 환자의 canonical 행에 그 mention의 값을 표시합니다. 비워 두면 추출 정보가 가장 많은 mention을 자동 선택합니다.

## 원문 확인 원칙

- 교수님·임상의 피드백이나 통합본과 값이 다르면, 반영 전에 논문 원문(본문 표, supplementary)으로 어느 쪽이 맞는지 먼저 확인합니다.
- 원문 근거는 정정 규칙과 변경 기록의 `evidence`에 인용문·위치로 남깁니다.
- 논문 자체의 오타는 `variant_c_original`에 원문 그대로 두고, 비교용 `variant_c_gnomad`/`variant_p`만 기준 transcript로 정정합니다.

## 수작업 정정 (`config/manual_value_corrections.csv`)

- 원문 확인으로 고친 값은 DB를 직접 고치는 것으로 끝내지 않고 이 규칙 파일에 한 행씩 남깁니다 (gene, mention_id, pmid, patient_id, field, value, previous_value, evidence, source_run).
- `scripts/curation/apply_manual_corrections.py`가 규칙을 DB에 적용합니다. 이미 같은 값이면 건너뛰므로 여러 번 실행해도 안전합니다.
- 재추출로 mention_id가 바뀌면 해당 규칙은 `unresolved`로 보고되며, 새 mention_id로 규칙을 고쳐야 합니다.
- `rebuild_canonical.py`는 빌드 전에 이 규칙을 자동 적용합니다.

## CP 번호 고정 (`config/canonical_id_locks/{GENE}.csv`)

- mention_id → canonical_patient_id를 고정해, 재빌드·재추출 후에도 기존 환자 번호가 유지되게 합니다.
- `rebuild_canonical.py --promote`가 교체 시 새 환자 번호를 lock에 추가합니다.
- 병합할 때는 canonical에 남는(표시되는) 행의 기존 CP를 유지하고 다른 쪽 CP만 없앱니다. 교수님 통합본과 CP 번호로 대조할 수 있게 하기 위해서입니다.
- 통합본에만 있던 환자를 복원할 때는 통합본이 부여한 CP 번호를 lock으로 지정합니다.

## 병합 규칙 (`config/manual_merge_rules.csv`)

- 같은 환자라는 원문 근거(같은 기관·저자·코호트 + 환자 정보 일치, 또는 "이전 보고" 명시)가 있을 때만 병합합니다. 같은 변이만으로는 병합하지 않습니다.
- `preferred_mention_id`를 지정하면 병합 환자의 canonical 행에 그 mention의 값을 표시합니다. 비우면 정보가 가장 많은 mention을 자동 선택합니다.
- 병합된 한 환자 안에서 같은 변이가 논문마다 다르게 표기된 경우(cDNA 누락, `c.` 유무, frameshift 표기, 삭제 염기 표기, 서술형 exon CNV)는 한 행으로 묶습니다.

## 변이 분류 기준 (2026-09-23 결정)

- VUS만 있는 환자: canonical에 포함하고 VUS임을 표시합니다. 표시 컬럼 `variant_classification_flag`(CRF export의 `VUS flag`)는 논문 분류(`variant_classification_in_paper`) 또는 ClinVar 원 분류(`clinvar_classification`)가 VUS/상충이면 채워집니다. 두 컬럼은 2026-09-23 이후 추출분부터 생기므로 기존 6개 유전자에는 아직 없습니다.
- 논문이 benign/polymorphism으로 분류한 변이 행은 자동 제외합니다 (`exclusion_reason=benign_or_polymorphism_in_paper`). 그 변이만 있는 환자는 결과적으로 제외됩니다.
- 원인 변이가 따로 있는 환자의 benign/다형성 변이 행: `config/manual_variant_exclusions.csv`로 제외합니다.
- 여러 유전자를 함께 포함하는 대형 CNV(다른 원인 유전자 포함)는 해당 유전자의 원인 변이로 보지 않습니다.

## 유전형 최초 보고 (2026-09-24 결정, 개선 항목 8)

임상 정보는 이전 논문(예: 연관 분석·임상 보고)에 실렸어도, 해당 환자의 변이를 처음 보고한 논문의 행을 original로 봅니다. 이전 임상 논문은 `cited_source`에 남기고 `originality_evidence`에 "first genotype report"를 적습니다 (이 문구가 있으면 이미 보고된 환자 필터에서 제외되지 않음). 이전 논문이 같은 환자의 변이를 이미 보고했다면 previously_reported입니다. 이유: 이전 논문에는 변이가 없어 그 행도 canonical에 들어갈 수 없으므로, 규칙이 없으면 환자가 DB에서 사라짐. 적용 예: 12243921 가족 B (Lewis 1996 임상 보고) — 교수님 통합본도 이 논문 환자로 수록.

## 증후군 정의 코호트 (2026-09-23 결정)

논문이 환자별로 변이만 제시하고, 코호트 전체가 하나의 증후군(예: 모든 환자 Dravet syndrome)으로 선정된 경우:
- 환자 행을 유지하고 `phenotype_review_category`는 그 증후군으로 둡니다 (교수님 통합본 기준과 동일).
- 코호트 통계·포함 기준에서 복사된 값(발병 나이 중앙값, 발작 유형, 치료 반응 등)은 "Not specified"로 비웁니다.
- `other_phenotypes`에 "Syndrome-defined cohort (...): variant listed per patient, no individual clinical details in paper"로 표시합니다.
- 논문 목록과 근거는 `config/syndrome_defined_cohorts.csv`.
- 여러 증후군이 섞인 코호트처럼 환자별 표현형을 정할 수 없으면 제외합니다 (예: 32853554, 41764152).

## 자폐·지적장애 코호트 (2026-09-25)

환자별 뇌전증 표현형 없이 자폐 또는 지적장애로만 선정된 코호트의 변이 행은 canonical에 넣지 않습니다 (41764152 결정과 동일, 25418537·25363768 적용). 통합본에 "Non-seizure"로 들어간 경우는 회신에서 차이로 알립니다.

## 다유전자 CNV 판단 (2026-09-25 보강)

결실·중복이 다른 뇌전증 원인 유전자(SCN1A·SCN2A·SCN3A·SCN9A·CHRNA4·EEF1A2 등)를 함께 포함하면 원인 변이로 보지 않고 제외합니다. 대상 유전자만 포함하거나 범위가 불명확한 단일 유전자 결실은 유지합니다.

## 같은 코호트 논문 간 중복 (2026-09-25)

같은 기관 코호트를 여러 논문이 다루면, 변이를 먼저 보고한 논문의 행을 원래 환자로 두고 나중 논문의 행은 previously_reported로 돌립니다 (예: CDKL5 Demarest 2019 = 31313283 기준으로 34028805·37771170 정리). 재발 변이라 한 명으로 특정할 수 없는 경우는 그대로 둡니다.

## 변이가 기록되지 않은 행

cDNA·protein이 모두 비어 있는 original 행은 canonical에서 빠지며, 예전에는 export 단계에서 조용히 사라졌습니다. 2026-09-23부터 excluded 파일에 `exclusion_reason=variant_not_recorded`로 남겨 원문 재확인 대상으로 봅니다 (당일 기준 6개 유전자 1,224행).

## ClinVar

- `clinvar_id`는 Variation ID(VCV accession의 숫자, 예: VCV000192252 → 192252)입니다. 2026-09-23 이전에는 Allele ID가 저장되어 있었고 전 유전자 정정했습니다.
- `clinvar_pathogenicity`는 Pathogenic / Likely pathogenic / Pathogenic/Likely pathogenic만 표시합니다.

## 변경 기록 (`outputs/current/curation_change_log.csv`)

- canonical 데이터를 바꿀 때마다 교수님 통합본 Patient ID 기준으로 행 단위 기록을 남깁니다 (`rebuild_canonical.py --promote --batch-id ...`가 자동 추가).
- 변경 유형: 값 정정, 중복 환자 병합, 제외, 추가(복원), 같은 변이 행 통합, 병합된 행의 DB 값 정정.
- 교수님 회신용 엑셀은 `runs/reply_packages/reply_1/build_change_log_workbook.py`로 만듭니다.
