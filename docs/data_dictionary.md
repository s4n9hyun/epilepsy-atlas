# 데이터 사전

DB 테이블, CSV 컬럼, 변이 표기 기준을 정리한 문서입니다.

## 핵심 개념

### Patient Mention

`patient_mentions`의 한 row는 “한 논문 안에서 언급된 한 환자-변이 조합”입니다.

예를 들어 한 논문이 환자 7명을 보고하면, 원칙적으로 7개 row가 생성됩니다. 같은 환자가 다른 review 논문에 다시 등장하면 그 mention도 별도 row로 저장합니다.

이렇게 하는 이유는 원문 출처를 잃지 않기 위해서입니다. 같은 환자인지 여부는 별도 deduplication 단계에서 판단합니다.

### Mention Source Type

논문에 등장하는 환자가 모두 그 논문에서 새로 보고된 환자는 아닙니다.

`mention_source_type`은 다음 값을 가집니다.

| 값 | 의미 |
| --- | --- |
| `original_patient` | 해당 논문에서 새로 보고한 환자 |
| `cited_patient` | 기존 문헌에서 가져오거나 다시 요약한 환자 |
| `aggregate_only` | 개별 환자 정보 없이 집계만 있는 경우 |
| `variant_only` | 환자별 임상정보 없이 variant list만 있는 경우 |
| `previously_reported_patient` | 기존 보고 환자를 다시 제시한 경우 |
| `reanalyzed_previous_patient` | 기존 환자를 재분석한 경우 |
| `manual_exclude` | 원문 확인 후 canonical에서 제외하기로 수동 표시한 경우 |
| `unclear` | 새 환자인지 인용 환자인지 불명확한 경우 |

모델 prompt에서 직접 요청하는 기본 값은 `original_patient`, `cited_patient`, `aggregate_only`, `variant_only`, `previously_reported_patient`, `reanalyzed_previous_patient`, `unclear`입니다. `manual_exclude`는 추출 이후 수동 QC 과정에서 canonical 제외/정리 목적으로 사용됩니다.

Canonical patient-variant output은 `original_patient` 중 strict source filter와 환자별 임상정보 조건을 통과한 row를 중심으로 만듭니다. 제외된 row는 삭제하지 않고 `{GENE}_excluded_mentions_YYYYMMDD_HHMMSS.csv`와 `all_genes_excluded_mentions_YYYYMMDD_HHMMSS.csv`에 보존합니다.

### External Dataset Link

논문 supplementary package가 아니라 본문에 Dryad, Figshare, Zenodo, institutional repository 같은 외부 dataset 링크가 있는 경우가 있습니다. 이런 링크는 DB의 `external_dataset_links` 테이블에 저장하고, 검토용 CSV export로도 관리합니다.

이 링크는 publisher supplementary material URL과 구분합니다. 저장 대상은 환자별 genotype, variant, phenotype 정보를 담을 가능성이 있는 dataset 링크입니다. 단순 홈페이지, ClinVar/LOVD variant page, code-only repository, RNA-seq/expression-only dataset은 제외합니다.

대표 예시는 PMID 30541864의 Dryad dataset입니다. 이 논문의 환자-변이 정보는 일반 supplementary file이 아니라 논문 본문에 연결된 Dryad Excel dataset으로 제공됩니다.

### Canonical Patient ID Lock

재추출이나 추가 논문 반영 후에도 기존 리뷰 결과의 CP 번호가 흔들리지 않도록 `{GENE}_canonical_id_locks.csv`를 사용합니다. 기존 mention이 이미 `CP0000123` 같은 canonical ID를 받았다면, rebuild 시 같은 ID를 유지합니다.

현재 `outputs/current/canonical_patients/`는 교수님 검토용 최종 canonical CSV 6개만 남기도록 정리되어 있습니다. CP lock과 dedup summary 같은 보조 파일은 `outputs/auxiliary/`에 보존합니다. Rebuild 전에는 해당 lock 파일을 작업 output directory에 복사하거나, archive의 보조 파일을 기준으로 작업해야 합니다.

### Timestamped CSV

앞으로 새로 생성되는 결과 CSV는 기본적으로 파일명 끝에 `YYYYMMDD_HHMMSS` timestamp를 붙입니다.

예:

```text
SYNGAP1_canonical_patients_20260616_161440.csv
multi_variant_patient_review_20260616_161440.csv
```

`scripts/core/build_final_outputs.py`와 `scripts/core/audit_canonical_qc.py`는 timestamp가 붙은 최신 CSV를 자동으로 찾아 읽습니다. 고정 파일명이 필요한 예외 상황에서는 각 스크립트에 `--no-timestamp`를 붙이면 기존 방식의 파일명을 쓸 수 있습니다.

### Manual Review Configs

수동 검토로 확정된 예외는 코드에 직접 박지 않고 `config/` 아래 CSV로 관리합니다.

| 파일 | 의미 |
| --- | --- |
| `manual_merge_rules.csv` | 서로 다른 mention이 같은 환자임이 원문 근거로 확인된 경우 canonical patient로 병합 |
| `manual_variant_exclusions.csv` | 같은 환자에게 여러 변이가 보고되었지만, causal/P/LP variant가 아닌 secondary variant를 canonical causal output에서 제외 |
| `gene_inheritance_models.csv` | 유전자별 inheritance model과 허용 causal variant 수. 별도 `multi_variant_patient_review_YYYYMMDD_HHMMSS.csv` 생성에 사용 |
| `manual_article_inclusions.csv` | `Gene AND epilepsy` PubMed 검색식에서 빠졌지만 교수님/manual review로 추가해야 하는 논문 목록 |
| `gene_phenotype_review_guidance.csv` | 유전자별 phenotype category, onset/follow-up 단위, 치료 반응 추출 포인트 |
| `crf_phenotype_label_map.csv` | canonical phenotype 라벨을 교수님 CRF 라벨로 바꾸는 매핑. SYNGAP1은 Phenotype + Phenotype_subcategory 2단 구조 |

`manual_variant_exclusions.csv`로 제외된 row는 삭제하지 않습니다. 원본 `patient_mentions`와 `{GENE}_excluded_mentions_YYYYMMDD_HHMMSS.csv`에 보존하고, canonical patient-variant output에서만 제외합니다.

`manual_article_inclusions.csv`에 들어간 논문은 DB의 `pubmed_articles`에도 추가하며, 다운로드/추출 전에는 `source_csv=config/manual_article_inclusions.csv`, `extraction_status=manual_inclusion_pending_download`로 추적합니다.

Inheritance model 기준은 canonical patient CSV 안에 직접 표시하지 않습니다. `scripts/core/audit_canonical_qc.py`가 `config/gene_inheritance_models.csv`를 읽어 별도 QC 파일인 `multi_variant_patient_review_YYYYMMDD_HHMMSS.csv`를 생성합니다.

### Confidence

`confidence`는 추출 row의 신뢰도입니다.

| 값 | 의미 |
| --- | --- |
| `High` | 환자 ID, variant, phenotype이 같은 표/문단에서 명확히 연결됨 |
| `Medium` | 정보가 분산되어 있지만 연결이 비교적 명확함 |
| `Low` | 연결 가능성은 있으나 불확실성이 있음 |


## DB 구조

각 유전자 DB에는 주로 다음 테이블이 있습니다.

### `pubmed_articles`

유전자별 PubMed 논문 목록입니다. 다운로드와 추출 진행 상황을 관리합니다.

주요 컬럼:

| 컬럼 | 의미 |
| --- | --- |
| `pmid` | PubMed ID. 파일 폴더 관리의 기본 ID |
| `pmcid` | PMC ID. 없는 논문도 있을 수 있음 |
| `gene` | 유전자명 |
| `title` | 논문 제목 |
| `downloaded` | 논문 파일 다운로드 여부 |
| `extracted` | 추출 완료 여부 |
| `extraction_status` | 추출 상태 |
| `extracted_at` | 추출 시각 |

### `processed_articles`

논문 단위 처리 로그입니다. 한 논문을 성공적으로 처리했는지, 몇 명의 patient mention이 나왔는지 기록합니다.

주요 컬럼:

| 컬럼 | 의미 |
| --- | --- |
| `article_id` | 처리 기준 ID. 보통 PMID |
| `pmid` | PubMed ID |
| `pmcid` | PMC ID |
| `gene` | 유전자명 |
| `status` | `success` 또는 실패 상태 |
| `patient_count` | 해당 논문에서 추출된 patient mention 수 |
| `files_processed` | 실제 추출에 사용한 파일명 목록 |
| `error` | 실패한 경우 오류 내용 |
| `processed_at` | 처리 시각 |
| `supplement_note` | 본문이 가리키는 환자별 supplement가 첨부 파일에 없을 때 모델이 남긴 설명 (2026-09-23 추가) |

### `article_files`

`article_files/{PMID}/`에 저장된 실제 파일을 DB에서 추적하기 위한 index 테이블입니다.

주요 컬럼:

| 컬럼 | 의미 |
| --- | --- |
| `pmid` | PubMed ID |
| `pmcid` | PMC ID |
| `file_name` | 저장된 파일명 |
| `relative_path` | `article_files/` 기준 상대 경로 |
| `file_type` | PDF, NXML, HTML, XLSX, DOCX 등 파일 유형 |
| `is_supplementary` | supplementary file 여부 |
| `file_size` | 파일 크기 |
| `indexed_at` | 파일 index 시각 |

### `external_dataset_links`

논문 본문에 연결된 외부 dataset 링크입니다. PMC supplementary file과 별도로 관리합니다.

주요 컬럼:

| 컬럼 | 의미 |
| --- | --- |
| `pmid` | PubMed ID |
| `pmcid` | PMC ID |
| `gene` | 유전자명 |
| `repository` | Dryad, Figshare, Zenodo 등 repository 이름 |
| `url` | dataset URL |
| `doi` | dataset DOI |
| `link_text` | 논문 안에서 링크 주변에 보인 텍스트 |
| `evidence_location` | 링크가 나온 본문 위치 |
| `expected_content` | 환자/변이/phenotype 정보가 있을 것으로 기대되는 내용 |
| `access_status` | 접근 가능 여부 |
| `manual_download_required` | 수동 다운로드 필요 여부 |
| `extraction_status` | 추출 반영 상태 |
| `local_file_path` | 다운로드한 외부 dataset 파일 경로 |

### `patient_mentions`

가장 중요한 추출 결과 테이블입니다. 한 row가 한 patient mention입니다.

DB 내부 핵심 컬럼:

| 컬럼 | 의미 |
| --- | --- |
| `variant_c_original` | 논문에서 읽은 원문 variant. 유효 HGVS cDNA뿐 아니라 genomic coordinate, CNV, 비표준 cDNA-like 표기도 보존할 수 있음 |
| `variant_c_gnomad` | 비교/ClinVar용 cDNA variant. 기준 transcript를 붙여 검증한 HGVS 표기 |
| `variant_p` | HGVS protein 표기. 단일 값만 저장 |
| `variant_type` | `sequence_variant_or_small_indel`, `exon_level_cnv`, `large_cnv`, `unclear` 중 하나 |
| `clinvar` | DB 내부 trace용 legacy 컬럼. 예: `listed: VCV000192252.60; AlleleID=190049` (Variation ID는 VCV 번호 192252), 없으면 `not_listed`, 조회 대상이 아니면 `not_queried` |
| `clinvar_pathogenicity` | ClinVar germline classification 중 `Pathogenic`, `Likely pathogenic`, `Pathogenic/Likely pathogenic`만 표시. 그 외 classification은 빈칸 |
| `clinvar_id` | ClinVar Variation ID만 저장 (VCV accession 번호, 예: VCV000192252 → 192252). 2026-09-23 이전 값은 Allele ID였으며 전 유전자 정정함 |
| `patient_id` | 논문 안에서의 환자 ID, case 번호, proband ID |
| `pmid` | 출처 논문의 PMID |

환자 phenotype 컬럼:

| 컬럼 | 의미 |
| --- | --- |
| `sex` | `Male`, `Female`, `Not specified` |
| `korean_patient_source` | `Korean multicenter cohort`, `Korean patient from literature`, `Not Korean / Not specified` |
| `korean_patient_evidence` | 국내 환자 source 판단 근거 |
| `seizure_presence` | seizure 여부 |
| `age_of_onset` | seizure 또는 증상 발병 나이 |
| `last_followup_age` | 마지막 추적 나이 |
| `seizure_types` | seizure type |
| `phenotype_review_category` | 리뷰용 phenotype category |
| `response_to_asm` | ASM 반응 요약 |
| `asm_used` | 사용한 anti-seizure medication |
| `responsive_asm` | 효과가 있었다고 명시된 ASM |
| `targeted_treatment_response` | 특정 치료 반응 |
| `epilepsy_type` | epilepsy type |
| `other_phenotypes` | 발달, 행동, 신경학적 동반 phenotype |

출처와 검토용 컬럼:

| 컬럼 | 의미 |
| --- | --- |
| `mention_source_type` | source type 구분. 자세한 값은 [Mention Source Type](#mention-source-type) 참고 |
| `originality_evidence` | 왜 original 또는 cited로 판단했는지 근거 |
| `cited_source` | 인용 환자인 경우 원 출처 |
| `location_in_paper` | 표 번호, supplement 파일명, 문단 위치 |
| `confidence` | 추출 신뢰도 |
| `used_supplementary` | supplementary 사용 여부 |
| `files_processed` | 분석에 사용한 파일명 목록 |
| `raw_patient_json` | 모델이 반환한 원 JSON. 검토와 재처리용 |
| `extraction_timestamp` | 추출 시각 |

오른쪽 metadata 컬럼:

| 컬럼 | 의미 |
| --- | --- |
| `gnomad_conversion_status` | cDNA 변환/검증 상태 |
| `nm_transcript` | 논문에 명시된 RefSeq transcript |
| `gnomad_transcript` | 비교 기준 transcript |
| `pmcid` | PMC ID |
| `mention_id` | DB 내부 row ID |

2026-09-23 추가 컬럼 (이후 추출분부터 채워짐):

| 컬럼 | 의미 |
| --- | --- |
| `variant_classification_in_paper` | 논문이 붙인 변이 분류: pathogenic / likely pathogenic / VUS / likely benign / benign / polymorphism / not stated |
| `variant_evidence_source` | 변이를 읽은 위치: table / supplementary_table / main_text / figure_only |
| `variant_phase` | 한 환자에 대상 유전자 변이가 둘 이상일 때: cis / trans / unknown_phase (하나면 single) |
| `clinvar_classification` | ClinVar 원 분류 (P/LP 외 VUS·상충·benign 포함) |

## 변이 표기 기준

### cDNA Variant

DB에는 두 종류의 cDNA 관련 variant 컬럼이 있습니다.

| 컬럼 | 의미 |
| --- | --- |
| `variant_c_original` | 논문에서 읽은 원문 variant. 예: `c.2528A>G`, `chr6 g.33406569 CTGTATG>CTG`, `93kb deletion@6P21.32` |
| `variant_c_gnomad` | 비교 기준 transcript를 붙여 검증한 HGVS 표기. 예: `NM_006772.3:c.2528A>G` |

논문마다 사용하는 transcript가 다르면 같은 생물학적 변이라도 cDNA 좌표가 달라질 수 있습니다. 그래서 원문 표기와 비교용 표기를 분리합니다.

`variant_c_original`에는 가능한 경우 유효한 cDNA HGVS를 저장합니다. 다만 논문에 genomic coordinate, VEP prediction, CNV, exon-level deletion, free-text deletion처럼 표준 HGVS cDNA로 정규화되지 않은 genotype이 환자별로 명시되어 있으면 그 원문을 보존합니다. 이 값은 CSV에서 `variant_c_source = raw_text`로 표시됩니다.

`variant_c_gnomad`를 만들 때는 Mutalyzer API를 사용합니다. `c.`로 시작하더라도 공백, OCR 오타, 잘못된 range 표기처럼 느슨한 값이 있을 수 있기 때문에 기준 transcript를 붙여 검증합니다.

검토용 CSV에는 `variant_c_gnomad`와 `variant_c_original`을 따로 노출하지 않고 다음 두 컬럼으로 합쳐 표시합니다.

| CSV 컬럼 | 의미 |
| --- | --- |
| `variant_c_source` | `variant_c`가 어디서 온 값인지 표시. `converted_with_nm_number`, `raw_text`, `not_specified` 중 하나 |
| `variant_c` | `variant_c_gnomad`가 있으면 그 값을 사용하고, 없으면 `variant_c_original` 원문을 사용 |

### 기준 Transcript

비교용 transcript는 gnomAD/ClinVar 대조를 위해 MANE Select RefSeq 기준으로 맞춥니다.

```python
GNOMAD_TRANSCRIPTS = {
    "SCN1A":   "NM_001165963.4",
    "SCN2A":   "NM_001040142.2",
    "CDKL5":   "NM_001323289.2",
    "EEF1A2":  "NM_001958.5",
    "KCNQ2":   "NM_172107.4",
    "SYNGAP1": "NM_006772.3",
}
```

코드 변수명은 `GNOMAD_TRANSCRIPTS`이지만, 실제 의미는 “gnomAD/ClinVar 비교를 위한 MANE Select RefSeq 기준 transcript”입니다.

### `gnomad_conversion_status`

| 상태 | 의미 |
| --- | --- |
| `already_gnomad` | 논문 transcript가 이미 기준 transcript와 같고 Mutalyzer 확인을 통과함 |
| `assumed_gnomad` | 논문에 NM transcript가 없어서 기준 transcript로 가정했고 Mutalyzer 확인을 통과함 |
| `converted` | Mutalyzer가 다른 transcript 표기를 기준 transcript 표기로 변환함 |
| `no_gnomad_found` | Mutalyzer 결과에서 기준 transcript 표기를 찾지 못함 |
| `not_specified` | cDNA variant가 명시되지 않음 |
| `invalid_hgvs` | cDNA HGVS가 아니거나 복구 불가능한 표기 |
| `mutalyzer_rejected` | HGVS처럼 보이지만 Mutalyzer가 최종적으로 거절함 |
| `api_error` | API 호출 실패. 정상 저장 흐름에서는 가능하면 남기지 않음 |

Mutalyzer에서 오류가 나면 스크립트는 해당 row만 GPT로 다시 짧게 QC합니다. 명확한 OCR/형식 오류만 보정하고, 그래도 실패하면 `invalid_hgvs` 또는 `mutalyzer_rejected`로 분류합니다.

### Protein Variant

`variant_p`는 transcript 변환 대상은 아닙니다. 대신 비교하기 쉽도록 HGVS protein 3-letter amino acid 표기로 정리합니다.

예:

```text
p.E122K        -> p.Glu122Lys
p.R687X        -> p.Arg687Ter
p.Arg687*      -> p.Arg687Ter
p.(R687*)      -> p.Arg687Ter
p.K114Sfs*20   -> p.Lys114SerfsTer20
```

정책:

- `p.(...)`의 괄호는 저장하지 않습니다.
- 한 변이에 여러 splice transcript/protein product가 병기된 경우는 단일 protein variant로 보기 어려우므로 `Not specified`로 저장합니다.
- protein consequence가 명시되지 않거나 읽을 수 없으면 `Not specified`로 저장합니다.
- `p.?`는 그대로 둡니다.


## 결과 CSV별 컬럼 설명

#### `{GENE}_patient_mentions_YYYYMMDD_HHMMSS.csv`

DB의 `patient_mentions`를 그대로 검토하기 쉽게 내보낸 CSV입니다. 한 row는 한 논문에서 발견된 한 patient mention입니다. Deduplication 전 원본 mention이므로, 같은 환자가 여러 논문에서 반복 등장할 수 있습니다.

| 컬럼 | 의미 |
| --- | --- |
| `variant_c_source` | `variant_c`의 출처. `converted_with_nm_number`는 내부 `variant_c_gnomad`, `raw_text`는 내부 `variant_c_original`에서 온 값 |
| `variant_c` | 검토용 통합 cDNA/genotype 컬럼. 표준화된 값이 있으면 `variant_c_gnomad`, 없으면 원문 `variant_c_original`을 표시 |
| `variant_p` | HGVS protein 표기. 단일 protein consequence만 저장 |
| `variant_type` | sequence/small indel, exon-level CNV, large CNV 구분 |
| `clinvar_pathogenicity` | ClinVar germline classification 중 `Pathogenic`, `Likely pathogenic`, `Pathogenic/Likely pathogenic`만 표시. 없으면 빈칸 |
| `ClinVar ID` | ClinVar Variation ID만 표시 |
| `patient_id` | 논문 안의 환자 ID, case 번호, proband ID |
| `pmid` | 출처 논문의 PubMed ID |
| `sex` | 성별. `Male`, `Female`, `Not specified` |
| `korean_patient_source` | 국내 환자 source 구분 |
| `korean_patient_evidence` | 국내 환자 source 판단 근거 |
| `seizure_presence` | seizure가 보고되었는지 여부 |
| `age_of_onset` | seizure 또는 증상 발병 나이 |
| `last_followup_age` | 마지막 추적 관찰 나이 |
| `seizure_types` | 논문에 명시된 seizure type |
| `phenotype_review_category` | 유전자별 review 기준에 맞춘 phenotype category |
| `response_to_asm` | anti-seizure medication 반응 요약 |
| `asm_used` | 사용된 anti-seizure medication 목록 |
| `responsive_asm` | 효과가 있었다고 명시된 ASM |
| `targeted_treatment_response` | 특정 치료, 식이, targeted treatment 반응 |
| `epilepsy_type` | 논문에 명시된 epilepsy type 또는 syndrome |
| `other_phenotypes` | 발달지연, 지적장애, 행동 문제, MRI/EEG 등 기타 phenotype |
| `mention_source_type` | source type 구분. 자세한 값은 [Mention Source Type](#mention-source-type) 참고 |
| `originality_evidence` | original/cited 판단 근거. 예: “our patient”, “previously reported” |
| `cited_source` | cited patient인 경우 원 출처 문헌 |
| `location_in_paper` | 정보가 나온 위치. 표 번호, supplement 파일명, page, section 등 |
| `confidence` | 환자-변이-phenotype 연결 신뢰도. `High`, `Medium`, `Low` |
| `used_supplementary` | supplementary file 사용 여부 |
| `files_processed` | 추출에 사용된 파일명 목록 |
| `raw_patient_json` | 모델이 반환한 원본 JSON. 재검토와 디버깅용 |
| `extraction_timestamp` | 추출 시각 |
| `gnomad_conversion_status` | `variant_c_gnomad` 변환/검증 상태 |
| `nm_transcript` | 논문에 명시된 RefSeq transcript |
| `gnomad_transcript` | 비교 기준 transcript |
| `pmcid` | PMC ID. 없을 수 있음 |
| `mention_id` | DB 내부 patient mention row ID |

#### `outputs/current/canonical_patients/{GENE}_canonical_patients/{GENE}_canonical_patients_YYYYMMDD_HHMMSS.csv`

Manual review와 최종 unique patient DB 후보 생성을 위한 CSV입니다. `original_patient` row 중 strict source filter를 통과한 row만 canonical 후보로 포함합니다. 원본 `{GENE}_patient_mentions_YYYYMMDD_HHMMSS.csv`의 모든 컬럼을 보존하고, 앞뒤에 canonical 검토용 컬럼을 추가합니다.

이 파일은 **unique patient 1명당 1행이 아니라 patient-variant 1개 조합당 1행**입니다. 같은 환자에게 여러 변이가 있으면 같은 `canonical_patient_id`가 여러 행에 반복됩니다. 따라서 환자 수는 row count가 아니라 distinct `canonical_patient_id`로 계산해야 합니다.

앞쪽 추가 컬럼:

| 컬럼 | 의미 |
| --- | --- |
| `canonical_patient_id` | dedup 후 canonical patient 후보 ID |
| `gene` | 유전자명 |

원본에서 보존되는 컬럼:

```text
variant_c_source, variant_c, variant_p, patient_id, pmid, sex, seizure_presence,
age_of_onset, last_followup_age, seizure_types, phenotype_review_category,
response_to_asm, asm_used, responsive_asm, targeted_treatment_response,
epilepsy_type, other_phenotypes, mention_source_type, originality_evidence,
cited_source, location_in_paper, confidence, used_supplementary,
files_processed, raw_patient_json, extraction_timestamp,
gnomad_conversion_status, nm_transcript, gnomad_transcript, pmcid, mention_id
```

뒤쪽 추가 컬럼:

| 컬럼 | 의미 |
| --- | --- |
| `supporting_pmids` | 같은 canonical patient 후보를 지지하는 PMID 목록 |
| `supporting_pmcids` | 같은 canonical patient 후보를 지지하는 PMCID 목록 |
| `source_patient_ids` | 원 논문들에서 사용된 patient ID 목록 |
| `mention_count` | 해당 canonical patient 후보에 연결된 original mention 수 |
| `canonical_patient_mention_count` | 해당 canonical patient 후보에 연결된 전체 original mention 수 |
| `canonical_patient_variant_count` | 해당 canonical patient 후보 안의 patient-variant row 수 |
| `dedup_key` | 자동 병합 판단에 사용한 내부 key |
| `match_confidence` | 자동 병합 신뢰도. 현재는 대부분 `singleton` |
| `match_reason` | canonical 후보가 만들어진 이유 또는 병합/비병합 사유 |

#### `outputs/auxiliary/{GENE}/{GENE}_patient_mentions_YYYYMMDD_HHMMSS.csv`

원본 mention 전체에 canonical 관련 컬럼을 추가한 파일입니다. `original_patient`뿐 아니라 `cited_patient`, `aggregate_only`, `variant_only`, `previously_reported_patient`, `reanalyzed_previous_patient`, `manual_exclude`, `unclear`도 모두 보존합니다.

추가 컬럼:

| 컬럼 | 의미 |
| --- | --- |
| `gene` | 유전자명 |
| `canonical_patient_id` | original patient로 canonical 후보에 포함된 경우 연결된 ID |
| `canonical_status` | `included` 또는 `excluded` |
| `exclusion_reason` | canonical 후보에서 제외된 이유 |

#### `outputs/auxiliary/{GENE}/{GENE}_patient_dedup_links_YYYYMMDD_HHMMSS.csv`

canonical patient 후보와 원 mention을 연결하는 표입니다. 어떤 mention이 어떤 canonical patient 후보에 연결되었는지 추적할 때 사용합니다.

| 컬럼 | 의미 |
| --- | --- |
| `canonical_patient_id` | canonical patient 후보 ID |
| `mention_id` | 원 patient mention ID |
| `pmid` | 출처 PMID |
| `pmcid` | 출처 PMCID |
| `patient_id` | 논문 안의 patient ID |
| `mention_source_type` | source type 구분. 자세한 값은 [Mention Source Type](#mention-source-type) 참고 |
| `match_confidence` | 연결 신뢰도 |
| `match_reason` | 연결 이유 |

#### `outputs/auxiliary/{GENE}/{GENE}_excluded_mentions_YYYYMMDD_HHMMSS.csv`

canonical patient 후보에는 포함하지 않은 mention입니다. 주로 `cited_patient`, `aggregate_only`, `variant_only`, `previously_reported_patient`, `reanalyzed_previous_patient`, `manual_exclude`, `unclear` 또는 strict source filter에서 제외된 row가 들어갑니다. 삭제된 것이 아니라 별도 파일에 보존됩니다.

#### `outputs/auxiliary/{GENE}/dedup_summary.json`

Deduplication 결과 요약입니다. 전체 mention 수, canonical 후보 수, 제외 사유별 count, 자동 병합 cluster 수, 적용한 dedup 정책을 기록합니다.

#### `outputs/auxiliary/{GENE}/{GENE}_canonical_id_locks.csv`

기존 canonical patient ID를 유지하기 위한 lock 파일입니다. 추가 추출이나 dedup rule 변경 후 rebuild하더라도 이미 검토된 mention은 기존 `CP...` 번호를 유지합니다.
