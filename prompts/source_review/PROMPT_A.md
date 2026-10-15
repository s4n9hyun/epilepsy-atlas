You are verifying extracted patient data against the original articles for an epilepsy gene–patient dataset.
Do NOT call any OpenAI/LLM API; read the files yourself (pdftotext -layout, python, unzip, openpyxl, pandas).
Do not modify any file except the result files named below. Put ALL temporary files (text extracts, helper scripts) only in <WORKSPACE>/runs/seven_genes/source_check/tmp/TASK_ID/ (create it); other agents run in parallel, so never write to shared scratch folders.

Your task: TASK_ID. Papers (PMIDs): PMID_LIST
For each PMID:
- Our rows for that paper: <WORKSPACE>/runs/seven_genes/source_check/input/<PMID>.csv
  (one row per patient–variant; columns: gene, canonical_patient_id, mention_id, pmid, patient_id (in-paper ID), location_in_paper,
  variant_c, variant_p, variant_type, sex, seizure_presence, age_of_onset, last_followup_age, seizure_types,
  phenotype_review_category, response_to_asm, targeted_treatment_response, mention_source_type, cited_source)
- Article files: <WORKSPACE>/article_files/<PMID>/ (main.pdf or main.fulltext.html, plus supplementary files; unzip archives if needed; ignore the _excluded/ subfolder, which holds files known to be wrong)
- Skip a PMID whose result file already exists (the job may be resuming).

For every row, find that individual in the article and check these fields against the source for THAT individual:
variant_c, variant_p, sex, age_of_onset, last_followup_age, seizure_presence, phenotype_review_category, targeted_treatment_response, mention_source_type.
Rules:
- variant_c is on the MANE transcript; a paper using legacy numbering for the same variant is a match. variant_p may be translated from the cDNA when the paper gives none — that is fine if correct.
- A value supported only by a cohort-level statement (e.g. "5 of 13 patients worsened on carbamazepine") is NOT an individual value: verdict cohort_level.
- "Not specified" is correct when the paper does not report that field for this individual.
- phenotype_review_category is a gene-specific category (allowed values and hints per gene: <WORKSPACE>/config/gene_phenotype_review_guidance.csv); judge whether it fits what the paper says about this patient (syndrome names must be stated or clearly described, not inferred).
- mention_source_type: original_patient = first report of this patient's genotype; if the paper says the patient was previously reported, give the earlier reference.
- Also note patients reported individually in the paper (with a variant in the gene and clinical data) that are missing from our rows, and our rows that do not exist in the paper (e.g. animal models, cell lines without clinical data, cohort summaries, other genes).

- First confirm the main file really is this article (title/authors). If it is a different article or only an abstract, write one line field="(paper)", verdict wrong_file or abstract_only, evidence = what the file is, and still check what you can.
- last_followup_age = the latest age at which the patient is described (current age, last evaluation, re-examination, enrolment/examination if nothing later, age at death).

Write one CSV per PMID: <WORKSPACE>/runs/seven_genes/source_check/results/<PMID>.csv with columns exactly:
pmid,gene,canonical_patient_id,mention_id,patient_id,field,current_value,verdict,source_value,evidence
- For each of our rows write ONE line with field="(patient)" and verdict found | not_in_paper | not_a_patient.
- Then one line for each field that is NOT correct, verdict one of: wrong | not_in_source | cohort_level | partial ; source_value = the value the source supports for this individual (same short style as our value); evidence = short verbatim quote or exact table/figure location.
- For patients missing from our rows: field="(missing_patient)", patient_id = in-paper ID, source_value = variant; sex; onset; phenotype, evidence = location.
- If the article cannot be read at all: one line field="(paper)", verdict unreadable, evidence = reason.
Write each PMID's file as soon as you finish that paper (do not wait until the end).
Finish with a short summary: papers done, rows checked, number of wrong/not_in_source/cohort_level/partial fields, missing patients, not_in_paper rows.

MANE Select transcripts (use these for cDNA/protein checks): PCDH19 NM_001184880.2, SCN8A NM_001330260.2, STXBP1 NM_001032221.6, GABRA1 NM_001127644.2, KCNT1 NM_020822.3, FOXG1 NM_005249.5, SMC1A NM_006306.4. Our stored cDNA for PCDH19/SCN8A/STXBP1/GABRA1 may be on an older transcript (NM_001105243.2, NM_014191.4, NM_003165.6, NM_000806.5); a difference caused only by that transcript choice is not an extraction error, but note the MANE value in source_value.
