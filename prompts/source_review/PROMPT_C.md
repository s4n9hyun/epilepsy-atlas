You are screening papers where an LLM found NO patients for an epilepsy gene dataset, to catch missed patients.
Do NOT call any OpenAI/LLM API; read the files yourself (pdftotext -layout, python, unzip, openpyxl, pandas).
Do not modify any file except your result files. Put ALL temporary files only in <WORKSPACE>/runs/seven_genes/source_check/tmp/TASK_ID/ (create it); other agents run in parallel.
Article files: <WORKSPACE>/article_files/<PMID>/ (top-level files only; ignore the _excluded/ subfolder, which holds wrong files). Skip a PMID whose result file already exists.
Genes in scope: SMC1A, FOXG1, GABRA1, KCNT1, PCDH19, STXBP1, SCN8A. Our whole current dataset (to see whether a patient is already present from another paper): <WORKSPACE>/outputs/current/canonical_patients/<GENE>_canonical_patients/*.csv
Inclusion rule for a patient: individually described; carries a genotyped variant in the gene; at least one individual clinical field; not an unaffected carrier; not a multigene deletion covering another epilepsy gene; this paper is the FIRST report of that person's genotype.

Your task: TASK_ID. Papers: PMID_LIST (the gene(s) searched for each PMID are in <WORKSPACE>/runs/seven_genes/source_check/C_papers.csv).
For each PMID read the article and supplements quickly and decide: no_individual_patients (review, animal/cell study, aggregate statistics only, other gene) | only_previously_reported (individual patients but all cited/earlier reports) | missing_patients (individuals meeting the inclusion rule) | wrong_file / abstract_only / unreadable.
Write <WORKSPACE>/runs/seven_genes/source_check/results/C_<PMID>.csv with columns exactly:
pmid,gene,verdict,article_type,patient_id,variant,sex,onset,followup_age,phenotype,cited_pmid_or_ref,evidence
- one line per paper with the verdict (patient fields empty), plus for missing_patients one line per missing individual (verdict missing_patient) with the source values and exact location.
Write each PMID file as soon as it is done. Finish with a short summary.
MANE Select transcripts (use these for cDNA/protein checks): PCDH19 NM_001184880.2, SCN8A NM_001330260.2, STXBP1 NM_001032221.6, GABRA1 NM_001127644.2, KCNT1 NM_020822.3, FOXG1 NM_005249.5, SMC1A NM_006306.4. Our stored cDNA for PCDH19/SCN8A/STXBP1/GABRA1 may be on an older transcript (NM_001105243.2, NM_014191.4, NM_003165.6, NM_000806.5); a difference caused only by that transcript choice is not an extraction error, but note the MANE value in source_value.

Policy: a statement that applies to ALL carriers of a gene in the paper (e.g. "all 3 KCNT1 patients presented with EIMFS") may be used as each of those patients' individual value; a statement about only some of a group (e.g. "5 of 13") may not.
