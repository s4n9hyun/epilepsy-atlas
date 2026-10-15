# Data inputs and reproducibility scope

No patient-level records, source-article PDFs/supplements, clinical feedback, login profiles, API keys or raw extraction logs are tracked here.

The frozen canonical-file manifest is `snapshots/manuscript-2026-10-07/input_manifest.json`. SHA256 hashes and row totals identify the original inputs. The corrected family total is patient-derived, not patient–variant-row-derived.

Full analyses require a separate workspace containing:

- `outputs/current/canonical_patients/<GENE>_canonical_patients/*.csv`: canonical records. Use exactly the manifest version, not an arbitrary newer export. The legacy loader selects the last filename.
- `databases/<GENE>.sqlite`: source mentions and article processing status.
- `databases/reference_cds/` and ClinVar/gnomAD caches: versioned reference inputs; live requests can change results.
- `runs/structured_phenotypes/`: source-read comorbidities, inheritance and adjudicated additions/corrections.
- `runs/review_fixes/`, `runs/seven_genes/source_check/`, `runs/crf_0904/`, and `runs/paper/`: review decisions, reconciliation and validation inputs referenced by the analysis modules.
- `outputs/current/curation_change_log.csv`: curation provenance.
- Verified treatment-label CSVs at `paper/analysis/scb_source_verified.csv` and `paper/analysis/treatment_verified_7genes.csv` in the code checkout (untracked).
- Patient-specific manual rules/ID locks for canonical rebuilding.

These artifacts have not yet been packaged for distribution. No public data download URL is claimed. The current repository supports code inspection and an offline aggregate snapshot check, not a complete independent clinical-data reconstruction. Historical model calls and human/source-review decisions are not deterministic reproducible computations.
