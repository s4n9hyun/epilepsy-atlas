# Epilepsy atlas

Code for literature-based extraction, curation and analysis of patients with variants in 13 epilepsy-associated genes.

This is a **private preparation repository**, not the final publication release. It contains the extraction implementation (including its embedded prompts and output schema), curation scripts, source-review prompt templates, and manuscript analysis code. Original operational folder names are retained so imports and relative paths remain understandable.

## Quick start: no credentials or clinical data needed

Use Python 3.13:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python -m unittest discover -s tests -v
python tools/reproduce_snapshot.py
```

The last command validates the frozen aggregate results and regenerates Figure 1 and a denominator/treatment consistency report under `outputs/snapshot/`. It does **not** call an LLM or re-estimate results from patient records. It is an offline check of the archived analysis outputs.

## Repository map

| Location | Contents |
|---|---|
| `scripts/core/` | Extraction, embedded prompts/schema, normalization, deduplication, QC, SQLite utilities |
| `scripts/curation/` | Correction, canonical rebuilding, transcript checks and review utilities |
| `prompts/source_review/` | Actual source-review templates; replace `<WORKSPACE>`, task and PMID placeholders before use |
| `config/` | Gene phenotype and inheritance guidance; no patient-specific correction rules |
| `paper/analysis/` | Statistical analyses, treatment coding, figure/table and manuscript rendering code |
| `snapshots/manuscript-2026-10-07/` | Aggregate results, selected tables and hashed canonical-input manifest |
| `docs/` | Data dictionary, historical curation policy and source-code provenance |
| `tests/` | Regression tests for denominators, outcome partitions, family counting and data-free imports |

## Historical agent execution

See [Agent execution environment, concurrency, time and cost](docs/agent_execution.md) for the recovered Claude Code CLI/model versions, eight-agent source-review batches, 23-process extraction launch, observed event windows and recorded cost counters. It distinguishes client cost estimates from actual bills and documents unrecovered costs. Supporting metadata and an audit script are included; raw conversations remain local.

## Frozen manuscript versus current code

The manuscript uses the **2026-10-07** data snapshot: 11,982 canonical patients, 12,041 patient–variant rows, 1,437 articles and 11,291 patient-derived families. Aggregate results include the **2026-10-16** consistency corrections. The source code was copied from the working project on 2026-10-16 and may contain later curation utilities; this is not a claim that the current extraction model or every script exactly recreates the historical extraction run.

`input_manifest.json` records filenames and SHA256 hashes for the 13 exact canonical inputs. Patient records are not distributed in this repository. The hash manifest identifies inputs; it does not supply the SQLite databases, source-review decisions or external-reference caches required for a complete rerun. See `data/README.md`.

## Extraction and full analysis

The actual extraction prompt and JSON schema are embedded in `scripts/core/batch_extract_patient_mentions.py`; gene-specific instructions are in `config/gene_phenotype_review_guidance.csv`.

```bash
python scripts/core/batch_extract_patient_mentions.py --help
python scripts/curation/rebuild_canonical.py --help
```

For fresh extraction, install the OpenAI SDK separately (`python -m pip install openai`) and set `OPENAI_API_KEY` in the environment or a local `.env`. Extraction uploads supplied documents to the configured API and incurs API costs; no extraction is triggered by setup, tests or the offline snapshot command. Choose the model explicitly for a new run and record the SDK/model versions. The SDK is deliberately not pinned to an invented historical version.

Full analysis requires the data layout documented in `data/README.md`. Set `EPILEPSY_WORKSPACE` to that data workspace for analysis. Repository source imports remain local. Some operational curation utilities expect inputs under the repository root; consult their CLI help before using them. Do not run mutating curation tools against your only copy of a dataset.

Full paper rendering also needs the manuscript source, journal template and pandoc/tectonic; those assets are not included here. Complete raw-data-to-manuscript reproduction has **not** been validated from this code-only repository.

## Publication preparation

Before a public release: select a software license with the authors, provide the distributable data/review artifacts (or a documented access route), validate a clean full-data run, finalize citation metadata, and tag the exact manuscript version. A GitHub release can then be archived with Zenodo for a version-specific DOI. No publication DOI or final `v1.0.0` release has been created.

Subsequent changes can be committed with `git add <changed-files>`, `git commit -m "Describe the change"`, and `git push`. Edit this repository for shareable code changes; the original research workspace is a separate copy and is not automatically synchronized.
