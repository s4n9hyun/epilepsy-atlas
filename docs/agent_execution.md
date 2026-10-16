# Agent execution environment, concurrency, time and cost

This record was reconstructed on 2026-10-16 from the retained **project-specific Claude Code JSONL sessions**, the extraction launcher/process logs, task manifests, and the research workflow log. It describes the historical database work, not the GitHub packaging session. Machine-readable evidence is in [runtime/observed_runtime.json](runtime/observed_runtime.json). Raw conversations, article text, patient data, tool arguments and credentials are not uploaded.

## Two separate execution systems

| Stage | Execution environment and model | Concurrency | Evidence and qualification |
|---|---|---|---|
| Initial six-gene extraction | Python batch extraction through the OpenAI Responses API; `gpt-5.4` recorded in the historical preparation notes | Not recovered for the complete historical phase | Do not infer concurrency from a later launch script |
| Seven-gene extraction, 2026-09-28 | Bash launcher + Python processes + OpenAI Python SDK, Responses API; `gpt-6-sol` | **23 initial worker processes**: SMC1A 1, FOXG1 2, GABRA1 4, KCNT1 4, PCDH19 4, STXBP1 4, SCN8A 4 | 23 matching initial worker logs and the workflow log agree. Workers are not autonomous review agents, and their count is not a measurement of simultaneous API requests. Remaining papers were later repartitioned and restarted. |
| Six-gene source-review batch, 2026-09-27 | **Claude Code CLI 2.1.283**, `claude-opus-5-5`, with parent-session orchestration and the `Agent` tool | **8 agents per dispatch batch**, as documented in the review run instructions | 54 planned tasks covering 415 papers/1,966 input rows; this was the remaining-paper review batch, not the whole six-gene project |
| Seven-gene source-review batch, 2026-09-28–29 | **Claude Code CLI 2.1.283**, `claude-opus-5-5`, using the `Agent` tool | **8 agents per dispatch batch** | 91 planned tasks covering 1,359 papers: A=493 retained-patient papers, B=362 papers whose extracted rows were excluded, C=504 papers with no extracted patients |
| Later curation/manuscript work | The retained main session includes multiple Claude Code versions; e.g. 2.1.290 on 2026-10-07 | Phase-specific counts not established | Do not apply the eight-agent setting to every subsequent task |

The CLI identity is established by session metadata (`entrypoint: cli`, version fields, model fields and `Agent` calls) in the local Claude project session. This is stronger evidence than the manuscript's generic phrase “Claude-based agents.” The named review subagent traces themselves record `claude-opus-5-5` and version 2.1.283.

The extraction logs show a macOS-style Python 3.9 user-library environment and an urllib3 warning reporting LibreSSL 2.8.3. An exact historical macOS release, machine CPU/RAM specification, OpenAI SDK version and full dependency lock were not recovered. The repository's Python 3.13 environment is the **tested sharing/analysis environment**, not a claim about the original extraction runtime. Neither model was hosted locally; no local GPU requirement is established by the records.

## Tools and orchestration

The retained main session records `Bash`, `Agent`, `Read`, `Write`, `Edit`, `SendMessage`, `TaskStop`, `ScheduleWakeup`, `Monitor`, and `CronCreate`/`CronList`/`CronDelete`, plus search/fetch tools. Deduplicated tool-call counts are in the evidence JSON; they cover mixed project work, not just source review.

The source-review prompts instruct agents to inspect local main articles and supplements with `pdftotext -layout`, Python, ZIP extraction, pandas and openpyxl. They explicitly prohibit calling another extraction/LLM API from the review task. The parent assigns task IDs/PMIDs; workers write per-paper result CSVs. Tasks that already have results are skipped on resume. The run instructions describe eight-agent batches, periodic session checks, and re-queuing tasks recorded as running for more than three hours. The seven-gene dispatcher uses `tasks.csv`, `running.csv` and `next.py`; it is a queue bookkeeping helper, not a standalone agent runtime.

Source-review concurrency of eight is the **documented dispatch policy**, not a reconstructed continuous utilization trace. Seventy matching six-gene and 101 timestamped seven-gene subagent traces were retained. These counts include attempts/restarts and cannot be substituted for the 54 and 91 planned task counts. Historical prompt paths were later reorganized: `20260927_full_source_check` maps to `runs/six_gene_source_check/full_source_check`, and `20260928_seven_gene_source_check` maps to `runs/seven_genes/source_check`.

For extraction, the saved launcher partitions each gene's paper list into up to four initial chunks and launches Python workers in the background. Requests within each worker are handled by the batch extraction script. The recorded setup uses a 120-second SQLite lock wait and four OpenAI-client retries; later slow-paper retries used 600–900-second timeouts. The current source does not explicitly set temperature, seed or reasoning effort for the main extraction request. Those unspecified settings must not be represented as fixed historical parameters.

## Time: what can actually be reported

| Measurement | Recorded or reconstructed value | Meaning |
|---|---|---|
| Seven-gene extraction launch to recorded completion | **2026-09-28 11:19–15:35**, about **4 h 16 min** | Workflow-log wall-clock interval; includes credit/quota interruptions and restarts. Log timestamps have no UTC offset, so their historical timezone is not independently established. |
| Extraction interruption checkpoints | **12:22** and **14:47** | Credit/quota exhaustion; the first interruption had 1,185/1,539 gene–article entries processed according to the workflow log |
| Six-gene review matched event window | **2026-09-26 18:05:12.934 UTC → 2026-09-27 14:35:29.277 UTC**, **20.50 h** | First-to-last event in matched review subagent traces; includes idle time, failed/restarted attempts and scheduling gaps |
| Seven-gene review matched event window | **2026-09-28 11:08:30.411 UTC → 2026-09-29 01:13:23.289 UTC**, **14.08 h** | Same trace-window definition; not time to finish every later correction or the final database |

The historical date labels use local run names; JSONL timestamps above are explicitly UTC. The extraction and review measurements have different boundaries and must not be added to claim total human time or end-to-end project duration. Review-trace spans are not billable compute time. No controlled measurement of clinician labor saved is available.

## Cost: recorded estimate versus actual payment

The retained main Claude Code session contains **13 cumulative `cost-state` records**. The final retained counter reports **USD 2,428.1402894**, rounded to **USD 2,428.14**. This is a **Claude Code-reported cumulative usage-cost estimate for the mixed project session**, not an invoice, not a source-review-only cost, and not the total cost of the study.

The final counter's model buckets are:

| Logged model bucket | Reported cost, USD |
|---|---:|
| `claude-opus-5-5[1m]` | 2,058.2706224 |
| `claude-opus-5-5` | 355.5882140 |
| `claude-opus-5` | 13.3953820 |
| `claude-haiku-4-5-20251001` | 0.8860710 |

These bucket labels are reproduced as logged; the bracket suffix is not treated as a separate confirmed review model. The named six-/seven-gene source-review traces record Opus 5.5. Other buckets reflect broader session activity whose exact role was not independently established.

The final counter also reports 914,972,225 ms total duration (254.16 h), 252,118,571 ms API duration (70.03 h), and 36,121,496 ms tool duration (10.03 h). These are client counters with potentially overlapping/accumulated scopes; they are **not** independent wall-clock measurements or additive measures of labor. The cost records lack their own timestamps. The evidence JSON retains the preceding event timestamp as a location marker only, not as a billing-period endpoint. Session events span 2026-09-23–2026-10-12; the counter must not be extrapolated to earlier extraction or later work.

**Do not sum the 13 checkpoints**: they are cumulative. The `hasUnknownModelCost` flag is false in the final record, but that does not establish invoice accuracy, subscription payments, completeness of subagent accounting, or complete project coverage. We do not add message-level estimates to the session counter, because that could double-count usage.

For **OpenAI extraction**, retained code/logs confirm separate API calls and quota exhaustion, but a reconciled per-run billed amount was not recovered. The inspected extraction implementation does not persist `response.usage`; a logged planning estimate of input tokens is not measured usage. Consequently, OpenAI extraction cost, total study cost, and monetary savings are **not reported as numeric totals**. No historical costs were back-calculated using today's model prices.

## Evidence and audit reproduction

- `runtime/observed_runtime.json`: metadata-only aggregate, source hashes, trace windows, CLI versions, tools and cumulative cost checkpoints.
- `runtime/workflow_evidence.json`: bounded facts from historical run instructions and the extraction workflow log; no patient data or raw conversations.
- `tools/audit_agent_runtime.py`: rebuilds the metadata summary locally. It needs original project session logs, which are not distributed.

```bash
python tools/audit_agent_runtime.py \
  --claude-project /path/to/claude/project-session-directory \
  --workspace /path/to/original/research-workspace \
  --out /tmp/observed_runtime.json
```

The collector deduplicates assistant messages/tool calls by their IDs, selects the long-lived project session, and matches review traces by their historical task-folder names. It emits metadata only. It does not run agents, call APIs, infer missing bills, or export conversation contents. Cost fields are client-provided estimates whose pricing/accounting implementation was not independently audited.
