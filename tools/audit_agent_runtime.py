"""Summarize local Claude project metadata without exporting prompts or clinical rows.

Only metadata is emitted. Raw session logs and extraction logs stay local.
Usage: python tools/audit_agent_runtime.py --claude-project PATH --workspace PATH --out FILE
"""
import argparse
from collections import Counter, defaultdict
from datetime import datetime
import hashlib
import json
from pathlib import Path
import re
from statistics import median


def events(path):
    with path.open() as stream:
        for line in stream:
            try:
                yield json.loads(line)
            except json.JSONDecodeError:
                continue


def stamp(s):
    return datetime.fromisoformat(s.replace('Z', '+00:00'))


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def summarize(project, workspace):
    # Select the long-lived project session by real assistant event count, not a hardcoded ID.
    candidates = list(project.glob('*.jsonl'))
    main = max(candidates, key=lambda p: sum(r.get('type') == 'assistant' for r in events(p)))
    versions = defaultdict(set)
    tools = Counter()
    models = Counter()
    entrypoints = Counter()
    seen_messages = set()
    seen_tools = set()
    checkpoints = []
    times = []
    previous_timestamp = None
    for row in events(main):
        ts = row.get('timestamp')
        if ts:
            times.append(ts)
            previous_timestamp = ts
            if row.get('version'):
                versions[ts[:10]].add(row['version'])
        if row.get('entrypoint'):
            entrypoints[row['entrypoint']] += 1
        message = row.get('message', {})
        if row.get('type') == 'assistant' and isinstance(message, dict):
            mid = message.get('id')
            if mid and mid not in seen_messages:
                seen_messages.add(mid)
                models[message.get('model', 'unknown')] += 1
            content = message.get('content', [])
            for block in content if isinstance(content, list) else []:
                if isinstance(block, dict) and block.get('type') == 'tool_use':
                    tid = block.get('id')
                    if tid and tid not in seen_tools:
                        seen_tools.add(tid)
                        tools[block['name']] += 1
        if row.get('type') == 'cost-state':
            checkpoints.append({
                'preceding_event_timestamp_utc': previous_timestamp,
                'reported_cost_usd': row.get('totalCostUSD'),
                'reported_total_duration_ms': row.get('totalDuration'),
                'reported_api_duration_ms': row.get('totalAPIDuration'),
                'reported_tool_duration_ms': row.get('totalToolDuration'),
                'has_unknown_model_cost': row.get('hasUnknownModelCost'),
                'model_usage': row.get('modelUsage')})
    groups = defaultdict(list)
    for path in project.rglob('agent-*.jsonl'):
        prompt = None
        trace_times = []
        trace_models = set()
        trace_versions = set()
        for row in events(path):
            if row.get('timestamp'):
                trace_times.append(row['timestamp'])
            if row.get('version'):
                trace_versions.add(row['version'])
            message = row.get('message', {})
            if isinstance(message, dict):
                if row.get('type') == 'user' and prompt is None:
                    prompt = str(message.get('content', ''))
                if message.get('model') not in (None, '<synthetic>'):
                    trace_models.add(message['model'])
        phase = None
        for name, marker in [('six_gene_source_review', '20260927_full_source_check'),
                             ('seven_gene_source_review', '20260928_seven_gene_source_check')]:
            if marker in (prompt or ''):
                phase = name
        if phase and trace_times:
            first, last = min(trace_times), max(trace_times)
            groups[phase].append({'first': first, 'last': last,
                'span_seconds': (stamp(last)-stamp(first)).total_seconds(),
                'models': sorted(trace_models), 'versions': sorted(trace_versions)})
    phases = {}
    for name, traces in groups.items():
        first = min(t['first'] for t in traces)
        last = max(t['last'] for t in traces)
        phases[name] = {'matched_trace_count': len(traces), 'first_event_utc': first, 'last_event_utc': last,
            'calendar_span_hours': (stamp(last)-stamp(first)).total_seconds()/3600,
            'median_trace_span_minutes': median(t['span_seconds'] for t in traces)/60,
            'models': sorted({m for t in traces for m in t['models']}),
            'versions': sorted({v for t in traces for v in t['versions']}),
            'interpretation': 'Attempts/resumptions, not unique tasks. Event spans include waiting and are not compute or human labor time.'}
    initial = sorted((workspace/'runs/seven_genes/extraction/logs').glob('*_part*_202609281119.log'))
    per_gene = Counter(p.name.split('_part')[0] for p in initial)
    sources = []
    for rel in ['WORKFLOW_LOG.md', 'runs/seven_genes/source_check/README.md',
                'runs/six_gene_source_check/full_source_check/README.md',
                'runs/seven_genes/extraction/run_snapshot_202609281140.sh',
                'scripts/core/batch_extract_patient_mentions.py']:
        sources.append({'workspace_relative_path': rel, 'sha256': digest(workspace/rel)})
    return {'schema_version': 1, 'scope': 'metadata of one retained project session and matching review subagent traces',
        'main_session_sha256': digest(main), 'main_event_range_utc': [min(times), max(times)],
        'entrypoints': dict(entrypoints), 'versions_by_date': {d: sorted(v) for d,v in sorted(versions.items())},
        'unique_assistant_messages_by_model': dict(models), 'unique_tool_calls_by_name': dict(tools),
        'source_review_phases': phases, 'cost_checkpoints': checkpoints,
        'cost_interpretation': 'Claude Code reported cumulative estimates, not invoices. Checkpoints must not be summed. Scope includes mixed project work; excludes separately billed OpenAI extraction. Final retained counter is not total project cost.',
        'extraction_initial_process_logs': {'count': len(initial), 'per_gene': dict(per_gene)},
        'sources': sources,
        'privacy': 'No message text, prompts, tool arguments/results, patient records, credentials, or account identifiers exported.'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--claude-project', type=Path, required=True)
    parser.add_argument('--workspace', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    result = summarize(args.claude_project, args.workspace)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({'source_review_phases': result['source_review_phases'],
                      'extraction_initial_process_logs': result['extraction_initial_process_logs'],
                      'last_cost_checkpoint': result['cost_checkpoints'][-1]}, indent=2))


if __name__ == '__main__':
    main()
