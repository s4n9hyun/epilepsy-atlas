"""Protect privacy and avoid double-counting cumulative/streamed runtime records."""
import json
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'tools'))
from audit_agent_runtime import summarize


class RuntimeAuditTests(unittest.TestCase):
    def test_stream_duplicates_cost_checkpoints_and_private_content(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp); project = root/'project'; workspace = root/'workspace'
            project.mkdir(); workspace.mkdir()
            for name in ['WORKFLOW_LOG.md', 'runs/seven_genes/source_check/README.md',
                         'runs/six_gene_source_check/full_source_check/README.md',
                         'runs/seven_genes/extraction/run_snapshot_202609281140.sh',
                         'scripts/core/batch_extract_patient_mentions.py']:
                path = workspace/name; path.parent.mkdir(parents=True, exist_ok=True); path.write_text('fixture')
            logs = workspace/'runs/seven_genes/extraction/logs'; logs.mkdir()
            (logs/'SMC1A_part0_202609281119.log').write_text('fixture')
            event = {'type': 'assistant', 'timestamp': '2026-09-28T01:00:00Z', 'version': 'test',
                     'message': {'id': 'same-message', 'model': 'test-model', 'content': [
                         {'type': 'tool_use', 'id': 'same-tool', 'name': 'Bash', 'input': {'command': 'PRIVATE_CONTENT'}}]}}
            records = [event, event, {'type': 'cost-state', 'totalCostUSD': 1.0},
                       {'type': 'cost-state', 'totalCostUSD': 2.0}]
            (project/'main.jsonl').write_text(''.join(json.dumps(r)+'\n' for r in records))
            result = summarize(project, workspace)
            self.assertEqual(result['unique_assistant_messages_by_model'], {'test-model': 1})
            self.assertEqual(result['unique_tool_calls_by_name'], {'Bash': 1})
            self.assertEqual([r['reported_cost_usd'] for r in result['cost_checkpoints']], [1.0, 2.0])
            self.assertEqual(result['extraction_initial_process_logs']['count'], 1)
            self.assertNotIn('PRIVATE_CONTENT', json.dumps(result))
            self.assertNotIn(str(root), json.dumps(result))


if __name__ == '__main__':
    unittest.main()
