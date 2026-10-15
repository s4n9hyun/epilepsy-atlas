import copy
import importlib
import json
from pathlib import Path
import socket
import sys
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT/'tools'), str(ROOT/'paper/analysis')]
from reproduce_snapshot import validate


class SnapshotTests(unittest.TestCase):
    def setUp(self):
        self.results = json.loads((ROOT/'snapshots/manuscript-2026-10-07/results.json').read_text())

    def test_saved_aggregate_arithmetic(self):
        checks = validate(self.results)
        row = next(r for r in checks if r['gene'] == 'SCN8A')
        self.assertEqual((row['DEE'], row['classified_patients'], row['DEE_percent']), (396, 529, 74.9))

    def test_wrong_denominator_is_rejected(self):
        self.results['table1']['SCN8A']['DEE label, n/N (%)'] = '396/805 (74.9)'
        with self.assertRaisesRegex(ValueError, 'denominator mismatch'):
            validate(self.results)

    def test_missing_treatment_subgroup_is_rejected(self):
        self.results['scb']['SCN8A other/not specified']['Benefit'] = '54 (81.8%)'
        with self.assertRaisesRegex(ValueError, 'outcome partition'):
            validate(self.results)

    def test_family_count_uses_one_row_per_patient(self):
        import pandas as pd
        from common import family_ids, load_patients
        rows = pd.DataFrame([
            dict(gene='SCN1A', canonical_patient_id='SYNTHETIC1', pmid='SYNTHETIC', patient_id='II-1', variant_c='c.1A>G', mention_id='M1'),
            dict(gene='SCN1A', canonical_patient_id='SYNTHETIC1', pmid='SYNTHETIC', patient_id='II-1', variant_c='c.2A>G', mention_id='M2')])
        self.assertEqual(family_ids(rows)[0].nunique(), 2)
        self.assertEqual(family_ids(load_patients(rows))[0].nunique(), 1)

    def test_import_does_not_load_data_or_use_network(self):
        with patch('socket.socket.connect', side_effect=AssertionError('Unexpected network access')):
            module = importlib.import_module('run_analysis')
        self.assertFalse(hasattr(module, 'rows'))


if __name__ == '__main__':
    unittest.main()
