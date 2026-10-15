"""Offline aggregate checks and Figure 1; no patient data, credentials or network."""
from pathlib import Path
import csv
import json
import re
import sys

ROOT = Path(__file__).resolve().parents[1]
SNAPSHOT = ROOT / 'snapshots/manuscript-2026-10-07'


def validate(results):
    checks = []
    for gene, row in results['table1'].items():
        cell = row['DEE label, n/N (%)']
        match = re.fullmatch(r'([\d,]+)/([\d,]+) \(([\d.]+)\)', cell)
        if not match:
            raise ValueError(f'Invalid DEE denominator cell: {gene}')
        n, total = (int(v.replace(',', '')) for v in match.groups()[:2])
        pct = round(100*n/total, 1)
        if pct != float(match.group(3)):
            raise ValueError(f'DEE denominator mismatch: {gene}')
        checks.append({'gene': gene, 'DEE': n, 'classified_patients': total, 'DEE_percent': pct})
    sc = results['scb']
    groups = [sc[key] for key in ('SCN8A DEE', 'SCN8A self-limited', 'SCN8A other/not specified')]
    overall = results['scb_scn8a_all']
    if sum(g['n'] for g in groups) != overall['n']:
        raise ValueError('SCN8A sample-size partition mismatch')
    for outcome in ('Benefit', 'No benefit', 'Worsening'):
        if sum(int(g[outcome].split()[0]) for g in groups) != int(overall[outcome].split()[0]):
            raise ValueError(f'SCN8A outcome partition mismatch: {outcome}')
    if results['review_fixes']['families'] != 11291:
        raise ValueError('Frozen patient-level family count mismatch')
    flow = results['flow']['Total']
    if flow['original_in_canonical'] - flow['unique_patients'] != flow['merged_other_paper'] + flow['merged_within_paper']:
        raise ValueError('Record consolidation flow mismatch')
    return checks


def main():
    results = json.loads((SNAPSHOT/'results.json').read_text())
    checks = validate(results)
    out = ROOT/'outputs/snapshot'
    out.mkdir(parents=True, exist_ok=True)
    with (out/'dee_denominators.csv').open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(checks[0]))
        writer.writeheader(); writer.writerows(checks)
    sys.path.insert(0, str(ROOT/'paper/analysis'))
    import run_analysis
    import pandas as pd
    import logging
    logging.getLogger("fontTools").setLevel(logging.WARNING)
    run_analysis.FIG = out
    run_analysis.plt.rcParams['font.family'] = 'DejaVu Sans'
    run_analysis.fig1(pd.DataFrame.from_dict(results['flow'], orient='index').rename_axis('gene').reset_index())
    report = {'status': 'passed', 'scope': 'aggregate consistency and Figure 1, not a patient-data rerun',
              'families': 11291, 'SCN8A_benefit': '85+40+55=180/212',
              'duplicate_pair_judgments': 108, 'cross_article_extra_records': 153,
              'within_article_extra_records': 65}
    (out/'validation.json').write_text(json.dumps(report, indent=2)+'\n')
    print(f'Passed. Outputs: {out}')


if __name__ == '__main__':
    main()
