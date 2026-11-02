"""Independent inventory/provenance audit. Does not run or alter any fit."""
from pathlib import Path
import argparse
import datetime
import hashlib
import json
import math
import os

P = Path(__file__).resolve().parents[2]
R = P / 'runs/20261031_shared_sr_cr_scan'
MASSES = [f'{n / 10:.1f}' for n in range(60, 121)]


def read(path):
    return json.loads(path.read_text())


def main(final=False):
    issues, missing, numerical = [], {}, {}
    targets = {
        'observed': {f'observed/m{m}/results.json' for m in MASSES},
        'observed_outcomes': {f'observed/m{m}/completion.json' for m in MASSES},
        'observed_gof_qc': {f'observed/m{m}/gof_diagnostics.json' for m in MASSES},
        'bias': {f'bias/fits/m{m}/{label}/toy{toy:03d}/result.json'
                 for m in MASSES for label in ['r0', 'rExpected'] for toy in range(1, 51)},
        'gof': {f'gof/toy{toy:03d}/m{m}/result.json' for m in MASSES for toy in range(1, 101)},
        'profiles': {f'bias/profiles/m{m}/{label}/toy001/result.json'
                     for m in MASSES for label in ['r0', 'rExpected']},
    }
    globs = dict(observed='observed/m*/results.json', observed_outcomes='observed/m*/completion.json',
                 observed_gof_qc='observed/m*/gof_diagnostics.json', bias='bias/fits/m*/*/toy*/result.json',
                 gof='gof/toy*/m*/result.json', profiles='bias/profiles/m*/*/toy*/result.json')
    counts = {}
    for kind, expected in targets.items():
        actual = {str(p.relative_to(R)) for p in R.glob(globs[kind])}
        missing[kind] = sorted(expected - actual)
        counts[kind] = len(actual & expected)
        for path in sorted(actual - expected):
            issues.append({'path': path, 'issue': 'unexpected_output_ID'})
        for relative in sorted(actual & expected):
            path = R / relative
            try:
                q = read(path)
                if kind in ['observed', 'bias', 'gof', 'profiles']:
                    mass_part = next(v for v in path.parts if v.startswith('m') and v[1:] in MASSES)
                    if float(q['mass']) != float(mass_part[1:]):
                        raise ValueError('saved mass differs from path')
                if kind in ['bias', 'gof', 'profiles']:
                    toy_part = next(v for v in path.parts if v.startswith('toy') and v[3:].isdigit())
                    if q['toy'] != int(toy_part[3:]):
                        raise ValueError('saved toy ID differs from path')
                if kind == 'bias':
                    label = path.parent.parent.name
                    if q.get('injection', q.get('label')) != label:
                        raise ValueError('saved injection differs from path')
                    if q.get('fit_valid'):
                        if not (math.isfinite(q.get('edm', math.nan)) and 0 <= q['edm'] < .01):
                            raise ValueError('valid point has invalid EDM')
                        if not q.get('CR_valid'):
                            raise ValueError('valid point has invalid CR')
                    numerical['bias_invalid'] = numerical.get('bias_invalid', 0) + (not q.get('fit_valid'))
                if kind == 'gof':
                    qc = q.get('nominal_diagnostics', {})
                    if not (qc.get('outcome_complete') or qc.get('available')):
                        missing.setdefault('gof_qc', []).append(relative)
                    numerical['gof_invalid'] = numerical.get('gof_invalid', 0) + (not q.get('valid'))
            except Exception as error:
                issues.append({'path': relative, 'issue': str(error)})
    manifest = read(R / 'workspace_manifest.json')
    unchanged = []
    for q in manifest['workspaces']:
        path = Path(q['workspace'])
        good = hashlib.sha256(path.read_bytes()).hexdigest() == q['workspace_sha256']
        unchanged.append(good)
        if not good:
            issues.append({'path': str(path), 'issue': 'seed_workspace_hash_changed'})
    if final:
        summary = read(R / 'report_summary.json')
        if not summary.get('completed'):
            issues.append({'issue': 'report_not_complete'})
        if summary.get('rendering_errors'):
            issues.append({'issue': 'report_rendering_errors', 'details': summary['rendering_errors']})
        for name in ['toys_audit.json', 'generation_strength_audit.json']:
            if not read(R / name).get('passed'):
                issues.append({'issue': 'supporting_audit_not_passed', 'path': name})
    result = dict(updated_at=datetime.datetime.now().astimezone().isoformat(), final_audit=final,
                  counts=counts, targets={k: len(v) for k, v in targets.items()}, missing=missing,
                  issues=issues, seed_hashes_unchanged=all(unchanged), numerical_flags=numerical,
                  inventory_complete=not any(missing.values()), audit_passed=not issues,
                  complete=not issues and not any(missing.values()),
                  meaning='Execution completeness and provenance only; numerical/statistical validity remains separately flagged.')
    destination = R / 'completion_audit.json'
    tmp = destination.with_suffix(f'.{os.getpid()}.tmp')
    tmp.write_text(json.dumps(result, indent=2)); tmp.replace(destination)
    print(json.dumps({k: result[k] for k in ['counts', 'numerical_flags', 'seed_hashes_unchanged', 'inventory_complete', 'audit_passed', 'complete', 'issues']}, indent=2))
    return int(bool(issues) or (final and not result['inventory_complete']))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--final', action='store_true')
    raise SystemExit(main(parser.parse_args().final))
