"""Check completeness, provenance, and retained numerical flags of the toy study."""
from pathlib import Path
from collections import Counter
import hashlib
import json
import math

P = Path(__file__).resolve().parents[2]
R = P / 'runs/20261030_upper_bias_cheb5_toys'


def main():
    config = json.loads((R / 'configuration.json').read_text())
    selected = config['selected_toy_ids']
    assert selected == list(range(1, 51))
    sources_changed = [name for name, digest in config['source_hashes'].items()
                       if hashlib.sha256((P / name).read_bytes()).hexdigest() != digest]
    missing, mismatches, results, seeds, data_hashes = [], [], [], [], []
    interval_methods = Counter()
    flags = Counter()
    counts = []
    for mass in config['masses']:
        mass = f'{mass:.1f}'
        for injection in config['injections']:
            generation = R / f'generation/m{mass}/{injection}'
            truth = json.loads((generation / 'truth.json').read_text())
            seeds.append(truth['seed'])
            audit = {a['toy']: a for a in truth['toy_audit']}
            assert len(audit) == truth['ntoys'] == 100
            data_hashes += [audit[k]['data_sha256'] for k in selected]
            expected_r = 0. if injection == 'r0' else config['injected_expected_medians'][mass]
            if abs(truth['rtrue'] - expected_r) > 1e-10:
                mismatches.append(f'{mass}/{injection}: generation rtrue')
            cell = []
            for toy in selected:
                directory = R / f'fits/m{mass}/{injection}/toy{toy:03d}'
                path = directory / 'result.json'
                if not path.exists():
                    missing.append(str(path.relative_to(R)))
                    continue
                row = json.loads(path.read_text())
                # Keep only the counting fields; the full ROOT-fit metadata is
                # validated one record at a time rather than retained in RAM.
                compact = {k: row.get(k) for k in ['valid', 'fit_valid', 'protocol']}
                cell.append(compact)
                results.append(compact)
                interval_methods[row.get('interval_method')] += 1
                flags.update(row.get('failure_flags', []))
                if (row['toy'] != toy or row['mass'] != float(mass)
                        or row['injection'] != injection or row.get('rtrue') != expected_r):
                    mismatches.append(str(path.relative_to(R)) + ': key/truth mismatch')
                if row.get('toy_events') != audit[toy]['events']:
                    mismatches.append(str(path.relative_to(R)) + ': event count mismatch')
                if row.get('valid'):
                    if not (row.get('CR_valid') and row.get('fit_valid') and row.get('interval_valid')):
                        mismatches.append(str(path.relative_to(R)) + ': validity mismatch')
                    lo, hi, fit = row['interval_low'], row['interval_high'], row['r']
                    if not (all(math.isfinite(v) for v in [lo, hi, fit]) and lo < fit < hi):
                        mismatches.append(str(path.relative_to(R)) + ': interval ordering')
                    if row['covered68'] != (lo <= expected_r <= hi):
                        mismatches.append(str(path.relative_to(R)) + ': coverage mismatch')
                    if toy == 1 and row.get('interval_method') != 'profile68':
                        mismatches.append(str(path.relative_to(R)) + ': unexpected interval method')
            counts.append(dict(mass=float(mass), injection=injection, completed=len(cell),
                               target=50, valid=sum(bool(a.get('valid')) for a in cell),
                               point_valid=sum(bool(a.get('fit_valid')) for a in cell),
                               protocol_v2=sum(a.get('protocol') == 'two_fixed_anchors_v2' for a in cell),
                               failures=sum(not a.get('valid') for a in cell)))
    protocol_complete = all(a.get('protocol') == 'two_fixed_anchors_v2' for a in results)
    output = dict(target=1600, completed=len(results), valid=sum(bool(a.get('valid')) for a in results),
                  point_valid=sum(bool(a.get('fit_valid')) for a in results),
                  protocol_v2=sum(a.get('protocol') == 'two_fixed_anchors_v2' for a in results),
                  complete=len(results) == 1600 and not missing and protocol_complete,
                  provenance_valid=not sources_changed and not mismatches,
                  sources_changed=sources_changed, mismatches=mismatches,
                  missing=missing, by_cell=counts, failure_flags=dict(flags),
                  interval_methods=dict(interval_methods),
                  interval_policy=config['interval_policy'], ensemble_coverage_evaluated=False,
                  independent_seed_count=len(set(seeds)), expected_seed_count=32,
                  selected_unique_dataset_hashes=len(set(data_hashes)),
                  selected_toys_per_cell=selected, unused_generated_ids=list(range(51, 101)),
                  interpretation='Failures are retained and counted; numerical validity does not itself establish statistical unbiasedness.')
    assert len(set(seeds)) == 32
    assert len(set(data_hashes)) == 1600
    (R / 'audit.json').write_text(json.dumps(output, indent=2))
    print(json.dumps({k: v for k, v in output.items() if k not in ['missing', 'by_cell', 'selected_toys_per_cell', 'unused_generated_ids']}))


if __name__ == '__main__':
    main()
