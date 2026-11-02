"""Independent count/truth/quality audit of the new degree3 validation outputs."""
from pathlib import Path
import collections, csv, datetime, json, math, sys

P = Path.cwd()
R = P / 'runs/20261030_upper_envelope_scan'
partial = '--partial' in sys.argv
bias_paths = sorted((R / 'bias/degree3').glob('m*/*/m*/result.json'))
gof_paths = sorted((R / 'gof/degree3').glob('toy*/m*/result.json'))
bias = [json.loads(p.read_text()) for p in bias_paths]
gof = [json.loads(p.read_text()) for p in gof_paths]
assert all(q['degree'] == 3 for q in bias + gof)
assert all(math.isfinite(q['score']) for q in gof)
if not partial:
    assert len(bias) == 48, len(bias)
    assert len(gof) == 1600, len(gof)

bias_summary = []
median_changes = []
for path, q in zip(bias_paths, bias):
    D = path.parent.parent
    truth = json.loads((D / 'generation/truth.json').read_text())
    cr = json.loads((D / 'cr_result.json').read_text())
    meta = json.loads((D / 'degree3_interpretation.json').read_text())
    assert abs(q['rtrue'] - truth['rtrue']) < 1e-10 * max(1, abs(q['rtrue']))
    assert q['CR_valid'] == cr['valid']
    rmax = max(100, 5 * q['rtrue'])
    valid = q['status'] == 0 and q['covQual'] == 3 and q['error'] > 0 and abs(q['r']) < .99 * rmax and cr['valid']
    assert valid == q['valid']
    if q['bias_sigma'] is not None:
        assert abs(q['bias_sigma'] - (q['r'] - q['rtrue']) / q['error']) < 1e-10
    if q['label'] != 'r0':
        current = json.loads((R / f'degree3/m{q["mass"]:.1f}/results.json').read_text())['limits']['expected_50']
        old = meta['expected_median_at_generation']
        if abs(current - old) > 1e-10:
            median_changes.append(dict(mass=q['mass'],label=q['label'],used=old,current=current,relative_change=current / old - 1))

for label in ['r0', 'rExpected', 'r10Expected']:
    rows = [q for q in bias if q['label'] == label]
    assert len({q['mass'] for q in rows}) == len(rows)
    valid = [q for q in rows if q['valid']]
    bias_summary.append(dict(injection=label,total=len(rows),valid=len(valid),
                             within_0p2=sum(abs(q['bias_sigma']) < .2 for q in valid),
                             bias_min=min((q['bias_sigma'] for q in valid),default=None),
                             bias_max=max((q['bias_sigma'] for q in valid),default=None)))

failed_cr = []
cr_means_at_boundary = []
for path in sorted((R / 'gof/degree3').glob('toy*/cr_result.json')):
    q = json.loads(path.read_text())
    if not q['valid']:
        failed_cr.append(dict(toy=int(path.parent.name[3:]),status=q['selected']['status'],covQual=q['selected']['covQual']))
    for v in q['means']:
        if min(v['value'] - v['low'], v['high'] - v['value']) < 1e-4:
            cr_means_at_boundary.append(dict(toy=int(path.parent.name[3:]),name=v['name'],value=v['value']))

gof_summary = []
for mass_int in range(105, 121):
    mass = mass_int / 10
    rows = [q for q in gof if q['mass'] == mass]
    ids = {q['toy'] for q in rows}
    assert len(ids) == len(rows)
    if not partial:
        assert ids == set(range(1, 101)), (mass, len(ids))
    observed_file = R / f'degree3/m{mass:.1f}/results.json'
    if not rows or not observed_file.exists():
        continue
    observed = json.loads(observed_file.read_text())['gof_observed']
    exceed = sum(q['score'] > observed for q in rows)
    valid_exceed = sum(q['CR_valid'] and q['score'] > observed for q in rows)
    invalid = sum(not q['CR_valid'] for q in rows)
    numerical_invalid = sum((not q['CR_valid']) or q['r_boundary'] for q in rows)
    numerical_valid_exceed = sum(q['CR_valid'] and not q['r_boundary'] and q['score'] > observed for q in rows)
    n = len(rows)
    gof_summary.append(dict(mass=mass,n=n,exceed=exceed,p=exceed / n,
                            observed=observed,CR_invalid=invalid,
                            p_failure_low=valid_exceed / n,p_failure_high=(valid_exceed + invalid) / n,
                            numerical_invalid=numerical_invalid,
                            p_numerical_low=numerical_valid_exceed / n,
                            p_numerical_high=(numerical_valid_exceed + numerical_invalid) / n,
                            r_boundaries=sum(q['r_boundary'] for q in rows)))
for path, q in zip(gof_paths, gof):
    cr = json.loads((path.parent.parent / 'cr_result.json').read_text())
    assert q['CR_valid'] == cr['valid']
    assert q['r_boundary'] == (abs(q['r']) > 99)

def export(name, rows):
    if not rows:
        return
    with (R / name).open('w', newline='') as h:
        writer = csv.DictWriter(h, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

export('degree3_bias_results.csv', bias)
export('degree3_gof_scores.csv', gof)
export('degree3_gof_summary.csv', gof_summary)
audit = dict(status='partial' if partial else 'completed_with_documented_flags',
             audited_at_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
             bias_count=len(bias),bias_target=48,bias_valid=sum(q['valid'] for q in bias),
             bias_summary=bias_summary,bias_interpretation='Chebyshev3 generator to degree3 envelope: closure only; no alternative-background validation',
             gof_score_count=len(gof),gof_score_target=1600,gof_distinct_datasets=100,
             gof_dataset_reuse='Same100 datasets reused across16 mass hypotheses; correlated scores',
             gof_summary=gof_summary,gof_failed_CR_datasets=failed_cr,
             gof_CR_mean_boundary_instances=cr_means_at_boundary,
             gof_r_boundaries=sum(q['r_boundary'] for q in gof),
             expected_median_changed_after_injection=median_changes,
             notes=['No invalid result is silently dropped.',
                    'Conditional Hesse errors after frozen CR-only Upsilon means.',
                    'CR stage is repeated on every generated sample.',
                    'Finite scores are not a global multistart minimization guarantee.',
                    'Original degree1/2 results and their three failed CR datasets are unchanged.'])
(R / 'degree3_validation_audit.json').write_text(json.dumps(audit, indent=2))
print(json.dumps({k:v for k,v in audit.items() if k not in ['gof_summary','gof_CR_mean_boundary_instances','notes']},indent=2))
if not partial:
    config_path = R / 'degree3_validation_configuration.json'
    config = json.loads(config_path.read_text())
    config.update(status='completed_with_documented_flags',completed_at_utc=audit['audited_at_utc'],audit='degree3_validation_audit.json')
    config_path.write_text(json.dumps(config, indent=2))
