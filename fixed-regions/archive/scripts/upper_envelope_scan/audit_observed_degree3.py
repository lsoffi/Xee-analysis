from pathlib import Path
import hashlib,json,math,re

P=Path.cwd();R=P/'runs/20261030_upper_envelope_scan/degree3';rows=[];sources={};acceptance=[]
for result in sorted(R.glob('m*/results.json')):
    a=json.loads(result.read_text());D=result.parent
    for region in ['signal','control']:
        b=json.loads((D/region/'audit.json').read_text());sources[b['source']]=b['source_sha256']
        for e in b['eras']:
            for proc in e['processes'].values():
                acceptance.extend(abs(q['analytic']-q['root_integral']) for q in proc['acceptance_checks'])
    L=a['limits'];sigma=L['expected_50']/1.959963984540054
    hits=[q for q,k in [('expected_84',2),('expected_97p5',4)] if abs(L[q]/(L['expected_50']+k*sigma)-1)<.003]
    log=(D/'joint/limits/combine.log').read_text()
    expected_warnings=[line.strip() for line in log.splitlines() if 'Best fit of asimov dataset' in line]
    refinement=D/'joint/refinement_degree3.json';spread=None
    if refinement.exists():
        trials=json.loads(refinement.read_text())['trials'];good=[q for q in trials if q['status']==0 and q['covQual']==3 and q['B_reference_matches'] and q['error']>0]
        if len(good)>1 and max(q['nll_nll0'] for q in good)-min(q['nll_nll0'] for q in good)<1e-5:spread=max(q['error'] for q in good)/min(q['error'] for q in good)-1
    rows.append(dict(mass=a['mass'],numerical_ok=a['numerical_ok'],r=a['r'],Z=a['Z'],p0=a['p0'],fd={s:dict(status=q['status'],covQual=q['covQual'],edm=q['edm'],minimum=q['minimum']) for s,q in a['fit'].items()},observed_limit_check=a['observed_limit_check'],expected_bracket_hits=hits,expected_warnings=expected_warnings,expected_validated=False,refinement_sigma_relative_spread=spread,refinement_sigma_unstable=spread is not None and spread>.1,significance_consistent=a['significance_consistent']))
hashes=[dict(source=p,sha256=h,unchanged=hashlib.sha256(Path(p).read_bytes()).hexdigest()==h) for p,h in sources.items()]
assert all(q['unchanged'] for q in hashes)
out=dict(completed=len(rows),target=16,numerical_ok=sum(q['numerical_ok'] for q in rows),numerical_ok_definition='FD status/covariance/positivity, observed CLs crossing, FD-vs-native significance; excludes validation of expected limits/bands',expected_validated=False,all_observed_crossings=all(q['observed_limit_check']['verified'] for q in rows),all_continua_positive=all(v['minimum']>0 for q in rows for v in q['fd'].values()),max_acceptance_error=max(acceptance,default=None),acceptance_checks=len(acceptance),input_hashes=hashes,results=rows,notes=['Expected bands retained native and provisional.','Only degree3; prior degree1/2 outputs untouched.','FitDiagnostics r may be negative; native significance quoted separately.'])
(R/'observed_audit.json').write_text(json.dumps(out,indent=2));print(json.dumps({k:v for k,v in out.items() if k not in ['input_hashes','results','notes']},indent=2));print('flags',[q['mass'] for q in rows if not q['numerical_ok']])
