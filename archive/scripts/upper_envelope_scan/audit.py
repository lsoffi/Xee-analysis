from pathlib import Path
import json,hashlib,math,csv,collections
import numpy as np,uproot
from scipy.stats import norm
P=Path.cwd();R=P/'runs/20261030_upper_envelope_scan'
a=[json.loads(p.read_text()) for p in sorted(R.glob('degree*/m*/results.json'))]
b=[json.loads(p.read_text()) for p in sorted(R.glob('bias/degree*/m*/*/m*/result.json'))]
g=[json.loads(p.read_text()) for p in sorted(R.glob('gof/degree*/toy*/m*/result.json'))]
sources={};checks=[];prediction_min=[];bands=[];physical=[];acceptance=[]
for q in a:
    D=R/f'degree{q["degree"]}/m{q["mass"]:.1f}'
    for reg in ['signal','control']:
        x=json.loads((D/reg/'audit.json').read_text());sources[x['source']]=x['source_sha256']
        for era in x['eras']:
            s=era['processes']['Zd'];acceptance.append(dict(degree=q['degree'],mass=q['mass'],region=reg,era=era['era'],full_rate=s['full_rate_r1'],fraction=s['fraction'],regional_rate=s['local_rate_r1']))
            for p in era['processes'].values():checks.extend(abs(z['analytic']-z['root_integral']) for z in p['acceptance_checks'])
    with uproot.open(D/'joint/fitdiagnostics/fitDiagnostics_upper.root') as f:
        prediction_min.append(min(float(min(f[f'shapes_fit_{s}/{reg}_ch{i}/total'].values())) for s in ['b','s'] for reg in ['SR','CR'] for i in range(1,5)))
    L=q['limits'];med=L['expected_50'];sigma=med/norm.ppf(.975)
    hit=[key for key,mul in [('expected_84',2),('expected_97p5',4)] if abs(L[key]/(med+mul*sigma)-1)<.003]
    bands.append(dict(degree=q['degree'],mass=q['mass'],suspected_bracket_hits=hit,ordered=all(L[x]<L[y] for x,y in zip(['expected_2p5','expected_16','expected_50','expected_84'],['expected_16','expected_50','expected_84','expected_97p5'])),status='provisional: upper-quantile crossing not independently validated'))
    for p in (D/'joint').glob('physical_branch*/results.json'):
        z=json.loads(p.read_text());z['directory']=str(p.parent.relative_to(R));physical.append(z)
hashes=[dict(source=p,expected=h,actual=hashlib.sha256(Path(p).read_bytes()).hexdigest()) for p,h in sources.items()]
shape_audit=json.loads((R/'extended_signal_shapes/audit.json').read_text())
hashes.append(dict(source=shape_audit['source'],expected=shape_audit['source_sha256'],actual=hashlib.sha256(Path(shape_audit['source']).read_bytes()).hexdigest()))
assert all(x['expected']==x['actual'] for x in hashes)
def export(name,rows):
    if not rows:return
    with (R/name).open('w',newline='') as h:
        writer=csv.DictWriter(h,fieldnames=rows[0].keys());writer.writeheader();writer.writerows(rows)
export('signal_acceptances.csv',acceptance)
export('bias_results.csv',b)
export('gof_toy_scores.csv',g)
rows=[]
for q in a:
    rows.append(dict(degree=q['degree'],mass=q['mass'],r=q['r'],r_error=q['r_error'],Z=q['Z'],p0=q['p0'],**q['limits'],FD_b_status=q['fit']['b']['status'],FD_b_cov=q['fit']['b']['covQual'],FD_s_status=q['fit']['s']['status'],FD_s_cov=q['fit']['s']['covQual'],CLs_verified=q['observed_limit_check']['verified'],unconstrained_FD_and_significance_agree=q['significance_consistent'],extrapolated=q['mass']>11))
export('observed_results.csv',rows)
stats=[]
for d in [1,2]:
    for label in ['r0','rExpected','r10Expected']:
        x=[q for q in b if q['degree']==d and q['label']==label];valid=[q for q in x if q['valid']]
        stats.append(dict(degree=d,injection=label,n=len(x),valid=len(valid),within_0p2=sum(abs(q['bias_sigma'])<.2 for q in valid),minimum=min([q['bias_sigma'] for q in valid],default=None),maximum=max([q['bias_sigma'] for q in valid],default=None)))
audit=dict(observed=len(a),bias=len(b),gof_scores=len(g),gof_target=3200,original_sources_unchanged=hashes,acceptance_checks=len(checks),max_acceptance_difference=max(checks),minimum_observed_bin_density=min(prediction_min),all_CLs_crossings_verified=all(q['observed_limit_check']['verified'] for q in a),all_observed_FD_status0_cov3=all(f['status']==0 and f['covQual']==3 for q in a for f in q['fit'].values()),expected_bands=bands,physical_branches=physical,bias_summary=stats,gof_CR_invalid=sum(not q['CR_valid'] for q in g),gof_r_boundaries=sum(q['r_boundary'] for q in g),notes=['Asimov bias, not ensemble mean or coverage.','Hesse errors conditional on frozen Upsilon means; numerical flags retained.','GoF common random numbers: 100 datasets per degree, 16 mass tests each; correlations across masses.','r-free saturated GoF can pass while background model produces spurious signal.','Expected upper bands are provisional, with suspected algorithm-bracket hits; failed alternative algorithm pilot not used.','Above 11 GeV relative lnN systematics held at 11; extrapolation model difference not included as nuisance.'])
(R/'validation.json').write_text(json.dumps(audit,indent=2));print({k:v for k,v in audit.items() if k not in ['original_sources_unchanged','expected_bands','physical_branches','notes']})
