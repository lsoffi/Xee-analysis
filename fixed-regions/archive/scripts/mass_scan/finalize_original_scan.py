#!/usr/bin/env python3
import pathlib,json,csv,hashlib,math
import numpy as np,uproot
P=pathlib.Path('/Users/imb/Documents/Codex/Analisi-statistica-Xee');R=P/'runs/20261029_original_mass_scan';cache=json.loads((R/'summary/processed.json').read_text());rows=list(cache.values());assert len(rows)==357,len(rows)
assert not any(r.get('postprocess_error') for r in rows)
inv=json.loads((R/'input_inventory.json').read_text());checks=[dict(path=a['path'],unchanged=hashlib.sha256(pathlib.Path(a['path']).read_bytes()).hexdigest()==a['sha256']) for a in inv];assert all(c['unchanged'] for c in checks)
(R/'summary/input_hashes_verified.json').write_text(json.dumps(checks,indent=2))
limits=sorted([r for r in rows if r['role']=='limits'],key=lambda a:a['mass']);sig=sorted([r for r in rows if r['role']=='significance'],key=lambda a:a['mass']);fd=[r for r in rows if r['role']=='fitdiagnostics'];
for r in fd:
 if r.get('r') is not None:
  r['r_interval_hits_lower_bound']=r['r']-abs(r.get('rLoErr') or 0)<=-30+1e-3
  r['r_interval_hits_upper_bound']=r['r']+abs(r.get('rHiErr') or 0)>=30-1e-3
  r['r_bestfit_at_bound']=min(abs(r['r']+30),abs(r['r']-30))<1e-3
(R/'summary/processed.json').write_text(json.dumps(cache,indent=2))
def writecsv(name,rs):
 keys=list(dict.fromkeys(k for r in rs for k in r))
 with open(R/'summary'/name,'w') as f:
  w=csv.DictWriter(f,fieldnames=keys);w.writeheader();w.writerows(rs)
writecsv('all_jobs.csv',rows);writecsv('fit_diagnostics.csv',fd)
warnings=[r for r in fd if not r.get('valid') or r.get('r_interval_hits_lower_bound') or r.get('r_interval_hits_upper_bound') or r.get('r_bestfit_at_bound')]
writecsv('fit_warnings.csv',warnings) if warnings else None
assert len(fd)==255 and len(limits)==51 and len(sig)==51
with uproot.recreate(R/'summary/scan_results.root') as f:
 limkeys=['mass','observed','expected_2p5','expected_16','expected_50','expected_84','expected_97p5']
 f['limits_scan']={k:np.array([r.get(k) if r.get(k) is not None else np.nan for r in limits],dtype='float64') for k in limkeys};f['limits_scan_valid']={'mass':np.array([r['mass'] for r in limits]),'valid':np.array([r.get('valid',False) for r in limits])}
 f['significance_scan']={k:np.array([r.get(k) if r.get(k) is not None else np.nan for r in sig],dtype='float64') for k in ['mass','Z_local','p0_local','log10_p0_local']}
plotcounts=dict(png=len(list((R/'fits').rglob('fits.png'))),pdf=len(list((R/'fits').rglob('fits.pdf'))));assert plotcounts==dict(png=255,pdf=255),plotcounts
observations={(r['mass_label'],r['era']):r['sum_weights'] for r in inv};verified=0;anomalies=[]
for r in fd:
 wd=R/r['directory'];files=list(wd.glob('*_fit_s_shapes.csv'));assert len(files)==(4 if r['era']=='allYears' else 1)
 for f in files:
  a=np.genfromtxt(f,delimiter=',',names=True);assert len(a)==450;
  if not np.isfinite(a['total']).all() or not (a['total']>=0).all():anomalies.append(str(f.relative_to(R)))
  assert abs(a['low'][0]-6)<1e-9 and abs(a['high'][-1]-12)<1e-9
  if r['era']=='allYears':
   era={'ch1':'2022','ch2':'2022EE','ch3':'2023','ch4':'2023BPix'}[f.name.split('_')[0]]
  else:era=r['era']
  assert abs(a['data'].sum()-observations[r['mass_label'],era])<1e-5,(r['mass_label'],era,a['data'].sum())
  verified+=1
valids=[r for r in sig if r.get('valid')];peak=max(valids,key=lambda r:r['Z_local']) if valids else None
qc=dict(jobs=len(rows),fitdiagnostics=len(fd),limits=len(limits),significance=len(sig),plotcounts=plotcounts,shape_tables_verified=verified,input_files_unchanged=len(checks),shape_prediction_anomalies=anomalies,diagnostic_interval_boundary_warnings=sum(bool(r.get('r_interval_hits_lower_bound') or r.get('r_interval_hits_upper_bound')) for r in fd),execution_failures=sum(r['exit_code']!=0 for r in rows),diagnostic_quality_warnings=sum(not r.get('valid',False) for r in fd),limit_valid=sum(r.get('valid',False) for r in limits),significance_valid=sum(r.get('valid',False) for r in sig),max_local_significance={k:peak[k] for k in ['mass','Z_local','p0_local']} if peak else None)
(R/'summary/final_verification.json').write_text(json.dumps(qc,indent=2))
(P/'runs/LATEST_ORIGINAL_MASS_SCAN.txt').write_text(str(R)+'\n')
print(json.dumps(qc,indent=2))
