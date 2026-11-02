from pathlib import Path
import json,uproot,numpy as np
P=Path('/Users/imb/Documents/Codex/Analisi-statistica-Xee');R=P/'runs/20261029_common95_cheb2';summary=[]
for f in R.glob('m*/results.json'):
 r=json.loads(f.read_text())
 with uproot.open(f.parent/'fitdiagnostics/fitDiagnostics_local.root') as root:
  delta=float(root['tree_fit_sb']['nll_nll0'].array(library='np')[0]);z=float(np.sqrt(max(0,-2*delta))) if r['r']['r']>0 else 0.
 checks=dict(Z_from_FitDiagnostics=z,Z_from_Significance=r.get('Z'),agreement_within_0p1=abs(z-r.get('Z',0))<.1)
 seeds={}
 for stage in ['b','s']:
  files=list((f.parent/f'significance_seed_{stage}').glob('higgsCombine*.root'))
  if files:
   with uproot.open(files[0]) as root:seeds[stage]=float(root['limit']['limit'].array(library='np')[0])
 checks['additional_initializations']=seeds
 checks['initializations_agree']=all(abs(v-r.get('Z',0))<.1 for v in seeds.values())
 r['consistency']=checks;r['numerical_ok']=r['fit_ok'] and checks['agreement_within_0p1'] and checks['initializations_agree'];f.write_text(json.dumps(r,indent=2));summary.append(dict(mass=r['mass'],**checks,numerical_ok=r['numerical_ok']))
(R/'consistency_checks.json').write_text(json.dumps(sorted(summary,key=lambda r:r['mass']),indent=2))
print('Consistent and good covariance',sum(r['numerical_ok'] for r in summary),'/',len(summary))
