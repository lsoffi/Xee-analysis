import sys
MASS=sys.argv[1]
MODEL='chebyshev2'
import ROOT,pathlib,json,csv
ROOT.gSystem.Load('libHiggsAnalysisCombinedLimit');P=pathlib.Path('/Users/imb/Documents/Codex/Analisi-statistica-Xee');R=P/'runs/20261030_sr6_cr4_scan'/('m'+MASS)/'joint';R.mkdir(parents=True,exist_ok=True);f=ROOT.TFile.Open(str(R/'fitdiagnostics/fitDiagnostics_local.root'));out={};rows=[]
for st in ['b','s']:
 fit=f.Get('fit_'+st);pars=[]
 for v in fit.floatParsFinal():
  a=dict(stage=st,name=v.GetName(),value=v.getVal(),error=v.getError(),error_low=v.getAsymErrorLo(),error_high=v.getAsymErrorHi(),minimum=v.getMin(),maximum=v.getMax());pars.append(a);rows.append(a)
 out[st]=dict(status=fit.status(),covQual=fit.covQual(),edm=fit.edm(),parameters=pars)
(R/'quality.json').write_text(json.dumps(out,indent=2))
with (R/'parameters.csv').open('w') as h:
 w=csv.DictWriter(h,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
import numpy as np,hashlib
wf=ROOT.TFile.Open(str(R/'workspace.root'));w=wf.Get('w');audits={region:json.loads((R.parent/region/'audit.json').read_text()) for region in ['signal','control']}
for a in audits.values():assert hashlib.sha256(pathlib.Path(a['source']).read_bytes()).hexdigest()==a['source_sha256']
w.var('r').setVal(1)
for i,era in enumerate(['2022','2022EE','2023','2023BPix']):
 for region,suffix in [('signal',''),('control','_CR')]:
  expected=audits[region]['eras'][i]['processes']['Zd']['local_rate_r1'];actual=w.function('local_rate_Zd_'+era+suffix).getVal();assert abs(actual-expected)<1e-7*max(1,abs(expected)),(MASS,era,region,actual,expected)
 assert w.obj('ModelConfig').GetPdf().getPdf('CR_ch'+str(i+1)).dependsOn(w.var('r'))
for st in ['b','s']:
 for a in out[st]['parameters']:
  if w.var(a['name']):w.var(a['name']).setVal(a['value'])
 minima=[]
 for era in ['2022','2022EE','2023','2023BPix']:
  for suf in ['','_CR']:
   c=[1]+[w.var(f'local_c{j}_{era}{suf}').getVal() for j in range(1,7 if suf=='' else 5)];roots=np.polynomial.chebyshev.chebroots(np.polynomial.chebyshev.chebder(c));pts=[-1,1]+[r.real for r in roots if abs(r.imag)<1e-9 and -1<r.real<1];minima.append(float(min(np.polynomial.chebyshev.chebval(pts,c))))
 out[st]['continuum_min']=min(minima)
(R/'quality.json').write_text(json.dumps(out,indent=2))
