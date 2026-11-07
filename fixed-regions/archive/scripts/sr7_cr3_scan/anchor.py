from pathlib import Path
import ROOT,json,subprocess,numpy as np
ROOT.gSystem.Load('libHiggsAnalysisCombinedLimit');P=Path.cwd();R=P/'runs/20261030_sr7_cr3_scan';d=R/'anchor';d.mkdir(exist_ok=True);f=ROOT.TFile.Open(str(R/'m7.5/joint/workspace.root'));w=f.Get('w');base=json.loads((P/'runs/20261030_fixed_sr_scan/m7.5/joint/quality.json').read_text());gen=json.loads((P/'runs/20261030_cheb7_bias_scan/generator_fit/fit.json').read_text())
for p in base['b']['parameters']:
 if w.var(p['name']):w.var(p['name']).setVal(p['value'])
for k,v in gen['parameters'].items():
 if k.startswith('local_') and not k.endswith('_CR') and w.var(k):w.var(k).setVal(v)
w.var('r').setVal(0);w.writeToFile(str(d/'workspace.root'));cmd=[str(P/'bin/combine'),'-M','FitDiagnostics','-d',str(d/'workspace.root'),'-m','7.5','--skipSBFit','--saveShapes','--saveNormalizations','--keepFailures','--minos','none','--cminDefaultMinimizerStrategy','1','--cminDefaultMinimizerTolerance','0.001','--X-rtd','MINIMIZER_freezeDisassociatedParams','-n','_anchor','-v','1'];(d/'command.json').write_text(json.dumps(cmd))
with (d/'combine.log').open('w') as out:subprocess.run(cmd,cwd=d,stdout=out,stderr=subprocess.STDOUT)
ff=ROOT.TFile.Open(str(d/'fitDiagnostics_anchor.root'));res=ff.Get('fit_b');params={p.GetName():p.getVal() for p in res.floatParsFinal()};mins=[]
for era in ['2022','2022EE','2023','2023BPix']:
 for suf,deg in [('',7),('_CR',3)]:
  c=[1]+[params[f'local_c{j}_{era}{suf}'] for j in range(1,deg+1)];roots=np.polynomial.chebyshev.chebroots(np.polynomial.chebyshev.chebder(c));pts=[-1,1]+[z.real for z in roots if abs(z.imag)<1e-8 and abs(z.real)<1];mins.append(float(min(np.polynomial.chebyshev.chebval(pts,c))))
a=dict(status=res.status(),covQual=res.covQual(),edm=res.edm(),minimum=min(mins),parameters=params);(d/'fit.json').write_text(json.dumps(a,indent=2));print(a['status'],a['covQual'],a['edm'],a['minimum'])
