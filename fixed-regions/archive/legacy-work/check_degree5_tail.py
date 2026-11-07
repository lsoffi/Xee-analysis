import ROOT,json
from pathlib import Path
ROOT.gSystem.Load('libHiggsAnalysisCombinedLimit');ROOT.gROOT.SetBatch(True)
wf=ROOT.TFile.Open('outputs/Xee_degree5_native_combined.root');w=wf.Get('w');ff=ROOT.TFile.Open('outputs/fitDiagnostics_degree5_native_combined.root');fit=ff.Get('fit_s');x=w.var('mass');obs=ROOT.RooArgSet(x);initial=w.allVars().snapshot();w.allVars().assignValueOnly(fit.floatParsFinal());results=[]
for era in ['2022','2022EE','2023','2023BPix']:
 pdf=w.pdf('dy_degree5_'+era)
 def bend():
  ys=[]
  for mass in [11,11.5,12]:x.setVal(mass);ys.append(pdf.getVal(obs))
  return ys[1]-.5*(ys[0]+ys[2]),ys[1]
 base,scale=bend();rows=[]
 for i in range(5):
  name=f't{i}_{era}_envelope';v=w.var(name);post=v.getVal();pre=initial.find(name).getVal();v.setVal(pre);b,_=bend();v.setVal(post);rows.append(dict(parameter=name,prefit=pre,postfit=post,tail_curvature_change=(base-b)/scale))
 rows.sort(key=lambda z:abs(z['tail_curvature_change']),reverse=True);results.append(dict(era=era,tail_curvature=base/scale,sensitivity=rows));print(era,[(z['parameter'],round(z['tail_curvature_change'],4)) for z in rows])
Path('outputs/degree5_tail_sensitivity.json').write_text(json.dumps(results,indent=2))
