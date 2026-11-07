import ROOT,json
from pathlib import Path
p=Path('low_cheb3_totalC0/results/m0.5');f=ROOT.TFile.Open(str(next(p.glob('fitDiagnostics*.root'))));wf=ROOT.TFile.Open(str(p/'workspace.root'));w=wf.Get('w')
for stage in 'bs':
 fit=f.Get('fit_'+stage);print(stage,len(fit.floatParsFinal()),[x.GetName() for x in fit.floatParsFinal()][:18])
for name in json.loads((p/'model.json').read_text())['constraints'][:2]:
 pdf=w.pdf(name);print(name,[(v.GetName(),v.getVal()) for v in pdf.servers()])
