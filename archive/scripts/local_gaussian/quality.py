import ROOT,pathlib,json,csv
ROOT.gSystem.Load('libHiggsAnalysisCombinedLimit');P=pathlib.Path('/Users/imb/Documents/Codex/Analisi-statistica-Xee');R=P/'runs/20261029_local_gaussian_7_8';f=ROOT.TFile.Open(str(R/'fitdiagnostics/fitDiagnostics_local.root'));out={};rows=[]
for st in ['b','s']:
 fit=f.Get('fit_'+st);pars=[]
 for v in fit.floatParsFinal():
  a=dict(stage=st,name=v.GetName(),value=v.getVal(),error=v.getError(),error_low=v.getAsymErrorLo(),error_high=v.getAsymErrorHi(),minimum=v.getMin(),maximum=v.getMax());pars.append(a);rows.append(a)
 out[st]=dict(status=fit.status(),covQual=fit.covQual(),edm=fit.edm(),parameters=pars)
(R/'quality.json').write_text(json.dumps(out,indent=2))
with (R/'parameters.csv').open('w') as h:
 w=csv.DictWriter(h,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
