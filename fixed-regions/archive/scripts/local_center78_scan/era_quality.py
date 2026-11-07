import ROOT,json,pathlib
ROOT.gSystem.Load('libHiggsAnalysisCombinedLimit');R=pathlib.Path('/Users/imb/Documents/Codex/Analisi-statistica-Xee/runs/20261029_local_center78_scan');out=[]
for a in json.loads((R/'era_inputs.json').read_text()):
 d=pathlib.Path(a['path']).parent;f=ROOT.TFile.Open(str(d/'fitdiagnostics/fitDiagnostics_era.root'));a=dict(a)
 for st in ['b','s']:
  fit=f.Get('fit_'+st);a[st+'_status']=fit.status();a[st+'_covQual']=fit.covQual();a[st+'_edm']=fit.edm()
  if st=='s':
   r=fit.floatParsFinal().find('r');a['r']=r.getVal();a['r_error']=r.getError()
 f.Close();out.append(a)
(R/'era_quality.json').write_text(json.dumps(out,indent=2))
