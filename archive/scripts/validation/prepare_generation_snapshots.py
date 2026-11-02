import ROOT,pathlib,json
ROOT.gSystem.Load('libHiggsAnalysisCombinedLimit');ROOT.gROOT.SetBatch(True)
P=pathlib.Path('/Users/imb/Documents/Codex/Analisi-statistica-Xee');R=P/'runs/20261029_background_validation';audit=[]
for mass in ['7.5','9.5']:
 for model in ['degree5','cheb7','bernstein']:
  path=R/f'workspaces/m{mass}_{model}.root';f=ROOT.TFile.Open(str(path));w=f.Get('w');d=w.genobj('discreteParams')
  if d:d.removeAll()
  w.writeToFile(str(path.with_suffix('.clean.root')))
 for model in ['cheb7','bernstein']:
  f=ROOT.TFile.Open(str(R/f'generators/m{mass}_{model}/higgsCombine_validation.MultiDimFit.mH{mass}.root'));w=f.Get('w');assert w.loadSnapshot('MultiDimFit')
  ff=ROOT.TFile.Open(str(R/f'generators/m{mass}_{model}/multidimfit_validation.root'));fit=ff.Get('fit_mdf');assert fit.status()==0
  d=w.genobj('discreteParams')
  if d:d.removeAll()
  w.var('r').setConstant(False);w.var('r').setVal(0);w.saveSnapshot('validation_generator',w.allVars())
  path=R/f'workspaces/generator_m{mass}_{model}.root';w.writeToFile(str(path));audit.append(dict(mass=mass,model=model,status=fit.status(),covQual=fit.covQual(),path=str(path)))
(R/'generator_audit.json').write_text(json.dumps(audit,indent=2));print(json.dumps(audit,indent=2))
