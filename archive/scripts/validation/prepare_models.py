import ROOT,json,pathlib
ROOT.gROOT.SetBatch(True);ROOT.gSystem.Load('libHiggsAnalysisCombinedLimit')
P=pathlib.Path('/Users/imb/Documents/Codex/Analisi-statistica-Xee');out=P/'runs/20261029_background_validation/workspaces';audit=[]
for mass in ['7.5','9.5']:
 src=P/f'data/mass_scan/261028/cards_region2_data_envelope_allCorrections_binned/ee/{mass}/Xee_ee_4_allYears.root';f=ROOT.TFile.Open(str(src))
 for model in ['original','degree5','cheb7','bernstein']:
  f=ROOT.TFile.Open(str(src));w=f.Get('w');repl=[]
  for i,era in enumerate(['2022','2022EE','2023','2023BPix'],1):
   old=w.pdf(f'shapeBkg_dy_ch{i}')
   if model=='degree5':
    coeff=ROOT.RooArgList(*[w.var(f't{j}_{era}_envelope') for j in range(5)])
    pdf=ROOT.RooChebychev(f'validation_degree5_{era}','',w.var('mass'),coeff);getattr(w,'import')(pdf)
    repl.append(f'{old.GetName()}={pdf.GetName()}')
   elif model in ['cheb7','bernstein']:
    pdf=old.getPdf(1 if model=='cheb7' else 0);repl.append(f'{old.GetName()}={pdf.GetName()}')
  if repl:
   for label,config in [('model_s','ModelConfig'),('model_b','ModelConfig_bonly')]:
    new=w.factory(f'EDIT::{label}_validation({label},'+','.join(repl)+')');assert new
    mc=w.obj(config);mc.SetPdf(new)
    params=new.getParameters(mc.GetObservables());nuis=ROOT.RooArgSet()
    for v in mc.GetNuisanceParameters():
     if params.find(v.GetName()) and not v.isConstant():nuis.add(v)
    mc.SetNuisanceParameters(nuis)
   if w.obj('discreteParams'):w.defineSet('discreteParams',ROOT.RooArgSet())
  w.var('r').setRange(-100,100)
  dest=out/f'm{mass}_{model}.root';w.writeToFile(str(dest))
  audit.append(dict(mass=mass,model=model,path=str(dest),events=w.data('data_obs').sumEntries(),bins=w.var('mass').getBins(),nuisance_count=w.obj('ModelConfig').GetNuisanceParameters().getSize(),pdf=w.obj('ModelConfig').GetPdf().GetName()))
(out/'audit.json').write_text(json.dumps(audit,indent=2))
