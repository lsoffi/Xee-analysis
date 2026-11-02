import ROOT,pathlib,json
ROOT.gROOT.SetBatch(True);ROOT.gSystem.Load('libHiggsAnalysisCombinedLimit');P=pathlib.Path('/Users/imb/Documents/Codex/Analisi-statistica-Xee');R=P/'runs/20261029_local_center78_scan';eras=['2022','2022EE','2023','2023BPix'];tasks=[]
for model in ['exp','bernstein2','chebyshev2']:
 f=ROOT.TFile.Open(str(R/model/'m7.8/workspace.root'));w=f.Get('w')
 for i,era in enumerate(eras,1):
  nw=ROOT.RooWorkspace('w');keep=[];cat=ROOT.RooCategory('CMS_channel','');w.cat('CMS_channel').setLabel('ch'+str(i));cat.defineType('ch'+str(i),0)
  for mcname in ['ModelConfig','ModelConfig_bonly']:
   mc=w.obj(mcname);pdf=mc.GetPdf().getPdf('ch'+str(i));sim=ROOT.RooSimultaneous('single_'+mcname,'',cat);sim.addPdf(pdf,'ch'+str(i));getattr(nw,'import')(sim,ROOT.RooFit.RecycleConflictNodes());keep.append(sim);newmc=ROOT.RooStats.ModelConfig(mcname,nw);newmc.SetPdf(nw.pdf(sim.GetName()));newmc.SetObservables(ROOT.RooArgSet(nw.var('mass'),nw.cat('CMS_channel')))
   for getter,setter in [('GetParametersOfInterest','SetParametersOfInterest'),('GetNuisanceParameters','SetNuisanceParameters'),('GetGlobalObservables','SetGlobalObservables')]:
    nodes=ROOT.RooArgSet()
    for v in getattr(mc,getter)():
     target=nw.arg(v.GetName())
     if target:nodes.add(target)
    getattr(newmc,setter)(nodes)
   getattr(nw,'import')(newmc);keep.append(newmc)
  selected=w.data('data_obs').reduce(ROOT.RooFit.Cut('CMS_channel==CMS_channel::ch'+str(i)))
  weight=ROOT.RooRealVar('evt_weight','',0,1e9);obs=ROOT.RooArgSet(nw.var('mass'),nw.cat('CMS_channel'),weight);data=ROOT.RooDataSet('data_obs','',obs,ROOT.RooFit.WeightVar(weight))
  for j in range(selected.numEntries()):
   row=selected.get(j);nw.var('mass').setVal(row.getRealValue('mass'));nw.cat('CMS_channel').setLabel('ch'+str(i));weight.setVal(selected.weight());data.add(obs,selected.weight())
  assert abs(data.sumEntries()-selected.sumEntries())<1e-8
  getattr(nw,'import')(data);nw.var('mass').setVal(7.8);nw.var('r').setRange(-30,30);nw.var('r').setVal(0)
  dest=R/'individual_eras'/model/era;dest.mkdir(parents=True,exist_ok=True);nw.writeToFile(str(dest/'workspace.root'));tasks.append(dict(model=model,era=era,path=str(dest/'workspace.root'),events=data.sumEntries()))
 f.Close()
(R/'era_inputs.json').write_text(json.dumps(tasks,indent=2))
