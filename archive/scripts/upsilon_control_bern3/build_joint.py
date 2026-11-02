import ROOT,pathlib,json
ROOT.gROOT.SetBatch(True);ROOT.gSystem.Load('libHiggsAnalysisCombinedLimit')
P=pathlib.Path('/Users/imb/Documents/Codex/Analisi-statistica-Xee');R=P/'runs/20261029_upsilon_control_bern3';(R/'joint').mkdir(exist_ok=True)
sf=ROOT.TFile.Open(str(R/'signal/workspace.root'));sr=sf.Get('w');cf=ROOT.TFile.Open(str(R/'control/workspace.root'));cr=cf.Get('w');nw=ROOT.RooWorkspace('w');keep=[]
# Remove factory EDIT alias bookkeeping before a second graph customization.
for arg in list(cr.allFunctions())+list(cr.allVars())+list(cr.allPdfs()):
 for attr in list(arg.attributes()):
  if str(attr).startswith('ORIGNAME:'):arg.setAttribute(str(attr),False)
# Customize CR graph: separate observable and continuum, share the original physical parameters.
custom=ROOT.RooCustomizer(cr.obj('ModelConfig_bonly').GetPdf(),'CR')
crmass=ROOT.RooRealVar('mass_CR','mass_CR',9.8,8.8,10.8);crmass.setBins(150);custom.replaceArg(cr.var('mass'),crmass);keep.append(crmass)
for era in ['2022','2022EE','2023','2023BPix']:
 for prefix in ['local_b1_', 'local_b2_', 'local_b3_', 'local_yield_']:
  old=cr.var(prefix+era);new=old.clone(old.GetName()+'_CR');custom.replaceArg(old,new);keep.append(new)
 for proc in ['upsilon1s','upsilon2s']:
  old=cr.function('local_acceptance_'+proc+'_'+era);new=old.clone(old.GetName()+'_CR');custom.replaceArg(old,new);keep.append(new)
cp=custom.build();keep.extend([custom,cp]);print('CR extra',[(x.GetName(),x.ClassName()) for x in cp._extraConstraints])
cat=ROOT.RooCategory('CMS_channel','CMS_channel')
for region in ['SR','CR']:
 for i in range(1,5):cat.defineType(region+'_ch'+str(i))
for mcname in ['ModelConfig','ModelConfig_bonly']:
 source=sr.obj(mcname);sim=ROOT.RooSimultaneousOpt('joint_'+mcname,'',cat)
 for i in range(1,5):
  ch='ch'+str(i);sim.addPdf(source.GetPdf().getPdf(ch),'SR_'+ch);sim.addPdf(cp.getPdf(ch),'CR_'+ch)
 for c in source.GetPdf()._extraConstraints:sim.addExtraConstraint(c)
 for c in cp._extraConstraints:sim.addExtraConstraint(c)
 getattr(nw,'import')(sim,ROOT.RooFit.RecycleConflictNodes());keep.append(sim)
 mc=ROOT.RooStats.ModelConfig(mcname,nw);mc.SetPdf(nw.pdf(sim.GetName()));mc.SetObservables(ROOT.RooArgSet(nw.var('mass'),nw.var('mass_CR'),nw.cat('CMS_channel')))
 for getter,setter in [('GetParametersOfInterest','SetParametersOfInterest'),('GetNuisanceParameters','SetNuisanceParameters'),('GetGlobalObservables','SetGlobalObservables')]:
  args=ROOT.RooArgSet()
  for v in getattr(source,getter)():
   if nw.arg(v.GetName()):args.add(nw.arg(v.GetName()))
  getattr(mc,setter)(args)
 getattr(nw,'import')(mc);keep.append(mc)
# Preserve Combine process metadata for saveShapes.
for arg in list(nw.allFunctions())+list(nw.allVars()):
 channel=arg.getStringAttribute('combine.channel')
 if channel and channel.startswith('ch'):arg.setStringAttribute('combine.channel',('CR_' if arg.GetName().endswith('_CR') else 'SR_')+channel)
weight=ROOT.RooRealVar('evt_weight','',0,1e9);obs=ROOT.RooArgSet(nw.var('mass'),nw.var('mass_CR'),nw.cat('CMS_channel'),weight);data=ROOT.RooDataSet('data_obs','',obs,ROOT.RooFit.WeightVar(weight));counts={}
for region,w,massname in [('SR',sr,'mass'),('CR',cr,'mass_CR')]:
 original=w.data('data_obs');nw.var('mass').setVal(7.5);nw.var('mass_CR').setVal(9.8)
 for j in range(original.numEntries()):
  row=original.get(j);channel=region+'_'+row.find('CMS_channel').getLabel();nw.cat('CMS_channel').setLabel(channel);nw.var(massname).setVal(row.getRealValue('mass'));weight.setVal(original.weight());data.add(obs,original.weight());counts[channel]=counts.get(channel,0)+original.weight()
getattr(nw,'import')(data);assert abs(data.sumEntries()-sr.data('data_obs').sumEntries()-cr.data('data_obs').sumEntries())<1e-8
nw.var('mass').setVal(7.5);nw.var('mass_CR').setVal(9.8);nw.var('r').setRange(-100,100);nw.var('r').setVal(0)
# Verify no duplicate copies of original nuisance constraints and correct rate transfer.
constraints=nw.obj('ModelConfig').GetPdf().getAllConstraints(nw.obj('ModelConfig').GetObservables(),nw.obj('ModelConfig').GetNuisanceParameters(),False)
print('ALL CONSTRAINTS',[v.GetName() for v in constraints]);assert len([v for v in constraints if v.GetName().startswith('CMS_scale_e_Pdf')])==1
for i,era in enumerate(['2022','2022EE','2023','2023BPix'],1):
 for proc in ['upsilon1s','upsilon2s']:
  scale=nw.var('CMS_EXO25020_bkgScale_'+proc+'_'+era);assert scale and not scale.isConstant()
  srate=nw.function('local_rate_'+proc+'_'+era);crate=nw.function('local_rate_'+proc+'_'+era+'_CR');assert srate and crate
  assert abs(srate.getVal()-sr.function(srate.GetName()).getVal())<1e-8
  assert abs(crate.getVal()-cr.function(srate.GetName()).getVal())<1e-8
  print('RATES',era,proc,srate.getVal(),crate.getVal())
nw.writeToFile(str(R/'joint/workspace.root'));(R/'joint/audit.json').write_text(json.dumps(dict(SR_range=[sr.var('mass').getMin(),sr.var('mass').getMax()],CR_range=[8.8,10.8],counts=counts,events=data.sumEntries(),constraints=[v.GetName() for v in constraints],independent_continua=True,shared_upsilon_normalizations=True,CR_signal_omitted_negligible_gaussian_acceptance=True),indent=2));print('DONE')
