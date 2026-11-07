from pathlib import Path
import ROOT,json,hashlib,sys
ROOT.gROOT.SetBatch(True);ROOT.gSystem.Load('libHiggsAnalysisCombinedLimit');ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.WARNING)
P=Path.cwd();SRC=P/'runs/20261030_sr6_cr4_scan';OUT=P/'runs/20261030_gof_thesis48'
eras=['2022','2022EE','2023','2023BPix'];tasks=[];reference={}
for k in range(60,91):
 mass=k/10;source=SRC/f'm{mass:.1f}/joint/workspace.root';f=ROOT.TFile.Open(str(source));w=f.Get('w');ff=ROOT.TFile.Open(str(SRC/f'm{mass:.1f}/joint/fitdiagnostics/fitDiagnostics_local.root'));fit=ff.Get('fit_b');assert fit.status()==0
 for v in fit.floatParsFinal():
  if w.var(v.GetName()):w.var(v.GetName()).setVal(v.getVal())
 for i,era in enumerate(eras,1):
  nw=ROOT.RooWorkspace('w');keep=[];cat=ROOT.RooCategory('CMS_channel','');channels=[f'SR_ch{i}',f'CR_ch{i}']
  for j,ch in enumerate(channels):cat.defineType(ch,j)
  for mcname in ['ModelConfig','ModelConfig_bonly']:
   mc=w.obj(mcname);sim=ROOT.RooSimultaneousOpt('era_'+mcname,'',cat)
   for ch in channels:sim.addPdf(mc.GetPdf().getPdf(ch),ch)
   for name in ['positive_continuum_'+era,'positive_continuum_'+era+'_CR']:sim.addExtraConstraint(w.pdf(name))
   getattr(nw,'import')(sim,ROOT.RooFit.RecycleConflictNodes());keep.append(sim);nmc=ROOT.RooStats.ModelConfig(mcname,nw);nmc.SetPdf(nw.pdf(sim.GetName()));nmc.SetObservables(ROOT.RooArgSet(nw.var('mass'),nw.var('mass_CR'),nw.cat('CMS_channel')))
   for getter,setter in [('GetParametersOfInterest','SetParametersOfInterest'),('GetNuisanceParameters','SetNuisanceParameters'),('GetGlobalObservables','SetGlobalObservables')]:
    nodes=ROOT.RooArgSet()
    for v in getattr(mc,getter)():
     if nw.arg(v.GetName()):nodes.add(nw.arg(v.GetName()))
    getattr(nmc,setter)(nodes)
   getattr(nw,'import')(nmc);keep.append(nmc)
  weight=ROOT.RooRealVar('evt_weight','',0,1e9);obs=ROOT.RooArgSet(nw.var('mass'),nw.var('mass_CR'),nw.cat('CMS_channel'),weight);data=ROOT.RooDataSet('data_obs','',obs,ROOT.RooFit.WeightVar(weight));contents=[];source_data=w.data('data_obs')
  nw.var('mass').setVal(7.5);nw.var('mass_CR').setVal(9.8)
  for j in range(source_data.numEntries()):
   row=source_data.get(j);ch=row.find('CMS_channel').getLabel()
   if ch not in channels:continue
   observable='mass_CR' if ch.startswith('CR_') else 'mass';x=row.getRealValue(observable);wt=source_data.weight();nw.cat('CMS_channel').setLabel(ch);nw.var(observable).setVal(x);weight.setVal(wt);data.add(obs,wt);contents.append((ch,round(x,12),wt))
  getattr(nw,'import')(data);nw.var('mass').setVal(7.5);nw.var('mass_CR').setVal(9.8);nw.var('r').setVal(0)
  disc=ROOT.RooArgSet(nw.cat('pdfindex_'+era));disc.setName('discreteParams');getattr(nw,'import')(disc,'discreteParams');nw.saveSnapshot('background_fit',nw.allVars())
  assert data.numEntries()==338 and nw.genobj('discreteParams').getSize()==1
  for pr in ['upsilon1s','upsilon2s']:
   mu=nw.var(f'Ymean_{pr}_{era}');assert not mu.isConstant() and nw.function(f'Ysigma_{pr}_{era}').isConstant()
   for ch in channels:assert nw.obj('ModelConfig').GetPdf().getPdf(ch).dependsOn(mu)
  dh=hashlib.sha256(json.dumps(contents).encode()).hexdigest()
  if era not in reference:reference[era]=dh
  assert dh==reference[era]
  # Compare B-only bin intensities to the original workspace at identical parameters.
  max_delta=0.;normsrc=ROOT.RooArgSet(w.var('mass'),w.var('mass_CR'));normnew=ROOT.RooArgSet(nw.var('mass'),nw.var('mass_CR'))
  for ch,x,wt in contents:
   observable='mass_CR' if ch.startswith('CR_') else 'mass';w.var(observable).setVal(x);nw.var(observable).setVal(x)
   a=w.obj('ModelConfig_bonly').GetPdf().getPdf(ch);b=nw.obj('ModelConfig_bonly').GetPdf().getPdf(ch)
   va=a.getVal(normsrc)*a.expectedEvents(normsrc);vb=b.getVal(normnew)*b.expectedEvents(normnew);max_delta=max(max_delta,abs(va-vb)/max(1,abs(va)))
  assert max_delta<1e-9,(era,mass,max_delta)
  nw.var('mass').setVal(7.5);nw.var('mass_CR').setVal(9.8)
  d=OUT/era/f'm{mass:.1f}';d.mkdir(parents=True,exist_ok=True);nw.writeToFile(str(d/'workspace.root'))
  audit=dict(era=era,mass=mass,source=str(source),source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),data_hash=dh,events=data.sumEntries(),bins=338,channels=channels,max_relative_bin_intensity_difference=max_delta,signal_strength_fixed=0,SR=[6,9.093333333333334],CR=[9.093333333333334,10.506666666666668],discrete_indices=1)
  (d/'audit.json').write_text(json.dumps(audit,indent=2));tasks.append(dict(era=era,mass=mass,path=str(d),seed=483000+1000*i+k))
 ff.Close();f.Close();print('Prepared',mass,flush=True)
(OUT/'tasks.json').write_text(json.dumps(tasks,indent=2))
