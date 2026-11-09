import ROOT as R,json
from pathlib import Path
R.gSystem.Load('libHiggsAnalysisCombinedLimit')
base=Path('../bias_phi_jpsi_quartic/m1.4')
files=[R.TFile.Open(str(base/(k+'_1.4.root'))) for k in ['C4','L4']]
w,l=[f.Get('w') for f in files]
w.loadSnapshot('background');l.loadSnapshot('background')
meta=json.loads((base/'L4.json').read_text())
rf=R.RooFit
def put(x,*opts):
 assert not getattr(w,'import')(x,*opts,rf.Silence())
cat=R.RooCategory('pdfindex','pdfindex');put(cat);cat=w.cat('pdfindex')
shared=','.join(v.GetName() for v in l.allVars() if not v.GetName().startswith('SR_'))
def import_l(obj):
 opts=[rf.RenameAllNodes('L4alt'),rf.RenameAllVariablesExcept('L4alt',shared),rf.RecycleConflictNodes()]
 put(obj,*opts)
 return w.arg(obj.GetName()+'_L4alt')
for stage,mcname in [('s','ModelConfig'),('b','ModelConfig_bonly')]:
 sim=R.RooSimultaneousOpt('envelope_'+stage,'',w.cat('CMS_channel'))
 keep=[]
 for ch in meta['channels']:
  label=ch['label'];pdfs=R.RooArgList();rates=R.RooArgList()
  for proc in ch['processes']:
   if stage=='b' and proc['signal']:continue
   shape=w.pdf(proc['shape']);rate=w.function(proc['coefficient'])
   if proc['name']=='continuum' and label.startswith('SR'):
    sn='env_'+shape.GetName();rn='env_'+rate.GetName()
    if not w.pdf(sn):
     ls=import_l(l.pdf(proc['shape']));lr=import_l(l.function(proc['coefficient']))
     multi=R.RooMultiPdf(sn,sn,cat,R.RooArgList(shape,ls));multi.setCorrectionFactor(0);put(multi,rf.RecycleConflictNodes())
     switch=R.RooFormulaVar(rn,'(@0==0)?@1:@2',R.RooArgList(cat,rate,lr));put(switch,rf.RecycleConflictNodes())
    shape=w.pdf(sn);rate=w.function(rn)
   pdfs.add(shape);rates.add(rate)
  p=R.RooAddPdf('env_'+stage+'_'+label,'',pdfs,rates);put(p,rf.RecycleConflictNodes());sim.addPdf(w.pdf(p.GetName()),label)
 for k in meta['constraints']:sim.addExtraConstraint(w.pdf(k))
 put(sim,rf.RecycleConflictNodes());mc=w.obj(mcname);mc.SetPdf(w.pdf(sim.GetName()))
 nuis=R.RooArgSet(mc.GetNuisanceParameters())
 for v in w.allVars():
  if (v.GetName().endswith('_L4alt') or v.GetName().startswith('logpoly')) and not v.isConstant():nuis.add(v)
 mc.SetNuisanceParameters(nuis)
assert cat.numTypes()==2
assert all(w.obj('ModelConfig').GetNuisanceParameters().find(v.GetName()) for v in w.allVars() if v.GetName().startswith('logpoly') and not v.isConstant())
cat.setConstant(False)
discrete=R.RooArgSet('discreteParams');discrete.add(cat)
getattr(w,'import')(discrete,'discreteParams')
assert w.genobj('discreteParams')
assert not w.obj('ModelConfig').GetNuisanceParameters().find('pdfindex')
assert not w.var('CMS_eff_trigger_e').isConstant()
# Compare each branch to its source using its own fitted nuisance values.
for idx,src in enumerate([w,l]):
 cat.setIndex(idx)
 for v in src.allVars():
  target=w.var(v.GetName()+'_L4alt') if idx and w.var(v.GetName()+'_L4alt') else w.var(v.GetName())
  if target and not v.isConstant():target.setVal(v.getVal())
 for ch in meta['channels']:
  lab=ch['label'];a=w.pdf('env_s_'+lab);b=src.pdf('pdf_s_'+lab);obs=R.RooArgSet(w.var('CMS_th1x'));other=R.RooArgSet(src.var('CMS_th1x'))
  assert abs(a.expectedEvents(obs)/b.expectedEvents(other)-1)<1e-8
  for i in range(len(ch['counts'])):
   w.var('CMS_th1x').setVal(i+.5);src.var('CMS_th1x').setVal(i+.5)
   assert abs(a.getVal(obs)-b.getVal(other))<1e-8
w.loadSnapshot('background');cat.setIndex(0);w.writeToFile('envelope.root')
