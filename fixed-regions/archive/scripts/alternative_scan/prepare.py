import argparse,json,pathlib,hashlib,math
import ROOT,numpy as np
ROOT.gROOT.SetBatch(True);ROOT.gSystem.Load('libHiggsAnalysisCombinedLimit')
P=pathlib.Path('/Users/imb/Documents/Codex/Analisi-statistica-Xee');R=P/'runs/20261029_alternative_background_scan';SRC=P/'data/mass_scan/261028/cards_region2_data_envelope_allCorrections_binned/ee';ERAS=['2022','2022EE','2023','2023BPix'];FAMILIES=['logquad','expquad','doubleexp','bernstein7'];ap=argparse.ArgumentParser();ap.add_argument('--stage',choices=['seed','scan'],required=True);args=ap.parse_args();audit=[];transfers=[]

def create(mass,family,selection):
 src=SRC/mass/f'Xee_ee_4_{selection}.root';f=ROOT.TFile.Open(str(src));w=f.Get('w');repl=[];newnames=[];keep=[];obs=ROOT.RooArgSet(w.var('mass'))
 nuisance_names=[v.GetName() for v in w.obj('ModelConfig').GetNuisanceParameters()]
 pdfs=[p for p in w.allPdfs() if p.GetName().startswith('shapeBkg_dy_')]
 assert len(pdfs)==(4 if selection=='allYears' else 1)
 for old in pdfs:
  ch=old.GetName().removeprefix('shapeBkg_dy_');era=ERAS[int(ch[2:])-1] if selection=='allYears' else selection;prefix=f'ALT_{family}_{era}';x=w.var('mass');ref=old.getPdf(0)
  grid=np.linspace(6.001,11.999,300);ys=[]
  for xx in grid:x.setVal(float(xx));ys.append(ref.getVal(obs))
  ys=np.maximum(ys,1e-12);ys=ys/np.trapezoid(ys,grid);u=(grid-9)/3;ell=np.log(grid/9)
  def var(s,val,lo,hi):
   v=ROOT.RooRealVar(prefix+'_'+s,'',max(lo+1e-6,min(hi-1e-6,float(val))),lo,hi);keep.append(v);newnames.append(v.GetName());return v
  if family in ['logquad','expquad']:
   coord=ell if family=='logquad' else u;c=np.linalg.lstsq(np.c_[np.ones(len(grid)),coord,coord**2],np.log(ys),rcond=None)[0]
   a=var('c1',c[1],-100 if family=='logquad' else -30,50 if family=='logquad' else 30);b=var('c2',c[2],-200 if family=='logquad' else -30,200 if family=='logquad' else 30)
   term='log(@2/9.)' if family=='logquad' else '((@2-9.)/3.)'
   pdf=ROOT.RooGenericPdf(prefix+'_pdf','',f'exp(@0*{term}+@1*pow({term},2))',ROOT.RooArgList(a,b,x))
  elif family=='doubleexp':
   best=(1e300,None)
   def expnorm(k):return np.exp(-k*(grid-6)/6)*k/(6*(-np.expm1(-k)))
   for slow in [.2,1,2,3,4,6,8]:
    for fast in [4,6,8,12,18,24,36,48]:
     if fast<=slow:continue
     for mix in [.05,.2,.5,.8,.95]:
      yy=mix*expnorm(slow)+(1-mix)*expnorm(fast);score=np.mean((np.log(yy)-np.log(ys))**2)
      if score<best[0]:best=(score,(slow,fast,mix))
   slow,fast,mix=best[1];k=var('kslow',slow,.001,30);d=var('gap',fast-slow,.001,60);frac=var('fraction',mix,0,1)
   cslow=ROOT.RooFormulaVar(prefix+'_slope_slow','-@0/6.',ROOT.RooArgList(k));cfast=ROOT.RooFormulaVar(prefix+'_slope_fast','-(@0+@1)/6.',ROOT.RooArgList(k,d));e1=ROOT.RooExponential(prefix+'_slow','',x,cslow);e2=ROOT.RooExponential(prefix+'_fast','',x,cfast);pdf=ROOT.RooAddPdf(prefix+'_pdf','',ROOT.RooArgList(e1,e2),ROOT.RooArgList(frac));keep += [cslow,cfast,e1,e2]
  else:
   a0=w.var(f'a0_{era}_envelope').getVal();assert a0>0;one=ROOT.RooConstVar(prefix+'_one','',1);coeff=ROOT.RooArgList(one);keep.append(one)
   for j in range(1,8):coeff.add(var('ratio'+str(j),w.var(f'a{j}_{era}_envelope').getVal()/a0,-100,100))
   pdf=ROOT.RooBernstein(prefix+'_pdf','',x,coeff);keep.append(coeff)
  keep.append(pdf);getattr(w,'import')(pdf);repl.append(f'{old.GetName()}={pdf.GetName()}')
  # The seed curves must be finite, positive and normalizable.
  candidate=w.pdf(pdf.GetName());v=[]
  for xx in grid:x.setVal(float(xx));v.append(candidate.getVal(obs))
  assert np.isfinite(v).all() and min(v)>0,(mass,family,era,min(v))
 for label,mcname in [('model_s','ModelConfig'),('model_b','ModelConfig_bonly')]:
  new=w.factory(f'EDIT::{label}_alternative({label},'+','.join(repl)+')');assert new;mc=w.obj(mcname);mc.SetPdf(new);nuis=ROOT.RooArgSet()
  for name in nuisance_names+newnames:nuis.add(w.var(name))
  mc.SetNuisanceParameters(nuis)
 discrete=w.genobj('discreteParams')
 if discrete:discrete.removeAll()
 assert w.var('mass').getMin()==6 and w.var('mass').getMax()==12 and w.var('mass').getBins()==450
 w.var('r').setRange(-100,100);w.var('r').setVal(0);w.var('r').setError(.1);w.var('r').setConstant(False)
 if args.stage=='scan':
  for era in ERAS:
   sf=ROOT.TFile.Open(str(R/f'seeds/{family}/{era}/fitDiagnostics_seed.root'));fit=sf.Get('fit_b');assert fit and fit.status()==0,(family,era)
   for v in fit.floatParsFinal():
    name=v.GetName();target=w.var(name)
    if target and (name.startswith(f'ALT_{family}_{era}_') or name.startswith('CMS_EXO25020_bkgScale_') and name.endswith('_'+era)):
     target.setVal(v.getVal());target.setError(v.getError());transfers.append(dict(mass=mass,family=family,era=era,parameter=name,source_value=v.getVal(),seed_value=target.getVal(),frozen=target.isConstant()))
   sf.Close()
  w.saveSnapshot('background_initialized',w.allVars())
 dest=R/(f'workspaces/seed/{family}/{selection}.root' if args.stage=='seed' else f'workspaces/{family}/m{mass}.root');dest.parent.mkdir(parents=True,exist_ok=True);w.writeToFile(str(dest))
 active=w.obj('ModelConfig').GetPdf().getParameters(w.obj('ModelConfig').GetObservables());names=[v.GetName() for v in active];assert sum(n.startswith('ALT_') for n in names)==len(newnames)
 audit.append(dict(mass=float(mass),mass_label=mass,family=family,selection=selection,path=str(dest),source=str(src),source_sha256=hashlib.sha256(src.read_bytes()).hexdigest(),data_entries=w.data('data_obs').numEntries(),events=w.data('data_obs').sumEntries(),new_shape_parameters=newnames,range=[6,12],bins_per_era=450))
 f.Close()

if args.stage=='seed':
 for family in FAMILIES:
  for era in ERAS:create('9.0',family,era)
else:
 for family in FAMILIES:
  for mass in sorted([p.name for p in SRC.iterdir() if p.is_dir() and (p/'Xee_ee_4_allYears.root').exists()],key=float):create(mass,family,'allYears')
(R/f'audit_{args.stage}.json').write_text(json.dumps(audit,indent=2))
if transfers:(R/'background_seed_transfer.json').write_text(json.dumps(transfers,indent=2))
print('PREPARED',len(audit),'workspaces',flush=True)
