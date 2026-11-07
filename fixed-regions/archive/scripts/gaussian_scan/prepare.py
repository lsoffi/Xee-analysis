import ROOT,pathlib,json,argparse,hashlib
import numpy as np
ROOT.gROOT.SetBatch(True);ROOT.gSystem.Load('libHiggsAnalysisCombinedLimit')
P=pathlib.Path('/Users/imb/Documents/Codex/Analisi-statistica-Xee');R=P/'runs/20261029_gaussian_signal_scan';SRC=P/'data/mass_scan/261028/cards_region2_data_envelope_allCorrections_binned/ee';ERAS=['2022','2022EE','2023','2023BPix'];ap=argparse.ArgumentParser();ap.add_argument('--stage',choices=['seed','scan'],required=True);args=ap.parse_args();audit=[];shapes=[];transfer=[]
def create(mass,model,selection):
 src=SRC/mass/f'Xee_ee_4_{selection}.root';f=ROOT.TFile.Open(str(src));w=f.Get('w');replace=[];keep=[];channels=['ch'+str(i) for i in range(1,5)] if selection=='allYears' else ['Xee_ee_4_'+selection];obs=ROOT.RooArgSet(w.var('mass'))
 for ch in channels:
  era=ERAS[int(ch[2:])-1] if selection=='allYears' else selection;multi=w.pdf('shapeBkg_dy_'+ch);assert multi and multi.ClassName()=='RooMultiPdf';cheb=multi.getPdf(1);assert cheb.ClassName()=='RooChebychev';replace.append(multi.GetName()+'='+cheb.GetName());coeffs=[v for v in cheb.getParameters(obs) if v.GetName().startswith('t')];assert len(coeffs)==7
  if selection=='allYears':
   dcb=w.pdf('shapeSig_Zd_'+ch);assert dcb.ClassName()=='RooCrystalBall';mean=w.function(dcb.x0().GetName());sigma=w.function(dcb.sigmaL().GetName());assert dcb.sigmaL().GetName()==dcb.sigmaR().GetName();assert mean and sigma
   gauss=ROOT.RooGaussian('gaussian_signal_'+era,'',w.var('mass'),mean,sigma);getattr(w,'import')(gauss);keep.append(gauss)
   if model=='gaussian':replace.append(dcb.GetName()+'='+gauss.GetName())
   rate=w.function(f'n_exp_bin{ch}_proc_Zd');w.var('r').setVal(1);shapes.append(dict(mass=float(mass),model=model,era=era,mean_nominal=mean.getVal(),sigma_nominal=sigma.getVal(),mean_node=mean.GetName(),sigma_node=sigma.GetName(),signal_events_at_r1=rate.getVal()))
 for label,mcname in [('model_s','ModelConfig'),('model_b','ModelConfig_bonly')]:
  # The background-only graph does not contain the signal PDF; replacing it is unnecessary.
  reps=replace if label=='model_s' else [p for p in replace if not p.startswith('shapeSig_')]
  new=w.factory(f'EDIT::{label}_study({label},'+','.join(reps)+')');assert new;w.obj(mcname).SetPdf(new)
 d=w.genobj('discreteParams')
 if d:d.removeAll()
 w.var('r').setRange(-100,100);w.var('r').setVal(0);w.var('r').setError(.1);w.var('r').setConstant(False)
 if args.stage=='scan':
  for era in ERAS:
   sf=ROOT.TFile.Open(str(R/f'seeds/{era}/fitDiagnostics_seed.root'));fit=sf.Get('fit_b');assert fit.status()==0
   for v in fit.floatParsFinal():
    n=v.GetName();target=w.var(n)
    if target and (n.startswith('t') and n.endswith('_'+era+'_envelope') or n.startswith('CMS_EXO25020_bkgScale_') and n.endswith('_'+era)):
     target.setVal(v.getVal());target.setError(v.getError());transfer.append(dict(mass=mass,model=model,era=era,parameter=n,source_value=v.getVal(),seed_value=target.getVal(),frozen=target.isConstant()))
   sf.Close()
  w.saveSnapshot('background_initialized',w.allVars())
 dest=R/(f'workspaces/seed/{selection}.root' if args.stage=='seed' else f'workspaces/{model}/m{mass}.root');dest.parent.mkdir(parents=True,exist_ok=True);w.writeToFile(str(dest))
 assert w.var('mass').getMin()==6 and w.var('mass').getMax()==12 and w.var('mass').getBins()==450
 active=w.obj('ModelConfig').GetPdf().getParameters(w.obj('ModelConfig').GetObservables());assert sum(v.GetName().startswith('t') and v.GetName().endswith('_envelope') for v in active)==7*len(channels)
 audit.append(dict(mass=float(mass),mass_label=mass,model=model,selection=selection,path=str(dest),source=str(src),source_sha256=hashlib.sha256(src.read_bytes()).hexdigest(),events=w.data('data_obs').sumEntries(),entries=w.data('data_obs').numEntries(),range=[6,12],bins_per_era=450));f.Close()
if args.stage=='seed':
 for era in ERAS:create('7.5','dcb',era)
else:
 for model in ['dcb','gaussian']:
  for mass in sorted([p.name for p in SRC.iterdir() if p.is_dir() and (p/'Xee_ee_4_allYears.root').exists()],key=float):create(mass,model,'allYears')
(R/f'audit_{args.stage}.json').write_text(json.dumps(audit,indent=2))
if shapes:(R/'signal_shape_definitions.json').write_text(json.dumps(shapes,indent=2))
if transfer:(R/'seed_transfer.json').write_text(json.dumps(transfer,indent=2))
print('PREPARED',len(audit),'workspaces',flush=True)
