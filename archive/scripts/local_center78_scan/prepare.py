import sys
MASS=sys.argv[1]
MODEL=sys.argv[2]
import ROOT,pathlib,json,hashlib
ROOT.gROOT.SetBatch(True);ROOT.gSystem.Load('libHiggsAnalysisCombinedLimit')
P=pathlib.Path('/Users/imb/Documents/Codex/Analisi-statistica-Xee');R=P/'runs/20261029_local_center78_scan'/MODEL/('m'+MASS);R.mkdir(parents=True,exist_ok=True);src=P/'runs/20261029_gaussian_signal_scan/workspaces/dcb'/('m'+MASS+'.root')
f=ROOT.TFile.Open(str(src));w=f.Get('w');x=w.var('mass');o=ROOT.RooArgSet(x);x.setRange('full_reference',6,12);x.setRange('local_window',7.293333333333333,8.306666666666667);w.var('r').setVal(1)
data=w.data('data_obs').reduce(ROOT.RooFit.Cut('mass>=7.293333333333333 && mass<8.306666666666667'));audit=[];rep=[];keep=[]
for i,era in enumerate(['2022','2022EE','2023','2023BPix'],1):
 ch=f'ch{i}';d=data.reduce(ROOT.RooFit.Cut(f'CMS_channel==CMS_channel::{ch}'));n=d.sumEntries()
 if MODEL=='exp':
  slope=ROOT.RooRealVar('local_slope_'+era,'',-0.6,-10,10);pdf=ROOT.RooExponential('local_exp_'+era,'',x,slope);keep.append(slope)
 elif MODEL=='bernstein2':
  c0=ROOT.RooConstVar('local_b0_'+era,'',1)
  c1=ROOT.RooRealVar('local_b1_'+era,'',.75,1e-6,20)
  c2=ROOT.RooRealVar('local_b2_'+era,'',.55,1e-6,20)
  pdf=ROOT.RooBernstein('local_bernstein2_'+era,'',x,ROOT.RooArgList(c0,c1,c2));keep.extend([c0,c1,c2])
 elif MODEL=='chebyshev2':
  c1=ROOT.RooRealVar('local_c1_'+era,'',-.28,-2,2)
  c2=ROOT.RooRealVar('local_c2_'+era,'',.02,-1,1)
  pdf=ROOT.RooChebychev('local_chebyshev2_'+era,'',x,ROOT.RooArgList(c1,c2));keep.extend([c1,c2])
 else:raise ValueError(MODEL)
 getattr(w,'import')(pdf,ROOT.RooFit.RecycleConflictNodes());rep.append(w.pdf('shapeBkg_dy_'+ch).getPdf(1).GetName()+'='+pdf.GetName())
 yn=ROOT.RooRealVar('local_yield_'+era,'',n,0,5*n);yn.setStringAttribute('combine.channel',ch);yn.setStringAttribute('combine.process','dy');getattr(w,'import')(yn);rep.append(f'n_exp_bin{ch}_proc_dy='+yn.GetName())
 # The exponential is the total smooth local background, including resonance tails.
 for p in ['upsilon1s','upsilon2s']:
  z=ROOT.RooConstVar('local_zero_'+p+'_'+era,'',0);z.setStringAttribute('combine.channel',ch);z.setStringAttribute('combine.process',p);getattr(w,'import')(z);rep.append(f'n_exp_bin{ch}_proc_{p}='+z.GetName());keep.append(z)
 dcb=w.pdf('shapeSig_Zd_'+ch);mean=w.function(dcb.x0().GetName());sigma=w.function(dcb.sigmaL().GetName());assert dcb.sigmaL().GetName()==dcb.sigmaR().GetName()
 params=ROOT.RooArgList(mean,sigma,dcb.alphaL(),dcb.nL(),dcb.alphaR(),dcb.nR())
 # Primitive of the unnormalized symmetric-width double Crystal Ball in t=(x-mu)/sigma.
 # Absolute alpha follows RooCrystalBall conventions; log branch also supports n=1.
 primitives=[]
 for edge in [6,7.293333333333333,8.306666666666667,12]:
  t=f'(({edge}-@0)/@1)';al='abs(@2)';ar='abs(@4)';nl='@3';nr='@5'
  gl=f'sqrt(TMath::Pi()/2)*TMath::Erf(-{al}/sqrt(2))';gr=f'sqrt(TMath::Pi()/2)*TMath::Erf({ar}/sqrt(2))'
  bl=f'({nl}/{al}-{al})';br=f'({nr}/{ar}-{ar})'
  left=f'pow({nl}/{al},{nl})*exp(-0.5*{al}*{al})*((abs({nl}-1)<1e-8)?log(({bl}+{al})/({bl}-{t})):(pow({bl}-{t},1-{nl})-pow({bl}+{al},1-{nl}))/({nl}-1))'
  right=f'pow({nr}/{ar},{nr})*exp(-0.5*{ar}*{ar})*((abs({nr}-1)<1e-8)?log(({br}+{t})/({br}+{ar})):(pow({br}+{t},1-{nr})-pow({br}+{ar},1-{nr}))/(1-{nr}))'
  expr=f'({t}<-{al})?({gl}+{left}):(({t}>{ar})?({gr}+{right}):(sqrt(TMath::Pi()/2)*TMath::Erf({t}/sqrt(2))))'
  primitive=ROOT.RooFormulaVar(f'dcb_primitive_{str(edge).replace(chr(46),chr(112))}_{era}',expr,params);primitives.append(primitive);keep.append(primitive)
 acceptance=ROOT.RooFormulaVar('local_signal_acceptance_'+era,'(@2-@1)/(@3-@0)',ROOT.RooArgList(*primitives))
 # Independent check against ROOT integration before cropping; also vary the original nuisance parameters.
 checks=[]
 for nuisance in [None,'CMS_scale_e',f'CMS_EXO25020_signalModelSigmaNuisance_{era}']:
  for shift in ([0] if nuisance is None else [-1,1]):
   v=w.var(nuisance) if nuisance else None;oldval=v.getVal() if v else 0
   if v:v.setVal(shift)
   integ=dcb.createIntegral(o,ROOT.RooFit.NormSet(o),ROOT.RooFit.Range('local_window'));truth=integ.getVal();calc=acceptance.getVal();assert abs(calc-truth)<1e-8,(MASS,era,nuisance,shift,calc,truth)
   checks.append(dict(nuisance=nuisance,shift=shift,analytic=calc,root_integral=truth))
   if v:v.setVal(oldval)
 rate=w.function(f'n_exp_bin{ch}_proc_Zd');fullrate=rate.getVal();frac=acceptance.getVal();newrate=ROOT.RooFormulaVar('local_signal_rate_'+era,'@0*@1',ROOT.RooArgList(rate,acceptance));newrate.setAttribute('combine.signal',True);newrate.setStringAttribute('combine.channel',ch);newrate.setStringAttribute('combine.process','Zd');getattr(w,'import')(newrate,ROOT.RooFit.RecycleConflictNodes());rep.append(rate.GetName()+'='+newrate.GetName());keep.extend([acceptance,newrate,pdf,yn])
 audit.append(dict(acceptance_checks=checks,era=era,events=n,signal_full_r1=fullrate,signal_fraction=frac,signal_local_r1=fullrate*frac,mean=w.pdf("shapeSig_Zd_"+ch).x0().getVal(),sigma=w.pdf("shapeSig_Zd_"+ch).sigmaL().getVal()))
# Crop both the data and the observable; keep named full range only for signal acceptance.
x.setRange(7.293333333333333,8.306666666666667);x.setBins(76);x.setVal(7.8);nw=ROOT.RooWorkspace('w')
for old,mcname in [('model_s_study','ModelConfig'),('model_b_study','ModelConfig_bonly')]:
 reps=rep if mcname=='ModelConfig' else [p for p in rep if '_proc_Zd=' not in p]
 pdf=w.factory('EDIT::'+old+'_local('+old+','+','.join(reps)+')');assert pdf
 mc=w.obj(mcname);getattr(nw,'import')(pdf,ROOT.RooFit.RecycleConflictNodes())
 newmc=ROOT.RooStats.ModelConfig(mcname,nw);newmc.SetPdf(nw.pdf(pdf.GetName()))
 for getter,setter in [('GetParametersOfInterest','SetParametersOfInterest'),('GetObservables','SetObservables'),('GetNuisanceParameters','SetNuisanceParameters'),('GetGlobalObservables','SetGlobalObservables')]:
  nodes=ROOT.RooArgSet()
  for v in getattr(mc,getter)():
   target=nw.arg(v.GetName())
   if target:nodes.add(target)
  getattr(newmc,setter)(nodes)
 getattr(nw,'import')(newmc);keep.append(newmc)
getattr(nw,'import')(data,ROOT.RooFit.Rename('data_obs'))
nw.var('r').setVal(0);nw.var('r').setRange(-5,5);nw.var('r').setConstant(False)
for a in audit:
 era=a['era'];nw.var('r').setVal(1);a['local_rate_after_range_change']=nw.function('local_signal_rate_'+era).getVal();assert abs(a['local_rate_after_range_change']/a['signal_local_r1']-1)<1e-9
nw.var('mass').setVal(7.8);nw.var('r').setVal(0);nw.writeToFile(str(R/'workspace.root'))
(R/'audit.json').write_text(json.dumps(dict(source=str(src),source_sha256=hashlib.sha256(src.read_bytes()).hexdigest(),range=[7.293333333333333,8.306666666666667],bins=76,eras=audit,total_events=data.sumEntries(),definition=MODEL+' per era for total local smooth background, including Upsilon tails'),indent=2))
print('DONE',data.sumEntries())
