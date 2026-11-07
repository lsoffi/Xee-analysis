import sys
MASS=sys.argv[1]
import ROOT,pathlib,json,hashlib
ROOT.gROOT.SetBatch(True);ROOT.gSystem.Load('libHiggsAnalysisCombinedLimit')
P=pathlib.Path('/Users/imb/Documents/Codex/Analisi-statistica-Xee');R=P/'runs/20261029_local_gaussian_scan_7_8'/('m'+MASS);R.mkdir(parents=True,exist_ok=True);src=P/'runs/20261029_gaussian_signal_scan/workspaces/gaussian'/('m'+MASS+'.root')
f=ROOT.TFile.Open(str(src));w=f.Get('w');x=w.var('mass');o=ROOT.RooArgSet(x);x.setRange('full_reference',6,12);x.setRange('local_window',7,8);w.var('r').setVal(1)
data=w.data('data_obs').reduce(ROOT.RooFit.Cut('mass>=7 && mass<8'));audit=[];rep=[];keep=[]
for i,era in enumerate(['2022','2022EE','2023','2023BPix'],1):
 ch=f'ch{i}';d=data.reduce(ROOT.RooFit.Cut(f'CMS_channel==CMS_channel::{ch}'));n=d.sumEntries()
 slope=ROOT.RooRealVar('local_slope_'+era,'',-1,-10,10);pdf=ROOT.RooExponential('local_exp_'+era,'',x,slope);getattr(w,'import')(pdf,ROOT.RooFit.RecycleConflictNodes());rep.append(w.pdf('shapeBkg_dy_'+ch).getPdf(1).GetName()+'='+pdf.GetName())
 yn=ROOT.RooRealVar('local_yield_'+era,'',n,0,5*n);yn.setStringAttribute('combine.channel',ch);yn.setStringAttribute('combine.process','dy');getattr(w,'import')(yn);rep.append(f'n_exp_bin{ch}_proc_dy='+yn.GetName())
 # The exponential is the total smooth local background, including resonance tails.
 for p in ['upsilon1s','upsilon2s']:
  z=ROOT.RooConstVar('local_zero_'+p+'_'+era,'',0);z.setStringAttribute('combine.channel',ch);z.setStringAttribute('combine.process',p);getattr(w,'import')(z);rep.append(f'n_exp_bin{ch}_proc_{p}='+z.GetName());keep.append(z)
 sig=w.pdf('gaussian_signal_'+era);dcb=w.pdf('shapeSig_Zd_'+ch);mean=w.function(dcb.x0().GetName());sigma=w.function(dcb.sigmaL().GetName())
 acceptance=ROOT.RooFormulaVar('local_signal_acceptance_'+era,'(TMath::Erf((8-@0)/(sqrt(2)*@1))-TMath::Erf((7-@0)/(sqrt(2)*@1)))/(TMath::Erf((12-@0)/(sqrt(2)*@1))-TMath::Erf((6-@0)/(sqrt(2)*@1)))',ROOT.RooArgList(mean,sigma))
 rate=w.function(f'n_exp_bin{ch}_proc_Zd');fullrate=rate.getVal();frac=acceptance.getVal();newrate=ROOT.RooFormulaVar('local_signal_rate_'+era,'@0*@1',ROOT.RooArgList(rate,acceptance));newrate.setAttribute('combine.signal',True);newrate.setStringAttribute('combine.channel',ch);newrate.setStringAttribute('combine.process','Zd');getattr(w,'import')(newrate,ROOT.RooFit.RecycleConflictNodes());rep.append(rate.GetName()+'='+newrate.GetName());keep.extend([acceptance,newrate,pdf,slope,yn])
 audit.append(dict(era=era,events=n,signal_full_r1=fullrate,signal_fraction=frac,signal_local_r1=fullrate*frac,mean=w.pdf("shapeSig_Zd_"+ch).x0().getVal(),sigma=w.pdf("shapeSig_Zd_"+ch).sigmaL().getVal()))
# Crop both the data and the observable; keep named full range only for signal acceptance.
x.setRange(7,8);x.setBins(75);nw=ROOT.RooWorkspace('w')
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
nw.var('r').setVal(0);nw.writeToFile(str(R/'workspace.root'))
(R/'audit.json').write_text(json.dumps(dict(source=str(src),source_sha256=hashlib.sha256(src.read_bytes()).hexdigest(),range=[7,8],bins=75,eras=audit,total_events=data.sumEntries(),definition='one exponential per era for total local smooth background, including Upsilon tails'),indent=2))
print('DONE',data.sumEntries())
