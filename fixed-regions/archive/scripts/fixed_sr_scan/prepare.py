import sys,math,csv
MASS=sys.argv[1];REGION=sys.argv[2];MODEL='chebyshev3'
import ROOT,pathlib,json,hashlib
ROOT.gROOT.SetBatch(True);ROOT.gSystem.Load('libHiggsAnalysisCombinedLimit')
P=pathlib.Path('/Users/imb/Documents/Codex/Analisi-statistica-Xee');R=P/'runs/20261030_fixed_sr_scan'/('m'+MASS)/REGION;R.mkdir(parents=True,exist_ok=True)
src=P/'data/mass_scan/261028/cards_region2_data_envelope_allCorrections_binned/ee'/MASS/'Xee_ee_4_allYears.root'
win={'definition':'fixed SR 6-9.0933333333 GeV'};LO=6.;HI=9.093333333333334;NB=232
if REGION=='control':LO=9.093333333333334;HI=10.506666666666668;NB=106
f=ROOT.TFile.Open(str(src));w=f.Get('w');x=w.var('mass');o=ROOT.RooArgSet(x);x.setRange('full_reference',6,12);x.setRange('local_window',LO,HI);w.var('r').setVal(1)
data=w.data('data_obs').reduce(ROOT.RooFit.Cut(f'mass>={LO:.16g} && mass<{HI:.16g}'));audit=[];rep=[];keep=[];constraints=[]
for i,era in enumerate(['2022','2022EE','2023','2023BPix'],1):
 ch=f'ch{i}';d=data.reduce(ROOT.RooFit.Cut(f'CMS_channel==CMS_channel::{ch}'));n=d.sumEntries()
 c1=ROOT.RooRealVar('local_c1_'+era,'',-.3,-2,2);c2=ROOT.RooRealVar('local_c2_'+era,'',.02,-1,1)
 c3=ROOT.RooRealVar('local_c3_'+era,'',0,-1,1)
 pdf=ROOT.RooChebychev('local_chebyshev3_'+era,'',x,ROOT.RooArgList(c1,c2,c3));keep.extend([c1,c2,c3,pdf])
 getattr(w,'import')(pdf,ROOT.RooFit.RecycleConflictNodes());rep.append(w.pdf('shapeBkg_dy_'+ch).GetName()+'='+pdf.GetName())
 # Exact positivity: endpoints and every stationary point inside [-1,1].
 def poly(t):return f'(1+@0*({t})+@1*(2*({t})*({t})-1)+@2*(4*({t})*({t})*({t})-3*({t})))'
 disc='(16*@1*@1-48*@2*(@0-3*@2))'
 roots=[f'((-4*@1{sgn}sqrt({disc}))/(24*@2))' for sgn in ['+','-']]
 checks=[f'((abs({t})>=1)?1:({poly(t)}>1e-8))' for t in roots]
 cubic=f'(({disc}<=0)?1:('+ '&&'.join(checks)+'))'
 quadratic='((@1<=0)?1:((abs(@0)>=4*@1)?1:(1-@1-@0*@0/(8*@1)>1e-8)))'
 expr=f'({poly(-1)}>1e-8)&&({poly(1)}>1e-8)&&((@2==0)?{quadratic}:{cubic})'
 constraint=ROOT.RooGenericPdf('positive_continuum_'+era,expr,ROOT.RooArgList(w.var(c1.GetName()),w.var(c2.GetName()),w.var(c3.GetName())));getattr(w,'import')(constraint);constraints.append(constraint.GetName());keep.append(constraint)
 yn=ROOT.RooRealVar('local_yield_'+era,'',n,0,5*n);yn.setStringAttribute('combine.channel',ch);yn.setStringAttribute('combine.process','dy');getattr(w,'import')(yn);rep.append(f'n_exp_bin{ch}_proc_dy='+yn.GetName());keep.append(yn)
 era_audit=dict(era=era,events=n,processes={})
 for proc in ['Zd','upsilon1s','upsilon2s']:
  dcb=w.pdf(('shapeSig_' if proc=='Zd' else 'shapeBkg_')+proc+'_'+ch);mean=w.function(dcb.x0().GetName());sigma=w.function(dcb.sigmaL().GetName());assert dcb.sigmaL().GetName()==dcb.sigmaR().GetName()
  params=ROOT.RooArgList(mean,sigma,dcb.alphaL(),dcb.nL(),dcb.alphaR(),dcb.nR())
  # Primitive of the unnormalized symmetric-width double Crystal Ball in t=(x-mu)/sigma.
  # Absolute alpha follows RooCrystalBall conventions; log branch also supports n=1.
  primitives=[]
  for edge in [6,LO,HI,12]:
   t=f'(({edge}-@0)/@1)';al='abs(@2)';ar='abs(@4)';nl='@3';nr='@5'
   gl=f'sqrt(TMath::Pi()/2)*TMath::Erf(-{al}/sqrt(2))';gr=f'sqrt(TMath::Pi()/2)*TMath::Erf({ar}/sqrt(2))'
   bl=f'({nl}/{al}-{al})';br=f'({nr}/{ar}-{ar})'
   left=f'pow({nl}/{al},{nl})*exp(-0.5*{al}*{al})*((abs({nl}-1)<1e-8)?log(({bl}+{al})/({bl}-{t})):(pow({bl}-{t},1-{nl})-pow({bl}+{al},1-{nl}))/({nl}-1))'
   right=f'pow({nr}/{ar},{nr})*exp(-0.5*{ar}*{ar})*((abs({nr}-1)<1e-8)?log(({br}+{t})/({br}+{ar})):(pow({br}+{t},1-{nr})-pow({br}+{ar},1-{nr}))/(1-{nr}))'
   expr=f'({t}<-{al})?({gl}+{left}):(({t}>{ar})?({gr}+{right}):(sqrt(TMath::Pi()/2)*TMath::Erf({t}/sqrt(2))))'
   primitive=ROOT.RooFormulaVar(f'dcb_primitive_{str(edge).replace(chr(46),chr(112))}_{proc}_{era}',expr,params);primitives.append(primitive);keep.append(primitive)
  acceptance=ROOT.RooFormulaVar('local_acceptance_'+proc+'_'+era,'(@2-@1)/(@3-@0)',ROOT.RooArgList(*primitives))
  # Independent check against ROOT integration before cropping; also vary the original nuisance parameters.
  checks=[]
  for nuisance in [None,'CMS_scale_e',f'CMS_EXO25020_signalModelSigmaNuisance_{era}'] if proc=='Zd' else [None,'CMS_scale_e']:
   for shift in ([0] if nuisance is None else [-1,1]):
    v=w.var(nuisance) if nuisance else None;oldval=v.getVal() if v else 0
    if v:v.setVal(shift)
    integ=dcb.createIntegral(o,ROOT.RooFit.NormSet(o),ROOT.RooFit.Range('local_window'));truth=integ.getVal();calc=acceptance.getVal();assert abs(calc-truth)<1e-8,(MASS,era,nuisance,shift,calc,truth)
    checks.append(dict(nuisance=nuisance,shift=shift,analytic=calc,root_integral=truth))
    if v:v.setVal(oldval)
  rate=w.function(f'n_exp_bin{ch}_proc_{proc}');fullrate=rate.getVal();frac=acceptance.getVal();newrate=ROOT.RooFormulaVar('local_rate_'+proc+'_'+era,'@0*@1',ROOT.RooArgList(rate,acceptance));newrate.setAttribute('combine.signal',proc=='Zd');newrate.setStringAttribute('combine.channel',ch);newrate.setStringAttribute('combine.process',proc);getattr(w,'import')(newrate,ROOT.RooFit.RecycleConflictNodes());rep.append(rate.GetName()+'='+newrate.GetName());keep.extend([acceptance,newrate,pdf,yn])
  era_audit['processes'][proc]=dict(full_rate_r1=fullrate,fraction=frac,local_rate_r1=fullrate*frac,acceptance_checks=checks)
  if proc=='Zd':assert 0 < frac <= 1
 ysum=sum(a['local_rate_r1'] for proc,a in era_audit['processes'].items() if proc!='Zd');w.var(yn.GetName()).setVal(max(1,n-ysum));audit.append(era_audit)
x.setRange(LO,HI);x.setBins(NB);x.setVal((LO+HI)/2);nw=ROOT.RooWorkspace('w')
for old,mcname in [('model_s','ModelConfig'),('model_b','ModelConfig_bonly')]:
 reps=rep if mcname=='ModelConfig' else [p for p in rep if '_proc_Zd=' not in p]
 pdf=w.factory('EDIT::'+old+'_local('+old+','+','.join(reps)+')');assert pdf
 for name in constraints:pdf.addExtraConstraint(w.pdf(name))
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
nw.var('r').setVal(1);nw.var('r').setRange(-100,100);nw.var('r').setConstant(False)
for a in audit:
 for proc,vals in a['processes'].items():
  value=nw.function('local_rate_'+proc+'_'+a['era']).getVal();assert abs(value-vals['local_rate_r1'])<1e-7*max(1,abs(value))
nw.var('mass').setVal((LO+HI)/2);nw.var('r').setVal(0);nw.writeToFile(str(R/'workspace.root'))
(R/'audit.json').write_text(json.dumps(dict(source=str(src),source_sha256=hashlib.sha256(src.read_bytes()).hexdigest(),requested_window=win,range=[LO,HI],bins=NB,eras=audit,total_events=data.sumEntries(),definition='Degree3 Chebyshev continuum + original DCB signal + explicit original Upsilon1S and Upsilon2S; all original Y normalizations/nuisances floating; continuum constrained to physical positive support'),indent=2))
print('DONE',MASS,LO,HI,NB,data.sumEntries())
