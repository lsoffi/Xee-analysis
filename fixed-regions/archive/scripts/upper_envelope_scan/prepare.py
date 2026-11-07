import sys,math,csv
MASS=sys.argv[1];SRDEG=int(sys.argv[2]);REGION=sys.argv[3];DEG=SRDEG if REGION=='signal' else 4;MODEL=f'chebyshev{DEG}'
import ROOT,pathlib,json,hashlib
ROOT.gROOT.SetBatch(True);ROOT.gSystem.Load('libHiggsAnalysisCombinedLimit')
P=pathlib.Path('/Users/imb/Documents/Codex/Analisi-statistica-Xee');R=P/'runs/20261030_upper_envelope_scan'/f'degree{SRDEG}'/('m'+MASS)/REGION;R.mkdir(parents=True,exist_ok=True)
src=(P/'data/mass_scan/261028/cards_region2_data_envelope_allCorrections_binned/ee'/MASS/'Xee_ee_4_allYears.root') if float(MASS)<=11 else (P/'runs/20261030_upper_envelope_scan/extrapolated_inputs'/('m'+MASS)/'Xee_ee_4_allYears.root')
win={'definition':'upper SR 10.5066666667-12'};LO=10.506666666666668;HI=12.;NB=112
if REGION=='control':LO=9.093333333333334;HI=10.506666666666668;NB=106
f=ROOT.TFile.Open(str(src));w=f.Get('w');x=w.var('mass');o=ROOT.RooArgSet(x);x.setRange('full_reference',6,12);x.setRange('local_window',LO,HI);w.var('r').setVal(1)
data=w.data('data_obs').reduce(ROOT.RooFit.Cut(f'mass>={LO:.16g} && mass<{HI:.16g}'));audit=[];rep=[];keep=[];constraints=[]
for i,era in enumerate(['2022','2022EE','2023','2023BPix'],1):
 ch=f'ch{i}';d=data.reduce(ROOT.RooFit.Cut(f'CMS_channel==CMS_channel::{ch}'));n=d.sumEntries()
 coeff=[ROOT.RooRealVar(f'local_c{j}_'+era,'',(-.3 if j==1 else .02 if j==2 else 0),-2 if j==1 else -1,2 if j==1 else 1) for j in range(1,DEG+1)]
 pdf=ROOT.RooChebychev(f'local_chebyshev{DEG}_'+era,'',x,ROOT.RooArgList(*coeff));keep.extend(coeff+[pdf]);getattr(w,'import')(pdf,ROOT.RooFit.RecycleConflictNodes());rep.append(w.pdf('shapeBkg_dy_'+ch).GetName()+'='+pdf.GetName())
 if REGION=='signal':
  b0=ROOT.RooConstVar('local_b0_'+era,'',1.)
  bs=[ROOT.RooRealVar(f'local_b{j}_'+era,'',1.,1e-6,1000.) for j in range(1,DEG+1)]
  bern=ROOT.RooBernstein(f'local_bernstein{DEG}_'+era,'',x,ROOT.RooArgList(b0,*bs));idx=ROOT.RooCategory('pdfindex_'+era,'')
  env=ROOT.RooMultiPdf('local_envelope_'+era,'',idx,ROOT.RooArgList(pdf,bern));env.setCorrectionFactor(.5)
  keep.extend([b0,bern,idx,env]+bs);getattr(w,'import')(env,ROOT.RooFit.RecycleConflictNodes());rep[-1]=w.pdf('shapeBkg_dy_'+ch).GetName()+'='+env.GetName()
 # Dense-grid positivity during fit; exact extrema check required after fitting before generation.
 checks=[]
 for it in range(257):
  t=-1+2*it/256;terms=[1.,t]
  for j in range(2,DEG+1):terms.append(2*t*terms[-1]-terms[-2])
  checks.append('(1'+''.join(f'+@{j-1}*({terms[j]:.17g})' for j in range(1,DEG+1))+'>1e-8)')
 constraint_args=ROOT.RooArgList(*[w.var(c.GetName()) for c in coeff]);expr='&&'.join(checks)
 constraint=ROOT.RooGenericPdf('positive_continuum_'+era,expr,constraint_args);getattr(w,'import')(constraint);constraints.append(constraint.GetName());keep.append(constraint)
 yn=ROOT.RooRealVar('local_yield_'+era,'',n,0,5*n);yn.setStringAttribute('combine.channel',ch);yn.setStringAttribute('combine.process','dy');getattr(w,'import')(yn);rep.append(f'n_exp_bin{ch}_proc_dy='+yn.GetName());keep.append(yn)
 era_audit=dict(era=era,events=n,processes={})
 for proc in ['Zd','upsilon1s','upsilon2s']:
  dcb=w.pdf(('shapeSig_' if proc=='Zd' else 'shapeBkg_')+proc+'_'+ch);mean=w.function(dcb.x0().GetName());sigma=w.function(dcb.sigmaL().GetName());assert dcb.sigmaL().GetName()==dcb.sigmaR().GetName()
  if proc!='Zd':
   # Fixed MC widths/tails; new shared free means, independent of signal nuisances.
   mn=ROOT.RooRealVar(f'Ymean_{proc}_{era}','',mean.getVal(),mean.getVal()-.2,mean.getVal()+.2)
   sg=ROOT.RooConstVar(f'Ysigma_{proc}_{era}','',sigma.getVal())
   tails=[ROOT.RooConstVar(f'Ytail_{proc}_{era}_{j}','',v.getVal()) for j,v in enumerate([dcb.alphaL(),dcb.nL(),dcb.alphaR(),dcb.nR()])]
   newpdf=ROOT.RooCrystalBall(f'local_shape_{proc}_{era}','',x,mn,sg,sg,*tails)
   keep.extend([mn,sg,newpdf]+tails);getattr(w,'import')(newpdf,ROOT.RooFit.RecycleConflictNodes())
   rep.append(dcb.GetName()+'='+newpdf.GetName());dcb=w.pdf(newpdf.GetName());mean=w.var(mn.GetName());sigma=w.function(sg.GetName())
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
  rate=w.function(f'extrapolated_rate_Zd_{era}') if proc=='Zd' and float(MASS)>11 else w.function(f'n_exp_bin{ch}_proc_{proc}');fullrate=rate.getVal();frac=acceptance.getVal();newrate=ROOT.RooFormulaVar('local_rate_'+proc+'_'+era,'@0*@1',ROOT.RooArgList(rate,acceptance));newrate.setAttribute('combine.signal',proc=='Zd');newrate.setStringAttribute('combine.channel',ch);newrate.setStringAttribute('combine.process',proc);getattr(w,'import')(newrate,ROOT.RooFit.RecycleConflictNodes());rep.append(rate.GetName()+'='+newrate.GetName());keep.extend([acceptance,newrate,pdf,yn])
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
(R/'audit.json').write_text(json.dumps(dict(source=str(src),source_sha256=hashlib.sha256(src.read_bytes()).hexdigest(),requested_window=win,range=[LO,HI],bins=NB,eras=audit,total_events=data.sumEntries(),definition=f'Upper SR envelope Cheb{SRDEG}/Bern{SRDEG}, CR4, signal in both; Upsilon means fixed after CR-only fit at seed stage'),indent=2))
print('DONE',MASS,LO,HI,NB,data.sumEntries())
