from pathlib import Path
import ROOT,json,math,hashlib
ROOT.gROOT.SetBatch(True);ROOT.gSystem.Load('libHiggsAnalysisCombinedLimit');P=Path.cwd();R=P/'runs/20261030_upsilon_mu_sigma_cheb3';pars=json.loads((P/'runs/20261029_upsilon_peak_windows/parameters.json').read_text())['upsilon'];eras=['2022','2022EE','2023','2023BPix'];src=P/'data/mass_scan/261028/cards_region2_data_envelope_allCorrections_binned/ee/7.5/Xee_ee_4_allYears.root';f=ROOT.TFile.Open(str(src));old=f.Get('w');configs=[]
for mult in [1.5,2.,3.]:
 ones=[p for p in pars if p['process']=='upsilon1s'];lo=min(p['mean']-mult*p['sigmaL'] for p in ones);hi=max(p['mean']+mult*p['sigmaR'] for p in ones);LO=6+math.floor((lo-6)*75+1e-8)/75;HI=6+math.ceil((hi-6)*75-1e-8)/75;NB=round((HI-LO)*75)
 for i,era in enumerate(eras,1):
  yp={p['process']:p for p in pars if p['era']==era};source=old.data('data_obs').reduce(ROOT.RooFit.Cut(f'CMS_channel==CMS_channel::ch{i} && mass>={LO!r} && mass<{HI!r}'));N=source.sumEntries()
  for degree in [3]:
   tag=f'k{mult:g}_{era}_b{degree}';d=R/tag;d.mkdir(exist_ok=True);w=ROOT.RooWorkspace('w');x=ROOT.RooRealVar('mass','mass',yp['upsilon1s']['mean'],LO,HI);x.setBins(NB);data=ROOT.RooDataHist('data_obs','',ROOT.RooArgList(x));counts=[]
   for j in range(source.numEntries()):
    row=source.get(j);x.setVal(row.getRealValue('mass'));data.add(ROOT.RooArgSet(x),source.weight());counts.append(dict(mass=x.getVal(),events=source.weight()))
   assert abs(data.sumEntries()-N)<1e-8;assert data.numEntries()==NB
   mu=ROOT.RooRealVar('upsilon_mu','Upsilon1S mean',yp['upsilon1s']['mean'],9.0,9.9);sigma=ROOT.RooRealVar('upsilon_sigma','Upsilon1S core sigma',yp['upsilon1s']['sigmaL'],.02,.6);keep=[x,mu,sigma];pdfs=[]
   for proc in ['upsilon1s','upsilon2s']:
    r=yp[proc];cs={key:ROOT.RooConstVar(proc+'_'+key,'',r[key]) for key in ['mean','sigmaL','sigmaR','alphaL','nL','alphaR','nR']};keep.extend(cs.values());pdf=ROOT.RooCrystalBall(proc+'_shape','',x,mu if proc=='upsilon1s' else cs['mean'],sigma if proc=='upsilon1s' else cs['sigmaL'],sigma if proc=='upsilon1s' else cs['sigmaR'],cs['alphaL'],cs['nL'],cs['alphaR'],cs['nR']);sampled=ROOT.RooBinSamplingPdf(proc+'_bin','',x,pdf,1e-7);pdfs.append(sampled);keep.extend([pdf,sampled])
   oldx=old.var('mass');oldx.setRange('local_measurement',LO,HI);ypdf=old.pdf(f'shapeBkg_upsilon2s_ch{i}');integ=ypdf.createIntegral(ROOT.RooArgSet(oldx),ROOT.RooFit.NormSet(ROOT.RooArgSet(oldx)),ROOT.RooFit.Range('local_measurement'));nominal2=old.function(f'n_exp_binch{i}_proc_upsilon2s').getVal()*integ.getVal();scale2=ROOT.RooRealVar('scale_upsilon2s','',1,0,10);scale2.setError(.2);const2=ROOT.RooConstVar('nominal_upsilon2s','',nominal2);n2=ROOT.RooFormulaVar('N_upsilon2s','@0*@1',ROOT.RooArgList(scale2,const2));n1=ROOT.RooRealVar('N_upsilon1s','',.45*N,0,3*N);n1.setError(math.sqrt(N));norms=[n1,n2];nuisnames=['N_upsilon1s','scale_upsilon2s'];keep.extend([n1,n2,scale2,const2]);mu.setError(.005);sigma.setError(.01)
   coeff=[ROOT.RooRealVar('cheb_c1','',-.3,-2,2),ROOT.RooRealVar('cheb_c2','',.02,-1,1),ROOT.RooRealVar('cheb_c3','',0,-1,1)]
   for c in coeff:c.setError(.02)
   cont=ROOT.RooChebychev('continuum_shape','',x,ROOT.RooArgList(*coeff));nb=ROOT.RooRealVar('N_continuum','',max(1,.55*N-nominal2),0,3*N);nb.setError(math.sqrt(N));pdfs.append(cont);norms.append(nb);keep.extend(coeff+[cont,nb]);nuisnames.extend(['N_continuum']+[c.GetName() for c in coeff])
   model=ROOT.RooAddPdf('sum_model','',ROOT.RooArgList(*pdfs),ROOT.RooArgList(*norms));keep.append(model)
   def poly(t):return f'(1+@0*({t})+@1*(2*({t})*({t})-1)+@2*(4*({t})*({t})*({t})-3*({t})))'
   disc='(16*@1*@1-48*@2*(@0-3*@2))'
   roots=[f'((-4*@1{sgn}sqrt({disc}))/(24*@2))' for sgn in ['+','-']]
   checks=[f'((abs({t})>=1)?1:({poly(t)}>1e-8))' for t in roots]
   cubic=f'(({disc}<=0)?1:('+ '&&'.join(checks)+'))'
   quadratic='((@1<=0)?1:((abs(@0)>=4*@1)?1:(1-@1-@0*@0/(8*@1)>1e-8)))'
   expr=f'({poly(-1)}>1e-8)&&({poly(1)}>1e-8)&&((@2==0)?{quadratic}:{cubic})'
   condition=ROOT.RooFormulaVar('positivity_indicator',expr,ROOT.RooArgList(*coeff));positive=ROOT.RooWrapperPdf('positive_continuum','',condition,True);keep.append(condition);full=ROOT.RooProdPdf('model','',ROOT.RooArgList(model,positive));keep.extend([positive,full]);getattr(w,'import')(full);getattr(w,'import')(data);mc=ROOT.RooStats.ModelConfig('ModelConfig',w);mc.SetPdf(w.pdf('model'));mc.SetObservables(ROOT.RooArgSet(w.var('mass')));mc.SetParametersOfInterest(ROOT.RooArgSet(w.var('upsilon_mu'),w.var('upsilon_sigma')));nuis=ROOT.RooArgSet(*[w.var(name) for name in nuisnames])
   mc.SetNuisanceParameters(nuis);mc.SetGlobalObservables(ROOT.RooArgSet());getattr(w,'import')(mc);w.var('mass').setVal(yp['upsilon1s']['mean']);w.writeToFile(str(d/'workspace.root'));a=dict(tag=tag,era=era,multiplier=mult,background_degree=degree,range=[LO,HI],bins=NB,events=N,nominal_mu=yp['upsilon1s']['mean'],nominal_sigma=yp['upsilon1s']['sigmaL'],fixed_shapes=yp,Y2_nominal_local_yield=nominal2,Y2_scale_bounds=[0,10],source=str(src),source_sha256=hashlib.sha256(src.read_bytes()).hexdigest(),data=counts,definition='Upsilon1S mu/sigma free, original alphaL/nL/alphaR/nR fixed; Upsilon2S original shape fixed, Y1 and continuum local yields free; Y2 scale free in original datacard bounds [0,10]; Chebyshev degree 3 on local window, coefficients c1 [-2,2], c2/c3 [-1,1], exact positivity constraint; exact bin averaging via RooBinSamplingPdf');(d/'audit.json').write_text(json.dumps(a,indent=2));configs.append(a)
(R/'configs.json').write_text(json.dumps([{k:v for k,v in c.items() if k not in ['fixed_shapes','data']} for c in configs],indent=2));print('models',len(configs))
