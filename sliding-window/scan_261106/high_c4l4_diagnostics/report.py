import ROOT as R, pathlib,json,sys,math
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
m=sys.argv[1];run=pathlib.Path(__file__).resolve().parent;p=run/('m'+m)
meta=json.load(open(p/'C4.json'));f=R.TFile.Open(str(next(p.glob('fitDiagnostics*.root'))))
out=pathlib.Path('/eos/user/s/soffi/www/Xee/261106/scan_SR_CR/high_c4l4_diagnostics');out.mkdir(parents=True,exist_ok=True)
summary={'mass':float(m),'stages':{}}
for stage in ['b','s']:
 fit=f.Get('fit_'+stage);assert fit
 nchan=len(meta['channels']);fig,axs=plt.subplots(nchan,2,figsize=(13,2.5*nchan),gridspec_kw={'width_ratios':[3,1]})
 rows=[]
 for k,ch in enumerate(meta['channels']):
  label=ch['label'];d=f.Get('shapes_fit_'+stage+'/'+label);assert d
  n=len(ch['counts']);g=d.Get('data');assert g and g.GetN()==n or g and g.GetN()>=n
  y=np.array([g.GetPointY(i) for i in range(n)]);h=d.Get('total');assert h
  mu=np.array([h.GetBinContent(i+1) for i in range(n)])
  valid=np.all(np.isfinite(mu)) and np.all(mu>0)
  chi=float(np.sum((y-mu)**2/mu)) if valid else None
  pull=(y-mu)/np.sqrt(mu) if valid else np.full(n,np.nan)
  a,b=axs[k];edges=np.asarray(ch["edges"]);x=edges.mean(axis=1) if edges.ndim==2 else (edges[:-1]+edges[1:])/2
  a.errorbar(x,y,yerr=np.sqrt(y),fmt='.',color='black',ms=2,label='Data')
  for key in d.GetListOfKeys():
   name=key.GetName();hh=d.Get(name)
   if not hh.InheritsFrom('TH1') or name in ['total_background','total_signal','total_covar']:continue
   v=np.array([hh.GetBinContent(i+1) for i in range(n)])
   a.step(x,v,where='mid',lw=1.5 if name=='total' else 1,label=name)
  a.set_title(label+(' · Pearson χ²/Nbins = %.1f/%d = %.2f'%(chi,n,chi/n) if valid else ' · INVALID expected counts'),fontsize=10)
  a.set_ylabel('Events / bin');a.legend(fontsize=6,ncol=5)
  b.axhline(0,color='black',lw=.8);b.axhspan(-2,2,color='gray',alpha=.15);b.plot(x,pull,'.',ms=2);b.set_ylabel('(data − model)/√model')
  for ax in [a,b]:ax.set_xlabel('Dielectron mass [GeV]')
  rows.append({'channel':label,'bins':n,'pearson_chi2':chi,'chi2_per_bin':chi/n if valid else None})
 summary['stages'][stage]={'status':fit.status(),'covQual':fit.covQual(),'edm':fit.edm(),'channels':rows,'nfloat':fit.floatParsFinal().getSize()}
 fig.suptitle(f'm = {m} GeV · envelope C4+L4 · fit_{stage} · r free in S+B\nSimultaneous SR+selected CRs · total C0 at shared boundaries · status {fit.status()}, covQual {fit.covQual()}\nPer-channel χ²/Nbins is diagnostic, not χ²/ndf (parameters are shared).',fontsize=12)
 fig.tight_layout(rect=[0,0,1,.96]);fig.savefig(out/f'm{m}_{stage}.png',dpi=110);fig.savefig(out/f'm{m}_{stage}.pdf');plt.close(fig)
(p/'diagnostics.json').write_text(json.dumps(summary,indent=2));(out/f'm{m}.json').write_text(json.dumps(summary,indent=2))
