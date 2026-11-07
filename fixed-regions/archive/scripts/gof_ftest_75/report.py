from pathlib import Path
import json,re
import numpy as np,uproot
from scipy.stats import chi2,beta
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
P=Path.cwd();R=P/'runs/20261029_gof_ftest_75';N=1000

def read(d):
 with uproot.open(next(d.glob('higgsCombine*.root'))) as f:a=f['limit'].arrays(library='np')
 assert np.all(np.isfinite(a['limit'])) and np.all(a['limit']>=0)
 if d==R/'degree2/toys_1000/fit_degree2':
  candidates=[a['limit'].copy()]
  for start in [-.05,.05]:
   with uproot.open(next((d.parent/f'uniform_start{start}').glob('higgsCombine*.root'))) as f:b=f['limit'].arrays(library='np')
   assert np.array_equal(a['iToy'],b['iToy']) and len(b['limit'])==N
   candidates.append(b['limit'])
  best=np.min(candidates,axis=0);np.savez(R/'degree2_multistart_audit.npz',toy_id=a['iToy'],original=candidates[0],minus=candidates[1],plus=candidates[2],best=best)
  (R/'multistart_summary.json').write_text(json.dumps(dict(rule='Minimum deviance across original fit and two uniform coefficient shifts +/-0.05; all 1000 null-degree2 toys refitted with strategy2/tol1e-5.',improved_over_0p01=int(sum(a['limit']-best>.01)),maximum_improvement=float(max(a['limit']-best))),indent=2));a['limit']=best
 return a

def tail(a,obs):
 k=int(sum(a>=obs));n=len(a);lo=0 if k==0 else float(beta.ppf(.025,k,n-k+1));hi=1 if k==n else float(beta.ppf(.975,k+1,n-k));return dict(exceedances=k,toys=n,p_raw=k/n,p_plus_one=(k+1)/(n+1),interval95=[lo,hi],upper95_one_sided=1-.05**(1/n) if k==0 else None)
obs={k:float(read(R/k/'observed')['limit'][0]) for k in ['degree2','degree3','degree2SR_3CR']};out={'observed_deviance':obs,'gof':{},'comparisons':{},'regional_deviance':{},'regional_note':'Diagnostic Poisson deviance reconstructed from saved binned shapes, not the full Combine likelihood statistic; differences from total GoF are recorded, not treated as exact equality.','checks':{}}
samples={}
for k,fit in [('degree2','2'),('degree3','3'),('degree2SR_3CR','mixed')]:
 a=read(R/k/f'toys_{N}/fit_degree{fit}')['limit'];assert len(a)==N;samples[k]=a;out['gof'][k]=tail(a,obs[k])
for name,generator,low,high,dof in [('joint_2_to_3','degree2','2','3',8),('SR_2_to_3_given_CR3','degree2SR_3CR','mixed','3',4)]:
 a=read(R/generator/f'toys_{N}/fit_degree{low}');b=read(R/generator/f'toys_{N}/fit_degree{high}');assert np.array_equal(a['iToy'],b['iToy']) and np.array_equal(a['iSeed'],b['iSeed']);assert len(a['limit'])==N
 delta=a['limit']-b['limit'];assert min(delta)>-.001;observed=obs[generator]-obs['degree3'];fstat=(delta/dof)/(b['limit']/1024);fobs=(observed/dof)/(obs['degree3']/1024)
 out['comparisons'][name]=dict(delta_deviance=observed,added_parameters=dof,toy_LR=tail(delta,observed),asymptotic_chi2_p=float(chi2.sf(observed,dof)),F=fobs,toy_F=tail(fstat,fobs),minimum_toy_delta=float(min(delta)))
 samples[name]=delta;np.savez(R/f'{name}_toys.npz',delta=delta,F=fstat,toy_id=a['iToy'])
# Exact extended-likelihood regional contributions at the B-only fitted parameters.
out['regional_deviance']=json.loads((R/'regional_exact.json').read_text())
out['regional_note']='Extended likelihood contributions evaluated from PDF densities at observed bin centers, including analytic expected yields; sum verified against Combine saturated statistic. Regional contributions are not independent p-values.'
for key,vals in out['regional_deviance'].items():
 out['checks'][key+'_deviance_agreement']=sum(v['total'] for v in vals.values())-obs[key]
 assert abs(out['checks'][key+'_deviance_agreement'])<.001
# Audit logs and receipts; do not silently discard bad toys.
for p in R.glob('*/toys_1000/fit_degree*/receipt.json'):
 assert json.loads(p.read_text())['exit_code']==0
 log=(p.parent/'combine.log').read_text();bad=re.findall(r'^.*(?:Minimization failed|ERROR|Failed to|failed to converge).*$',log,re.M);assert not bad,(str(p),bad[:3])
out['definition']='B-only r=0, original shared Upsilon scales and nuisance constraints; frequentist conditional toys at each null B-only fit; same toys and auxiliary snapshots paired across degrees. Saturated GoF. Nominal F denominator df=1064-40=1024 (10 constrained nuisances and 10 auxiliary constraints cancel). Toy calibration used for conclusions.'
out['multistart']=json.loads((R/'multistart_summary.json').read_text())
(R/'summary.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2))
fig,axs=plt.subplots(1,2,figsize=(12,4.5))
for ax,key,title in zip(axs,['degree2','degree3'],['Chebyshev grado 2','Chebyshev grado 3']):
 ax.hist(samples[key],bins=35,color='#0072B2',alpha=.7);ax.axvline(obs[key],color='#D55E00',lw=2,label=f'Osservato: {obs[key]:.2f}');v=out['gof'][key];ax.set_title(f'{title}\np = {v["p_raw"]:.3f} ({v["exceedances"]}/1000)');ax.set_xlabel('Devianza rispetto al modello saturo');ax.set_ylabel('Pseudoesperimenti');ax.legend()
fig.suptitle('Goodness of fit B-only — SR + CR, quattro ere');fig.tight_layout();fig.savefig(R/'gof.png',dpi=160);fig.savefig(R/'gof.pdf');plt.close(fig)
fig,axs=plt.subplots(1,2,figsize=(12,4.5))
for ax,key,title in zip(axs,['joint_2_to_3','SR_2_to_3_given_CR3'],['Grado 2 → 3 in SR e CR','Grado 2 → 3 in SR; CR fissata al grado 3']):
 v=out['comparisons'][key];ax.hist(samples[key],bins=35,color='#009E73',alpha=.7);ax.axvline(v['delta_deviance'],color='#D55E00',lw=2,label=f'Osservato: {v["delta_deviance"]:.2f}');ax.set_title(title+f'\nSuperamenti: {v["toy_LR"]["exceedances"]}/1000');ax.set_xlabel('Δ devianza = 2 log(L grado alto / L grado basso)');ax.set_ylabel('Pseudoesperimenti sotto il modello di grado basso');ax.legend()
fig.suptitle('Confronto tra gradi — likelihood ratio calibrato con toy');fig.tight_layout();fig.savefig(R/'degree_test.png',dpi=160);fig.savefig(R/'degree_test.pdf');plt.close(fig)
body='<meta charset="utf-8"><style>body{font:17px system-ui;max-width:1200px;margin:30px auto}img{width:100%}pre{white-space:pre-wrap}</style><h1>GoF e confronto Chebyshev 2/3 a massa7,5</h1><p>SR6,373–7,920 e CR8,8–10,8GeV; quattro ere, binning originale, r=0, Υoriginali. 1000toy frequentist per ogni generatore, paired per confronto. Test globale e test mirato SR con CRgrado3. Soglia di lettura5%; calibrazione condizionata ai dati. Non è un test di bias o copertura del segnale e non stabilisce il grado ottimale fra tutti quelli possibili.</p>'
for im in ['gof','degree_test']:body+=f'<img src="{im}.png">'
body+='<h2>Risultati completi e contributi regionali</h2><pre>'+json.dumps(out,indent=2)+'</pre>';(R/'index.html').write_text(body)
