from pathlib import Path
import json,numpy as np,uproot
from scipy.stats import beta
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
P=Path.cwd();R=P/'runs/20261030_sr6_cr4_scan';d=R/'gof'
def read(directory):
 p=next(directory.glob('higgsCombine*.root'));return uproot.open(p)['limit']['limit'].array(library='np')
obs=float(read(d/'observed')[0]);toys=read(d/'toys');toys=toys[np.isfinite(toys)];N=len(toys);ex=int(sum(toys>=obs));p=ex/N if N else None;interval=[float(beta.ppf(.025,ex,N-ex+1)) if ex else 0,float(beta.ppf(.975,ex+1,N-ex)) if ex<N else 1];out=dict(observed_raw=obs,observed_displayed=obs,N_requested=100,N_finite=N,exceedances=ex,p=p,CP95=interval,definition='Combine saturated statistic, envelope; identical constant penalty 28 in data/toys; B-only common to all hypotheses',preliminary=True)
(d/'summary.json').write_text(json.dumps(out,indent=2));print(out)
fig,axs=plt.subplots(1,2,figsize=(13,4.5));axs[0].hist(toys,bins=18,color='steelblue',label=f'{N}/100 toy');axs[0].axvline(obs,color='firebrick',label=f'Dati: {obs:.2f}');axs[0].set_xlabel('Statistica saturated di Combine');axs[0].set_ylabel('Toy');axs[0].set_title(f'GoF B-only: p={p:.3f}, intervallo MC 95% [{interval[0]:.3f}, {interval[1]:.3f}]');axs[0].legend()
rows=[]
for dd in R.glob('m*/joint/gof_sb'):
 if not (dd/'receipt.json').exists():continue
 try:rows.append(dict(mass=float(dd.parent.parent.name[1:]),D=float(read(dd)[0])))
 except Exception:pass
rows.sort(key=lambda x:x['mass']);(d/'sb_observed.json').write_text(json.dumps(rows,indent=2));axs[1].axhline(obs,color='gray',ls='--',label='B-only comune');sbmap={a['mass']:a['D'] for a in rows};grid=np.arange(60,91)/10;axs[1].plot(grid,[sbmap.get(m,np.nan) for m in grid],'o-',label='S+B osservato');axs[1].set_xlim(5.95,9.05);axs[1].set_xlabel('Massa ipotizzata [GeV]');axs[1].set_ylabel('Statistica saturated di Combine');axs[1].set_title('Diagnostica S+B; nessun pGoF S+B calibrato per massa');axs[1].legend();axs[1].grid(alpha=.2);fig.suptitle('SR envelope Cheb6/Bern6 + CR Cheb4, fit simultaneo');fig.tight_layout()
for ext in ['png','pdf']:fig.savefig(R/f'gof.{ext}',dpi=150)
plt.close(fig)
fig,ax=plt.subplots(figsize=(9,4));x=np.arange(60,91)/10;ax.axhspan(interval[0],interval[1],alpha=.15,color='steelblue',label='Precisione Monte Carlo (95%)');ax.plot(x,np.full(len(x),p),color='steelblue',label='pGoF B-only: stesso test a ogni massa');ax.set_ylim(0,1);ax.set_xlabel('Ipotesi di massa [GeV]');ax.set_ylabel('p-value goodness of fit');ax.legend();ax.grid(alpha=.2);fig.tight_layout()
for ext in ['png','pdf']:fig.savefig(R/f'gof_vs_mass.{ext}',dpi=150)
plt.close(fig)
