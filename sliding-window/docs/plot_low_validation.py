import json,glob,csv
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import beta
base=Path(__file__).resolve().parent.parent
out=Path('/eos/user/s/soffi/www/Xee/261106/scan_SR_CR/low_validation');out.mkdir(exist_ok=True)
d=json.load(open(base/'task_monitor/low_aggregate.json'))
plt.rcParams.update({'font.size':11})
fig,ax=plt.subplots(figsize=(9,5.5))
for g,offset,color in [('C3',-.006,'#2166ac'),('L3',.006,'#d6604d')]:
 rr=[r for r in d['coverage'] if r['generator']==g];x=np.array([r['mass'] for r in rr])+offset;y=np.array([r['coverage'] for r in rr]);lo=np.array([r['CI95'][0] for r in rr]);hi=np.array([r['CI95'][1] for r in rr])
 ax.errorbar(x,y*100,yerr=np.array([y-lo,hi-y])*100,fmt='o',capsize=5,color=color,label=f'{g} generator; C3+L3 fit')
ax.axhline(95,color='black',ls='--',label='95% reference');ax.set(xlabel='Signal mass hypothesis [GeV]',ylabel='Toys not excluding the injected signal [%]',xticks=[.4,.5,.6],ylim=(88,100),xlim=(.37,.63),title='Low SR: CLs coverage at median-expected signal')
ax.grid(alpha=.2);ax.legend(loc='lower left',fontsize=10)
fig.text(.5,.025,'200 toys per mass and generator; exact 95% binomial intervals. Trigger profiled.',ha='center',fontsize=10)
fig.tight_layout(rect=(0,.06,1,1))
for ext in ['png','pdf']:fig.savefig(out/f'coverage_low.{ext}',dpi=160)
plt.close(fig)
fig,axes=plt.subplots(1,3,figsize=(15,4.8))
for ax,r in zip(axes,d['gof']):
 m=r['mass'];vals=np.array([json.load(open(f))['statistic'] for f in glob.glob(str(base/f'gof_low_sr/results/{m:.1f}_*.json'))]);obs=r['observed'];k=r['tail'];n=r['n'];lo=beta.ppf(.025,k,n-k+1) if k else 0;hi=beta.ppf(.975,k+1,n-k) if k<n else 1;r['p_CI95']=[float(lo),float(hi)]
 edges=np.linspace(min(vals.min(),obs)*.97,max(vals.max(),obs)*1.03,28)
 ax.hist(vals,bins=edges,color='#9ecae1',edgecolor='white',label='Frequentist toys')
 ax.hist(vals[vals>=obs],bins=edges,color='#ef8a62',edgecolor='white',label='Toy statistic >= data')
 ax.axvline(obs,color='#b2182b',lw=2,label='Observed')
 ax.set(title=f'm = {m:.1f} GeV',xlabel='SR-only saturated GoF statistic',ylabel='Toys')
 ax.text(.97,.96,f'p = {k}/{n} = {k/n:.3f}\n95% CI [{lo:.3f}, {hi:.3f}]',transform=ax.transAxes,ha='right',va='top',bbox=dict(facecolor='white',alpha=.9,edgecolor='none'));ax.grid(axis='y',alpha=.15)
axes[0].legend(fontsize=8,loc='upper left')
fig.suptitle('Low SR: C3+L3 envelope goodness of fit',fontsize=16)
fig.text(.5,.025,'Simultaneous SR+CR fit, free r, profiled trigger; statistic evaluated in SR only. 500 toys per mass.',ha='center',fontsize=11)
fig.tight_layout(rect=(0,.07,1,.93))
for ext in ['png','pdf']:fig.savefig(out/f'gof_low.{ext}',dpi=160)
plt.close(fig)
(out/'results.json').write_text(json.dumps(d,indent=2))
for key in ['coverage','gof']:
 with (out/f'{key}.csv').open('w') as f:
  w=csv.DictWriter(f,fieldnames=list(d[key][0]));w.writeheader();w.writerows(d[key])
(out/'index.html').write_text('<html><body><h1>Low SR validation</h1><p>Complete samples. Exact 95% binomial intervals. Coverage is conditional on the tested generators and median-expected injections. The local GoF at 0.6 GeV shows tension (p=0.006); this is not a signal significance.</p><h2>Coverage</h2><img width="900" src="coverage_low.png"><h2>Goodness of fit</h2><img width="1200" src="gof_low.png"><p><a href="coverage_low.pdf">Coverage PDF</a> | <a href="gof_low.pdf">GoF PDF</a> | <a href="results.json">Results</a></p></body></html>')
print(out)
