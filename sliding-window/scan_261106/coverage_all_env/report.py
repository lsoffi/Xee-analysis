import csv,json,math
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import beta
base=Path(__file__).resolve().parent
out=Path('/eos/user/s/soffi/www/Xee/261106/scan_SR_CR/coverage_all_env');out.mkdir(parents=True,exist_ok=True)
med=json.loads((base/'medians.json').read_text());summary=[];toys=[];seen=set()
for m,inj in med.items():
 rows=[]
 for p in (base/('m'+m)/'results').glob('*.json'):
  r=json.loads(p.read_text());assert r['seed'] not in seen;seen.add(r['seed'])
  assert r['injection']==inj and math.isfinite(r['CLs']) and 0<=r['CLs']<=1
  assert r['excluded']==(r['CLs']<.05)
  r={'mass':float(m),**r};rows.append(r);toys.append(r)
 for g in ('C4','L4'):
  rr=[r for r in rows if r['generator']==g];n=len(rr);assert n<=200
  if not n:continue
  k=sum(r['excluded'] for r in rr);covered=n-k
  lo=beta.ppf(.025,covered,n-covered+1) if covered else 0
  hi=beta.ppf(.975,covered+1,n-covered) if covered<n else 1
  summary.append(dict(mass=float(m),generator=g,injection=inj,valid=n,missing=200-n,excluded=k,coverage=covered/n,lo=lo,hi=hi))
for name,data in [('coverage',summary),('toys',toys)]:
 if data:
  with (out/(name+'.csv')).open('w') as f:
   w=csv.DictWriter(f,fieldnames=data[0]);w.writeheader();w.writerows(data)
complete=len(toys)==6400
fig,(ax,counts)=plt.subplots(2,1,figsize=(11,7),sharex=True,gridspec_kw={'height_ratios':[3,1]})
ax.axhline(.95,ls='--',color='black',label='Nominal 95%')
for g,off,color in [('C4',-.007,'#1766ad'),('L4',.007,'#d26b19')]:
 rr=[r for r in summary if r['generator']==g and r['valid']==200]
 x=[r['mass']+off for r in rr];y=np.array([r['coverage'] for r in rr])
 ax.errorbar(x,y,yerr=[y-[r['lo'] for r in rr],[r['hi'] for r in rr]-y],fmt='o-',capsize=3,color=color,label=g+' generator')
 counts.plot(x,[r['valid'] for r in rr],'o',color=color)
ax.set(ylabel='Fraction not excluding injected r',ylim=(.85,1.005),title='C4 + L4 envelope; trigger profiled\nMedian expected signal injected; exact binomial 95% intervals'+('' if complete else ' — preliminary'))
ax.legend();ax.grid(alpha=.2);counts.set(xlabel='Signal mass [GeV]',ylabel='Valid toys',xticks=np.arange(1.2,2.8,.1),ylim=(0,215));counts.grid(alpha=.2)
fig.tight_layout()
for ext in ('png','pdf'):fig.savefig(out/('coverage_vs_mass.'+ext),dpi=170)
(out/'status.json').write_text(json.dumps({'complete':complete,'results':len(toys),'expected':6400},indent=2))
(out/'index.html').write_text('<!doctype html><html lang="en"><meta charset="utf-8"><title>Coverage across masses</title><h1>Coverage of the C4+L4 CLs procedure</h1><p>Trigger profiled. Median expected signal injected. A toy excludes its true injection when CLs &lt; 0.05. Bars are pointwise exact Clopper-Pearson 95% intervals. Only groups with all 200 results appear on the plot.</p><p>'+str(len(toys))+'/6400 results.</p><img style="max-width:100%" src="coverage_vs_mass.png"><p><a href="coverage_vs_mass.pdf">PDF</a> | <a href="coverage.csv">CSV</a></p></html>')
print('RESULTS',len(toys),'complete',complete)
for r in summary:
 if r['valid']==200:print(r)
