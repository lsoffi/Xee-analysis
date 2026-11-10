"""Rebuild the middle-region GoF summary from individual toy receipts."""
import argparse,csv,json,math
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
from collections import defaultdict
from scipy.stats import beta
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
p=argparse.ArgumentParser();p.add_argument('--results',type=Path,default=Path(__file__).parent/'results');p.add_argument('--observed',type=Path,default=Path(__file__).parent.parent/'gof_sr_only/results');p.add_argument('--output',type=Path,required=True);a=p.parse_args();a.output.mkdir(parents=True,exist_ok=True)
def read(f):return json.loads(f.read_text())
groups=defaultdict(dict)
with ThreadPoolExecutor(max_workers=16) as ex:
 for d in ex.map(read,a.results.glob('*.json')):
  m=d['mass'];seed=d['seed'];v=d['statistic'];assert seed not in groups[m] and math.isfinite(v) and v>=0
  assert seed in range(91000000+round(m*10)*1000,91000000+round(m*10)*1000+500)
  groups[m][seed]=v
rows=[]
for i in range(12,28):
 m=i/10;v=list(groups[m].values());assert len(v)==500,(m,len(v))
 observed=read(a.observed/f'{m:.1f}_observed.json')['statistics'][0];k=sum(x>=observed for x in v);n=len(v)
 rows.append(dict(mass=m,n=n,expected_toys=500,tail=k,p=k/n,lo=0 if k==0 else beta.ppf(.025,k,n-k+1),hi=1 if k==n else beta.ppf(.975,k+1,n-k),observed=observed,complete=True))
with (a.output/'gof.csv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
fig,ax=plt.subplots(figsize=(10,5));ax.errorbar([r['mass'] for r in rows],[r['p'] for r in rows],yerr=[[r['p']-r['lo'] for r in rows],[r['hi']-r['p'] for r in rows]],fmt='o',capsize=3,label='500 toys per mass')
ax.axhline(.05,color='firebrick',ls='--',label='p = 0.05');ax.set(xlabel='Signal mass [GeV]',ylabel='GoF p-value',ylim=(-.03,1.08),title='SR-only saturated GoF · simultaneous SR+CR fit · C4+L4');ax.legend();ax.grid(alpha=.2)
fig.text(.12,.02,'Exact binomial 95% intervals; r and trigger profiled',fontsize=8);fig.tight_layout(rect=[0,.05,1,1])
for ext in ['png','pdf']:fig.savefig(a.output/f'gof.{ext}',dpi=160)
