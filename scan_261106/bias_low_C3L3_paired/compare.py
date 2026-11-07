import csv,json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
web=Path('/eos/user/s/soffi/www/Xee/261106/scan_SR_CR')
dest=web/'bias_low_C3L3_paired';dest.mkdir(exist_ok=True)
jobs=[l.split() for l in Path('all.txt').read_text().splitlines()]
assert len(jobs)==len({r[-1] for r in jobs})==4800
rng=np.random.default_rng(781)
rows=[];summary=[];pairs=[]
def stats(x,ix):
 b=x[ix];v=b.mean(1)/b.std(1,ddof=1)
 return x.mean()/x.std(ddof=1),v
for tag in ['env','expected']:
 for i in [4,5,6]:
  group=f'{tag}_{i}';m=i/10;idx={r['seed']:r for r in csv.DictReader((web/f'bias_{tag}_0{i}_200/toys.csv').open())}
  for g in ['C2','L2','C3','L3']:
   sample=[]
   for cfg,mass,gen,inj,seed in jobs:
    if cfg!=group or gen!=g:continue
    t=idx[seed];o=json.loads((Path(cfg)/'results'/(seed+'.json')).read_text())
    assert t['gen']==o['generator']==g and o['seed']==int(seed) and float(mass)==m
    c=o['branches']['C3'];l=t['L3_ok']=='True'
    win=('C3' if c['score']<=float(t['L3_score']) else 'L3') if c['ok'] and l else ''
    val=c['r'] if win=='C3' else float(t['L3_r']) if win else None
    row=dict(tag=tag,mass=m,gen=g,seed=seed,injection=float(inj),oldwin=t['winner'],newwin=win,old=float(t['r']) if t['winner'] else None,new=val)
    rows.append(row);sample.append(row)
   assert len(sample)==200
   for scope in ['own','paired']:
    reps=[]
    for key in ['old','new']:
     keep=[r for r in sample if r[key] is not None and (scope=='own' or None not in (r['old'],r['new']))]
     x=np.array([r[key]-r['injection'] for r in keep])
     if key=='old' or scope=='own':ix=rng.integers(len(x),size=(5000,len(x)))
     bias,b=stats(x,ix);reps.append((bias,b))
     lo,hi=np.quantile(b,[.025,.975])
     summary.append(dict(tag=tag,mass=m,gen=g,scope=scope,model=key,n=len(x),bias=bias,lo=lo,hi=hi,mean_residual=x.mean(),sd=x.std(ddof=1),L3wins=sum(r[key+'win']=='L3' for r in keep)))
    if scope=='paired':
     lo,hi=np.quantile(reps[1][1]-reps[0][1],[.025,.975]);pairs.append(dict(tag=tag,mass=m,gen=g,n=len(x),delta=reps[1][0]-reps[0][0],lo=lo,hi=hi))
for name,data in [('toys',rows),('summary',summary),('paired_delta',pairs)]:
 with (dest/(name+'.csv')).open('w') as f:
  w=csv.DictWriter(f,fieldnames=data[0]);w.writeheader();w.writerows(data)
fig,axes=plt.subplots(2,4,figsize=(14,6),sharex=True,sharey=True)
for j,tag in enumerate(['env','expected']):
 for k,g in enumerate(['C2','L2','C3','L3']):
  ax=axes[j,k];ax.axhspan(-.2,.2,color='#e4f2e4');ax.axhline(0,color='gray',lw=.6)
  for model,offset,label in [('old',-.002,'C2 + L3'),('new',.002,'C3 + L3')]:
   v=[r for r in summary if (r['scope'],r['tag'],r['gen'],r['model'])==('own',tag,g,model)]
   y=np.array([r['bias'] for r in v]);ax.errorbar(np.array([r['mass'] for r in v])+offset,y,yerr=[y-[r['lo'] for r in v],[r['hi'] for r in v]-y],marker='o',capsize=2,label=label)
  ax.set_xticks([.4,.5,.6]);ax.set_title('Generator '+g)
  if j:ax.set_xlabel('Mass [GeV]')
  if k==0:ax.set_ylabel(('r = 0' if not j else 'Median expected r')+'\nNormalized bias')
axes[0,0].legend(fontsize=8);fig.tight_layout()
for ext in ['png','pdf']:fig.savefig(dest/('comparison.'+ext))
for tag in ['env','expected']:
 for model in ['old','new']:
  v=[r for r in summary if (r['scope'],r['tag'],r['model'])==('own',tag,model)]
  print(tag,model,sum(r['n'] for r in v),sum(abs(r['bias'])>.2 for r in v),[(r['gen'],r['mass'],round(r['bias'],3)) for r in v if abs(r['bias'])>.2])
 v=[r for r in summary if r['scope']=='paired' and r['tag']==tag]
 a=np.array([r['sd'] for r in v if r['model']=='new'])/np.array([r['sd'] for r in v if r['model']=='old']);print('paired sd ratio',a.min(),a.max())
