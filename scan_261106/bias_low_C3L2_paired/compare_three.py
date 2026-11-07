import csv,json,math
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
web=Path('/eos/user/s/soffi/www/Xee/261106/scan_SR_CR')
dest=web/'bias_low_C3L2_paired'
jobs=[r.split() for r in Path('all.txt').read_text().splitlines()]
assert len(jobs)==len({r[-1] for r in jobs})==4800
def read(path):
 with path.open() as f:return list(csv.DictReader(f))
prev=read(web/'bias_low_C3L3_paired/summary.csv')
models=['C2+L3','C3+L3','C3+L2']
rows=[]; summary=[];rng=np.random.default_rng(781)
for tag in ['env','expected']:
 for i in (4,5,6):
  cfg=f'{tag}_{i}'
  orig={x['seed']:x for x in read(web/f'bias_{tag}_0{i}_200/toys.csv')}
  for g in ('C2','L2','C3','L3'):
   sample=[]
   for group,m,gen,inj,s in jobs:
    if (group,gen)!=(cfg,g):continue
    t=orig[s];assert t['gen']==g
    b={k:dict(ok=t[k+'_ok']=='True',score=t[k+'_score'],r=t[k+'_r']) for k in ('C2','L3')}
    for k,p in [('L2',Path('.')),('C3',Path('../bias_low_C3L3_paired'))]:
     obj=json.loads((p/cfg/'results'/(s+'.json')).read_text())
     assert obj['seed']==int(s) and obj['generator']==g
     b[k]=z=obj['branches'][k]
     if z['ok']:assert math.isclose(z['score'],z['D']+{'L2':8,'C3':12}[k],abs_tol=1e-8)
    r=dict(seed=s,tag=tag,mass=float(m),generator=g,injection=float(inj))
    for name in models:
     keys=name.split('+')
     win=min(keys,key=lambda k:float(b[k]['score'])) if all(b[k]['ok'] for k in keys) else ''
     r[name]=float(b[win]['r']) if win else None
     r[name+'_winner']=win
    rows.append(r);sample.append(r)
   assert len(sample)==200
   for scope in ['own','common']:
    for model in models:
     keep=[r for r in sample if r[model]!=None and (scope=='own' or all(r[k]!=None for k in models))]
     x=np.array([r[model]-r['injection'] for r in keep])
     if scope=='own' or model==models[0]:ix=rng.integers(len(x),size=(5000,len(x)))
     y=x[ix];lo,hi=np.quantile(y.mean(1)/y.std(1,ddof=1),[.025,.975])
     bias=x.mean()/x.std(ddof=1)
     if scope=='own' and model in models[:2]:
      old=next(z for z in prev if (z['tag'],float(z['mass']),z['gen'],z['scope'],z['model'])==(tag,i/10,g,'own',['old','new'][models.index(model)]))
      assert int(old['n'])==len(x) and abs(float(old['bias'])-bias)<1e-10
     summary.append(dict(tag=tag,mass=i/10,generator=g,model=model,scope=scope,n=len(x),bias=bias,lo=lo,hi=hi,sd=x.std(ddof=1),mean_residual=x.mean(),**{k:sum(r[model+'_winner']==k for r in keep) for k in ['C2','C3','L2','L3']}))
dest.mkdir(exist_ok=True)
for name,data in [('summary',summary),('toys',rows)]:
 with (dest/(name+'.csv')).open('w') as f:
  w=csv.DictWriter(f,fieldnames=data[0]);w.writeheader();w.writerows(data)
fig,axes=plt.subplots(2,4,figsize=(14,7),sharex=True,sharey=True)
for j,tag in enumerate(['env','expected']):
 for k,g in enumerate(['C2','L2','C3','L3']):
  ax=axes[j,k];ax.axhspan(-.2,.2,color='#e4f2e4')
  for h,model in enumerate(models):
   v=[z for z in summary if (z['scope'],z['tag'],z['generator'],z['model'])==('common',tag,g,model)]
   y=np.array([z['bias'] for z in v]);ax.errorbar(np.array([.4,.5,.6])+.004*(h-1),y,yerr=[y-[z['lo'] for z in v],[z['hi'] for z in v]-y],label=model,marker='o',capsize=2)
  ax.set_title('Generator '+g);ax.set_xticks([.4,.5,.6])
  if j:ax.set_xlabel('Mass [GeV]')
  if k==0:ax.set_ylabel(('r = 0' if not j else 'Median expected r')+'\nNormalized bias (95% CI)')
axes[0,0].legend();fig.tight_layout()
for ext in ['png','pdf']:fig.savefig(dest/('comparison.'+ext))
for model in models:
 v=[z for z in summary if z['scope']=='common' and z['model']==model]
 print(model,sum(abs(z['bias'])>.2 for z in v),max(abs(z['bias']) for z in v),sum(z['lo']>-.2 and z['hi']<.2 for z in v),min(z['n'] for z in v),max(z['n'] for z in v))
