import json,csv,math,datetime
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
from collections import defaultdict,Counter
import numpy as np
from scipy.stats import beta
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
base=Path('/afs/cern.ch/work/s/soffi/DIELE-WORK/Analisi-statistica-Xee/scan_261106')
import os
out=Path(os.environ['PLOT_OUTPUT']);out.mkdir(parents=True,exist_ok=True)
stamp=datetime.datetime.now(datetime.timezone.utc).strftime('%Y-%m-%d %H:%M UTC')
def read(p):
 try:return json.loads(p.read_text())
 except Exception:return {'read_error':str(p)}
def load(paths):
 with ThreadPoolExecutor(max_workers=16) as ex:return list(ex.map(read,paths))
def savecsv(name,rows):
 if rows:
  with (out/(name+'.csv')).open('w') as f:
   w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
def savefig(name,fig):
 fig.savefig(out/(name+'.png'),dpi=160);fig.savefig(out/(name+'.pdf'));plt.close(fig)
issues=[]

lim=[]
for f in (base/'high_c4l4_limits_bias/limits').glob('*.json'):
 d=read(f);lim.append({'mass':d['mass'],'median':next(x['limit'] for x in d['quantiles'] if abs(x['quantile']-.5)<1e-5)})
print('Reading bias',flush=True)
manifest=read(base/'high_c4l4_limits_bias/manifest.json');allowed={d['seed']:(d['mass'],d['generator'],d['tag']) for d in manifest}
paths=list((base/'high_c4l4_limits_bias').glob('m*/results/*.json'));data=load(paths);groups=defaultdict(list);seen=set()
for d in data:
 if 'read_error' in d:issues.append(d);continue
 seed=d['seed'];key=(d['mass'],d['generator'],d['tag']);assert allowed[seed]==key and seed not in seen;seen.add(seed)
 b=d['branches'];ok=all(v.get('ok',False) and v['status']==0 and math.isfinite(v['score']) and abs(v['r'])<79.999 and abs(v['score']-v['D']-16)<1e-6 for v in b.values())
 assert set(b)=={'C4','L4'}
 win=min(b,key=lambda k:b[k]['score']) if ok else None
 assert bool(d['valid'])==ok and d['winner']==win
 inj=d['injection'];assert inj==0 if d['tag']=='r0' else inj>0
 if d['tag']=='expected':
  expect=next(r['median'] for r in lim if abs(r['mass']-d['mass'])<1e-6);assert abs(inj-expect)<1e-8
 groups[key].append((b[win]['r']-inj if ok else None,win,inj))
summary=[];rng=np.random.default_rng(90110)
for key,rr in sorted(groups.items()):
 m,g,tag=key;x=np.array([r[0] for r in rr if r[0] is not None]);n=len(x);mean=sd=bias=lo=hi=float('nan')
 if n>=2:
  mean=float(x.mean());sd=float(x.std(ddof=1))
  if sd>0:bias=mean/sd
 if n>=20 and sd>0:
  rep=rng.choice(x,(2000,n));s=rep.std(axis=1,ddof=1);v=rep.mean(axis=1)[s>0]/s[s>0];lo,hi=np.quantile(v,[.025,.975])
 summary.append(dict(mass=m,generator=g,tag=tag,received=len(rr),valid=n,invalid=len(rr)-n,pending=200-len(rr),injection=rr[0][2],mean_residual=mean,sd=sd,bias=bias,lo=lo,hi=hi,C4wins=sum(r[1]=='C4' for r in rr),L4wins=sum(r[1]=='L4' for r in rr)))
savecsv('bias',summary)
for tag in ['r0','expected']:
 fig,(ax,ct)=plt.subplots(2,1,figsize=(13,7),sharex=True,gridspec_kw={'height_ratios':[3,1]})
 ax.axhspan(-.2,.2,color='#deefde');ax.axhline(0,color='gray',lw=.7)
 for g,color,offset in [('C4','#2878b5',-.014),('L4','#d07825',.014)]:
  rr=[r for r in summary if r['tag']==tag and r['generator']==g]
  for complete,marker in [(True,'o'),(False,'s')]:
   vv=[r for r in rr if (r['received']==200)==complete and r['valid']>=20 and math.isfinite(r['lo'])]
   if vv:ax.errorbar([r['mass']+offset for r in vv],[r['bias'] for r in vv],yerr=[[max(0,r['bias']-r['lo']) for r in vv],[max(0,r['hi']-r['bias']) for r in vv]],fmt=marker,color=color,ms=4,capsize=2,mfc=color if complete else 'white',label=g+(' · 200 received' if complete else ' · partial'))
  ct.scatter([r['mass']+offset for r in rr],[r['valid'] for r in rr],s=12,color=color,label=g)
 ax.set(ylabel='Mean(r̂ − injected r) / sample SD',title='C4+L4 envelope bias · '+('r = 0' if tag=='r0' else 'r = median expected limit')+' · provisional')
 ax.legend(fontsize=8,ncol=2);ax.grid(alpha=.2);ct.axhline(200,color='gray',ls='--',lw=.7);ct.set(ylim=(-5,215),ylabel='Valid toys',xlabel='Signal mass [GeV]',xlim=(4.,9.1));ct.grid(alpha=.2)
 fig.text(.1,.015,stamp+' · both branches must be valid · bootstrap 95% intervals · points with <20 valid toys omitted · trigger profiled',fontsize=8);fig.tight_layout(rect=[0,.04,1,1]);savefig('bias_'+tag,fig)

(out/'status.json').write_text(json.dumps(dict(updated=stamp,received=len(seen),expected=40000,valid=sum(r['valid'] for r in summary),invalid=sum(r['invalid'] for r in summary),issues=issues),indent=2))
print('BIAS',len(seen),flush=True)
