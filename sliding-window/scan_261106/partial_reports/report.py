import json,csv,math,datetime
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
from collections import defaultdict,Counter
import numpy as np
from scipy.stats import beta
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
base=Path(__file__).resolve().parent.parent
out=Path('/eos/user/s/soffi/www/Xee/261106/scan_SR_CR/partial_reports');out.mkdir(parents=True,exist_ok=True)
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
print('Reading GoF',flush=True)
data=load(list((base/'gof_sr_single_toy/results').glob('*.json')));groups=defaultdict(dict)
for d in data:
 if 'read_error' in d:issues.append(d);continue
 m=d['mass'];seed=d['seed'];v=d['statistic']
 assert seed not in groups[m] and math.isfinite(v) and v>=0
 groups[m][seed]=v
rows=[]
for i in range(12,28):
 m=i/10;v=list(groups[m].values());n=len(v)
 if not n:continue
 obs=read(base/f'gof_sr_only/results/{m:.1f}_observed.json')['statistics'][0];k=sum(t>=obs for t in v)
 lo=0 if k==0 else beta.ppf(.025,k,n-k+1);hi=1 if k==n else beta.ppf(.975,k+1,n-k)
 rows.append(dict(mass=m,n=n,expected_toys=500,tail=k,p=k/n,lo=lo,hi=hi,observed=obs,complete=n==500))
savecsv('gof',rows)
fig,ax=plt.subplots(figsize=(10,5))
for complete,marker in [(True,'o'),(False,'s')]:
 rr=[r for r in rows if r['complete']==complete]
 if rr:ax.errorbar([r['mass'] for r in rr],[r['p'] for r in rr],yerr=[[r['p']-r['lo'] for r in rr],[r['hi']-r['p'] for r in rr]],fmt=marker,ms=5,capsize=3,label='500 toys' if complete else 'Partial toy sample')
ax.axhline(.05,color='firebrick',ls='--',lw=1,label='p = 0.05');ax.set(xlabel='Signal mass [GeV]',ylabel='GoF p-value',ylim=(-.03,1.08),title='SR-only saturated GoF · simultaneous SR+CR fit · C4+L4')
ax.legend();ax.grid(alpha=.2);fig.text(.12,.02,stamp+' · exact binomial 95% intervals; r and trigger profiled',fontsize=8);fig.tight_layout(rect=[0,.05,1,1]);savefig('gof',fig)
print('GoF done',sum(r['n'] for r in rows),flush=True)
print('Reading limits',flush=True)
lim=[]
for folder in [base/'limits_recommended_envelope/results',base/'high_c4l4_limits_bias/limits']:
 for d in load(list(folder.glob('*.json'))):
  if 'read_error' in d:issues.append(d);continue
  q=d['quantiles'];vals=[]
  for level in [.025,.16,.5,.84,.975,-1]:
   hits=[r['limit'] for r in q if abs(r['quantile']-level)<1e-5];assert len(hits)==1 and math.isfinite(hits[0]) and hits[0]>=0;vals.append(hits[0])
  assert vals[:5]==sorted(vals[:5]);lim.append(dict(zip(['mass','minus2','minus1','median','plus1','plus2','observed'],[d['mass']]+vals)))
lim.sort(key=lambda r:r['mass']);savecsv('limits',lim);assert len({r['mass'] for r in lim})==len(lim)
fig,ax=plt.subplots(figsize=(12,6));chunks=[]
for r in lim:
 if not chunks or r['mass']-chunks[-1][-1]['mass']>.11:chunks.append([])
 chunks[-1].append(r)
for i,rr in enumerate(chunks):
 x=[r['mass'] for r in rr];ax.fill_between(x,[r['minus2'] for r in rr],[r['plus2'] for r in rr],color='#f6cc43',label='Expected ±2σ' if i==0 else None);ax.fill_between(x,[r['minus1'] for r in rr],[r['plus1'] for r in rr],color='#69bb63',label='Expected ±1σ' if i==0 else None)
 ax.plot(x,[r['median'] for r in rr],'k--',label='Median expected' if i==0 else None);ax.plot(x,[r['observed'] for r in rr],'ko-',ms=3,label='Observed' if i==0 else None)
for a,b,label in [(.65,1.155,'ω/φ CR'),(2.8,4.,'Charmonium CR')]:ax.axvspan(a,b,color='gray',alpha=.13);ax.text((a+b)/2,.97,label,transform=ax.get_xaxis_transform(),ha='center',va='top',fontsize=9)
ax.set(yscale='log',xlabel='Signal mass [GeV]',ylabel='95% CL upper limit on r',title=f'Sliding-window envelope limits · {len(lim)}/69 masses available')
ax.legend(fontsize=9);ax.grid(which='both',axis='y',alpha=.18);fig.text(.12,.02,stamp+' · C3+L3 below φ; C4+L4 above φ · missing masses left as gaps · high-mass model validation ongoing',fontsize=8);fig.tight_layout(rect=[0,.05,1,1]);savefig('limits',fig)
print('Limits done',len(lim),flush=True)
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
status=dict(updated=stamp,gof_toys=sum(r['n'] for r in rows),limits=len(lim),bias_received=len(seen),bias_valid=sum(r['valid'] for r in summary),bias_invalid=sum(r['invalid'] for r in summary),issues=issues)
(out/'status.json').write_text(json.dumps(status,indent=2));print(status,flush=True)
page='<html><body><h1>Xee results — '+stamp+'</h1><p>Partial results. Missing or invalid results are not replaced with zeros. Bias intervals use valid toys only; see CSV for exclusions and pending counts. High-mass envelope validation is ongoing.</p>'
for name in ['gof','limits','bias_r0','bias_expected']:
 page+=f'<h2>{name}</h2><a href="{name}.pdf">PDF</a><br><img style="max-width:100%" src="{name}.png">'
page+='<p><a href="gof.csv">GoF CSV</a> · <a href="limits.csv">Limits CSV</a> · <a href="bias.csv">Bias CSV</a></p></body></html>';(out/'index.html').write_text(page)
