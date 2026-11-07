import csv,json,math
from pathlib import Path
import numpy as np,ROOT
from scipy.integrate import quad
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT.gSystem.Load('libHiggsAnalysisCombinedLimit')
d=Path.cwd();orig=d.parent/'bias_phi_jpsi_quartic';out=Path('/eos/user/s/soffi/www/Xee/261106/scan_SR_CR/bias_phi_jpsi_quartic_r0')
geom={r['mass']:r for r in json.loads((orig/'geometry.json').read_text())}
seeds=set();groups={};rows=[];edges=[]
for line in (d/'all.txt').read_text().splitlines():
 m,g,inj,seed=line.split();assert float(inj)==json.loads((d/"medians.json").read_text())[m] and seed not in seeds and g in ('C4','L4');seeds.add(seed)
 q=json.loads((d/('m'+m)/'results'/(seed+'.json')).read_text());assert q['seed']==int(seed) and q['generator']==g
 b=q['branches'];assert set(b)=={'C4','L4'}
 for v in b.values():
  if v['ok']:assert v['status']==0 and abs(v['r'])<79.999 and math.isfinite(v['D']) and abs(v['score']-v['D']-16)<1e-7
 win=min(b,key=lambda k:b[k]['score']) if all(v['ok'] for v in b.values()) else ''
 row=dict(mass=float(m),gen=g,seed=seed,winner=win,r=b[win]['r'] if win else '')
 for k,v in b.items():row.update({k+'_'+a:h for a,h in v.items()})
 rows.append(row);groups.setdefault((float(m),g),[]).append(row)
assert len(seeds)==6400 and len(groups)==32
summary=[];rng=np.random.default_rng(781)
for (m,g),rr in sorted(groups.items()):
 assert len(rr)==200
 keep=[r for r in rr if r['winner']];injected=json.loads((d/'medians.json').read_text())[str(m)];x=np.array([r['r'] for r in keep])-injected;sd=x.std(ddof=1);rep=rng.choice(x,(5000,len(x)));lo,hi=np.quantile(rep.mean(1)/rep.std(1,ddof=1),[.025,.975])
 summary.append(dict(mass=m,gen=g,valid=len(x),failed=200-len(x),C4wins=sum(r['winner']=='C4' for r in keep),L4wins=sum(r['winner']=='L4' for r in keep),injection=injected,mean_residual=x.mean(),sd=sd,bias=x.mean()/sd,lo=lo,hi=hi))
 p=orig/f'm{m:.1f}';md=json.loads((p/(g+'.json')).read_text());r=geom[m]
 assert md['geometry']['sr_segments']==r['sr_segments'] and md['geometry']['selected_CRs']==r['selected_CRs']
 contacts=[(c['group'],x) for c in md['eras'][0]['CRs'] for x in c['range'] if any(abs(y-x)<1e-9 for seg in r['sr_segments'] for y in seg)]
 if contacts:
  f=ROOT.TFile.Open(str(p/f'{g}_{m:.1f}.root'));w=f.Get('w');w.loadSnapshot('background')
 for em in md['eras']:
  assert [(a['group'],a['edge']) for a in em['C1_contacts']]==contacts
  for group,e in contacts:
   era=em['era'];get=lambda n:w.arg(n).getVal()
   c=next(c for c in em['CRs'] if c['group']==group);a,b=c['range'];k=c['slope'];cr=get(c['continuum_norm'])*k*math.exp(k*(e-a))/math.expm1(k*(b-a))
   a,b=em['SR_reference'];coeff=[get(n) for n in em['continuum_coefficients']];sr=np.polynomial.chebyshev.chebval(2*(e-a)/(b-a)-1,coeff) if g=='C4' else coeff[0]
   if group=='omega_phi':
    for n in ('jpsi','psi2s'):
     mu=get('mu_'+n+'_'+era);sig,al,nl,ar,nr=[get(n+'_'+j+'_'+era) for j in ('sigma','alphaL','nL','alphaR','nR')]
     def density(x):
      t=(x-mu)/sig
      if t < -al:return math.exp(-al**2/2-nl*math.log1p(al*(-al-t)/nl))
      if t > ar:return math.exp(-ar**2/2-nr*math.log1p(ar*(t-ar)/nr))
      return math.exp(-t*t/2)
     sr+=get('N_'+n+'_'+era)*density(e)/quad(density,2.8,4,points=[v for v in (mu-al*sig,mu+ar*sig) if 2.8<v<4])[0]
   edges.append(abs(sr-cr)/cr)
 if contacts:f.Close()
assert len(edges)==32 and max(edges)<1e-6,max(edges)
out.mkdir(exist_ok=True)
for n,data in [('toys',rows),('bias',summary)]:
 with (out/(n+'.csv')).open('w') as f:
  w=csv.DictWriter(f,fieldnames=data[0]);w.writeheader();w.writerows(data)
fig,ax=plt.subplots(figsize=(10,5));ax.axhspan(-.2,.2,color='#deefde');ax.axhline(0,color='.5',lw=.7)
for g,delta in [('C4',-.007),('L4',.007)]:
 rr=[r for r in summary if r['gen']==g];y=np.array([r['bias'] for r in rr]);ax.errorbar(np.array([r['mass'] for r in rr])+delta,y,yerr=[y-[r['lo'] for r in rr],[r['hi'] for r in rr]-y],marker='o',capsize=2,label=g)
ax.set(xlabel='Mass [GeV]',ylabel='Mean (fitted - injected r) / toy standard deviation',title='C4 + L4 envelope, zero signal injected');ax.legend(title='Generator');fig.tight_layout()
for ext in ['png','pdf']:fig.savefig(out/('bias_vs_mass.'+ext))
print('VALID',sum(r['valid'] for r in summary),'EDGE',max(edges))
for g in ('C4','L4'):
 rr=[r for r in summary if r['gen']==g];print(g,'outside',[(r['mass'],round(r['bias'],3),r['valid']) for r in rr if abs(r['bias'])>.2],'range',min(r['valid'] for r in rr),max(r['valid'] for r in rr),'worst',max(rr,key=lambda r:abs(r['bias'])))
