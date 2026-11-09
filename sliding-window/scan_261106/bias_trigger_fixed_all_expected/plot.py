import csv,json,math
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
base=Path(__file__).resolve().parent.parent
out=Path('/eos/user/s/soffi/www/Xee/261106/scan_SR_CR/bias_trigger_fixed_all_expected')
out.mkdir(parents=True,exist_ok=True)
groups={};seeds=set();rows=[]
for directory in ('bias_trigger_fixed_all_expected','bias_trigger_fixed_m14'):
 d=base/directory
 for line in (d/'rows.txt').read_text().splitlines():
  m,g,inj,s=line.split();m=float(m);inj=float(inj);s=int(s)
  assert s not in seeds;seeds.add(s)
  q=json.loads((d/'results'/f'{s}.json').read_text())
  assert q['seed']==s and q['generator']==g
  b=q['branches'];assert set(b)=={'C4','L4'}
  for v in b.values():
   if v['ok']:
    assert v['status']==0 and abs(v['r'])<79.999 and math.isfinite(v['D'])
    assert abs(v['score']-v['D']-16)<1e-7
  win=min(b,key=lambda k:b[k]['score']) if all(v['ok'] for v in b.values()) else ''
  row=dict(mass=m,generator=g,injection=inj,seed=s,winner=win,rhat=b[win]['r'] if win else '',valid=bool(win))
  rows.append(row);groups.setdefault((m,g),[]).append(row)
assert len(seeds)==6400 and len(groups)==32
rng=np.random.default_rng(781);summary=[]
for (m,g),rr in sorted(groups.items()):
 assert len(rr)==200 and len({r['injection'] for r in rr})==1
 keep=[r for r in rr if r['valid']];x=np.array([r['rhat']-r['injection'] for r in keep]);assert len(x)>1
 sd=x.std(ddof=1);rep=rng.choice(x,(5000,len(x)));lo,hi=np.quantile(rep.mean(1)/rep.std(1,ddof=1),[.025,.975])
 summary.append(dict(mass=m,generator=g,injection=rr[0]['injection'],valid=len(x),excluded=200-len(x),C4wins=sum(r['winner']=='C4' for r in keep),L4wins=sum(r['winner']=='L4' for r in keep),mean_residual=float(x.mean()),sd=float(sd),bias=float(x.mean()/sd),lo=float(lo),hi=float(hi)))
for name,data in [('toys',rows),('bias',summary)]:
 with (out/f'{name}.csv').open('w') as f:
  w=csv.DictWriter(f,fieldnames=data[0]);w.writeheader();w.writerows(data)
fig,(ax,counts)=plt.subplots(2,1,figsize=(11,7),sharex=True,gridspec_kw={'height_ratios':[3,1]})
ax.axhspan(-.2,.2,color='#deefde',label='Reference |bias| < 0.2');ax.axhline(0,color='.5',lw=.7)
for g,offset,color in [('C4',-.008,'#1766ad'),('L4',.008,'#d26b19')]:
 rr=[r for r in summary if r['generator']==g];xx=np.array([r['mass'] for r in rr])+offset;y=np.array([r['bias'] for r in rr])
 ax.plot(xx,y,'o-',color=color,label=f'{g} generator')
 ax.vlines(xx,[r['lo'] for r in rr],[r['hi'] for r in rr],color=color)
 counts.plot(xx,[r['valid'] for r in rr],'o-',color=color)
ax.set_ylabel(r'$\mathrm{mean}(\hat r-r_{inj})/\mathrm{SD}(\hat r-r_{inj})$')
ax.set_title('C4 + L4 envelope — median expected signal injected\nTrigger nuisance fixed at 0; 95% bootstrap intervals')
ax.legend();counts.set(xlabel='Signal mass [GeV]',ylabel='Valid toys / 200',xticks=np.arange(1.2,2.8,.1));counts.grid(alpha=.25)
fig.tight_layout()
for ext in ('png','pdf'):fig.savefig(out/f'bias_vs_mass.{ext}',dpi=170)
(out/'index.html').write_text('<!doctype html><html lang="en"><meta charset="utf-8"><title>Bias with fixed trigger</title><h1>C4 + L4 envelope with fixed trigger</h1><p>Median expected signal injected. Both fit branches must pass the existing validity criterion. The smaller D+16 selects the envelope branch. Bias is the mean residual divided by its sample standard deviation. Intervals use 5000 bootstrap resamples at 95% confidence.</p><img style="max-width:100%" src="bias_vs_mass.png"><p><a href="bias_vs_mass.pdf">PDF</a> · <a href="bias.csv">Summary CSV</a> · <a href="toys.csv">Toy validity CSV</a></p></html>')
print('VALID',sum(r['valid'] for r in summary),'TOTAL',len(rows))
print('OUTSIDE',[(r['mass'],r['generator'],round(r['bias'],3)) for r in summary if abs(r['bias'])>.2])
print('OUTPUT',out)
