from pathlib import Path
import json,csv,statistics as st,math
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
web=Path('/eos/user/s/soffi/www/Xee/261106/scan_SR_CR/bias_expected_04_200');web.mkdir(exist_ok=True)
jobs=[l.split() for n in ('pilot','remaining') for l in Path(n+'.txt').read_text().splitlines()]
rows=[];summary=[]
for _,g,_,seed in jobs:
 primary=Path('results')/seed/'result.json';retry=Path('retry_results')/seed/'result.json'
 attempts=[json.loads(f.read_text()) for f in (primary,retry) if f.exists()]
 assert attempts,seed
 d=attempts[-1].copy()
 d['branches']={k:min((a['branches'][k] for a in attempts),key=lambda v:(not v['ok'],float('inf') if v['score'] is None else v['score'])) for k in ('C2','L3')}
 assert d['generator']==g and d['seed']==int(seed)
 b=d['branches'];win=min(b,key=lambda k:b[k]['score']) if all(v['ok'] for v in b.values()) else ''
 row=dict(seed=seed,gen=g,winner=win,r=b[win]['r'] if win else '')
 for k,v in b.items():row.update({k+'_'+key:value for key,value in v.items()})
 rows.append(row)
for g in ('C2','L2','C3','L3'):
 sub=[r for r in rows if r['gen']==g];assert len(sub)==200
 good=[r for r in sub if r['winner']];rs=np.array([r['r'] for r in good])-11.25;sd=rs.std(ddof=1);rng=np.random.default_rng(781);rep=rng.choice(rs,(5000,len(rs)));lo,hi=np.quantile(rep.mean(1)/rep.std(1,ddof=1),[.025,.975])
 summary.append(dict(generator=g,used=len(rs),failed=200-len(rs),cheb2_wins=sum(r['winner']=='C2' for r in good),r_injected=11.25,mean_residual=rs.mean(),sd=sd,mean_error=sd/len(rs)**.5,bias_norm=rs.mean()/sd,CI_low=lo,CI_high=hi))
for name,items in [('toys',rows),('bias',summary)]:
 with (web/(name+'.csv')).open('w') as out:
  wr=csv.DictWriter(out,fieldnames=list(items[0]));wr.writeheader();wr.writerows(items)
fig,ax=plt.subplots(figsize=(7,4))
ax.axhspan(-.2,.2,color='#dcf1e2');ax.axhline(0,color='gray',lw=1)
for i,r in enumerate(summary):
 y=r['bias_norm'];ax.errorbar(i,y,yerr=[[y-r['CI_low']],[r['CI_high']-y]],fmt='o',color='#115c82',capsize=4)
ax.set(xticks=range(4),xticklabels=['Cheb 2','Log 2','Cheb 3','Log 3'],xlabel='Generatore',ylabel='Bias / dispersione dei toy',title='0,4 GeV · r iniettato = 11.2500 · envelope Cheb 2 + Log 3')
ax.margins(x=.15);fig.tight_layout();fig.savefig(web/'bias.png',dpi=150)
print(json.dumps(summary))
