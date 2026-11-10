from pathlib import Path
import csv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
data=Path(__file__).resolve().parent/'data/high_partial'
lim=[{k:float(v) for k,v in r.items()} for r in csv.DictReader((data/'limits.csv').open()) if float(r['mass'])>=4.1]
fig,ax=plt.subplots(figsize=(11,5));chunks=[]
for r in lim:
 if not chunks or r['mass']-chunks[-1][-1]['mass']>.11:chunks.append([])
 chunks[-1].append(r)
for i,rr in enumerate(chunks):
 x=[r['mass'] for r in rr]
 ax.fill_between(x,[r['minus2'] for r in rr],[r['plus2'] for r in rr],color='#f6cc43',label='Expected +/-2 sigma' if i==0 else None)
 ax.fill_between(x,[r['minus1'] for r in rr],[r['plus1'] for r in rr],color='#69bb63',label='Expected +/-1 sigma' if i==0 else None)
 ax.plot(x,[r['median'] for r in rr],'k--',label='Median expected' if i==0 else None)
 ax.plot(x,[r['observed'] for r in rr],'ko-',ms=3,label='Observed' if i==0 else None)
ax.set(yscale='log',xlabel='Mass hypothesis [GeV]',ylabel='95% CL upper limit on r',title=f'PRELIMINARY high SR: C4+L4 | {len(lim)}/50 masses')
ax.grid(alpha=.2);ax.legend();fig.tight_layout();fig.savefig(data/'limits_high.png',dpi=180);plt.close(fig)
