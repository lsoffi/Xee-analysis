import csv
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
root=Path(__file__).resolve().parent
def rows(name):return list(csv.DictReader((root/name).open()))
for reg in ['low','middle','high']:
 f,axes=plt.subplots(1,2,figsize=(10,3),sharey=True)
 for ax,zero in zip(axes,[True,False]):
  if reg=='low':
   data=[r for r in rows('low.csv') if r['model']=='C3+L3' and r['scope']=='common' and r['tag']==('env' if zero else 'expected')]
  elif reg=='middle':data=rows('middle_zero.csv' if zero else 'middle_injected.csv')
  else:data=rows('high_zero.csv') if zero else [r for r in rows('high_snapshot.csv') if r['tag']=='expected']
  for i,g in enumerate(sorted({r.get('generator',r.get('gen')) for r in data})):
   ds=sorted([r for r in data if r.get('generator',r.get('gen'))==g],key=lambda r:float(r['mass']))
   offset=(i-1.5)*.003 if reg=='low' else 0
   ax.plot([float(r['mass'])+offset for r in ds],[float(r['bias']) for r in ds],['o','s','^','D'][i],ms=4,label=g)
  ax.axhspan(-.2,.2,color='#e6efe8',zorder=-2);ax.axhline(0,lw=.6,color='gray');ax.set(title='Background only' if zero else 'Median expected signal injected',xlabel='Signal mass [GeV]');ax.legend(fontsize=8,ncol=2)
 axes[0].set_ylabel('Normalised bias');f.tight_layout();f.savefig(root/(reg+'.png'),dpi=155);plt.close(f)
