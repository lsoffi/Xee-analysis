import numpy as np,json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path
p=Path('outputs');eras=[('2022','ch1'),('2022EE','ch2'),('2023','ch3'),('2023BPix','ch4')];comps=[('dy','Chebyshev'),('upsilon1s','Υ(1S)'),('upsilon2s','Υ(2S)'),('signal','Segnale 9 GeV')]
def read(degree,case,ch,stage):return np.genfromtxt(p/f'degree{degree}_allFloat_{case}_{ch}_{stage}_shapes.csv',delimiter=',',names=True)
for mode in ['combined','individual']:
 for stage in ['fit_b','fit_s']:
  fig,axs=plt.subplots(4,4,figsize=(19,12),layout='constrained',sharex=True)
  for j,(era,ch) in enumerate(eras):
   case='combined' if mode=='combined' else era;arrays={d:read(d,case,ch,stage) for d in [5,7]}
   for i,(key,title) in enumerate(comps):
    ax=axs[i,j];positive=[]
    for d,c,ls in [(5,'#1766ba','--'),(7,'#cc2d2d','-')]:
     a=arrays[d];y=a[key];lab=f'Grado {d}'
     if np.any(y<0):y=abs(y);lab+=', valore assoluto (r < 0)';ls='-.'
     positive.extend(y[y>0]);x=a['mass'];ed=np.r_[x-1/150,x[-1]+1/150]
     if np.any(y>0):ax.stairs(np.where(y>0,y,np.nan),ed,color=c,ls=ls,lw=1.6,label=lab,baseline=None)
     else:ax.plot([],[],color=c,ls=ls,label=lab+': nullo')
    if positive:ax.set_yscale('log');ax.set_ylim(min(positive)*.6,max(positive)*1.8)
    else:ax.text(.5,.5,'Segnale nullo (r = 0)',ha='center',transform=ax.transAxes)
    ax.set_title(era+' — '+title,fontsize=10);ax.set_xlim(6,12);ax.grid(alpha=.13,which='both');ax.legend(fontsize=6.5)
    if j==0:ax.set_ylabel('Eventi / 0.01333 GeV')
    if i==3:ax.set_xlabel('Massa [GeV]')
  model='Fit combinato' if mode=='combined' else 'Fit separati per era';label='solo fondo' if stage=='fit_b' else 'segnale + fondo'
  fig.suptitle(f'{model}, {label} — funzioni post-fit grado 5 / grado 7\nRange 6–12, binning originale; coefficienti e tutte le normalizzazioni del fondo liberi',fontsize=13)
  for ext in ['png','pdf']:fig.savefig(p/f'overlay_degree5_degree7_allFloat_{mode}_{stage}_components.{ext}',dpi=160)
  plt.close(fig)
# Focus on the Chebyshev, including actual pointwise ratio.
fig,axs=plt.subplots(3,4,figsize=(17,10),layout='constrained')
for j,(era,ch) in enumerate(eras):
 arrays={d:read(d,'combined',ch,'fit_s') for d in [5,7]}
 for i,(lo,hi) in enumerate([(6,12),(10.8,12)]):
  ax=axs[i,j];vv=[]
  for d,c,ls in [(5,'#1766ba','--'),(7,'#cc2d2d','-')]:
   a=arrays[d];mask=(a['mass']>=lo)&(a['mass']<=hi);ax.plot(a['mass'][mask],a['dy'][mask],color=c,ls=ls,lw=1.6,label=f'Grado {d}');vv.extend(a['dy'][mask])
  ax.set_yscale('log');ax.set_ylim(min(vv)*.8,max(vv)*1.25);ax.set_xlim(lo,hi);ax.set_title(era+(' — Chebyshev S+B' if i==0 else ' — zoom coda'));ax.set_ylabel('Eventi / bin');ax.grid(alpha=.15,which='both');ax.legend(fontsize=8)
 a,b=arrays[5],arrays[7];axs[2,j].plot(a['mass'],b['dy']/a['dy'],color='#6a3d9a',lw=1.5);axs[2,j].axhline(1,color='gray',ls='--',lw=.8);axs[2,j].set_xlim(6,12);axs[2,j].set_ylabel('Post grado 7 / grado 5');axs[2,j].set_xlabel('Massa [GeV]');axs[2,j].grid(alpha=.2)
fig.suptitle('Fit combinato S+B — confronto del Chebyshev post-fit, tutto floating\nCurve con le normalizzazioni del fit; rapporti non normalizzati a integrale 1',fontsize=13)
for ext in ['png','pdf']:fig.savefig(p/f'overlay_degree5_degree7_allFloat_chebyshev_focus.{ext}',dpi=170)
# Short result report.
ms={d:json.load(open(p/f'degree{d}_allFloat_all_eras_summary.json')) for d in [5,7]};lines=['# Post-fit grado5 / grado7, tutto floating','','Range6–12GeV,450 bin per era. Floatano coefficienti, normalizzazioni dei3 fondi, r e nuisance originali; coefficienti MC delle DCB restano costanti. Stesse opzioni Combine: strategia0,robustFit1,r tra−30 e30, nessun freezeParameters.','', '| Era | r grado5 (68%) | r grado7 (68%) | covQual B/S+B grado5 | covQual B/S+B grado7 |','|---|---:|---:|---|---|']
for a,b in zip(ms[5],ms[7]):
 def fmt(r):return f"{r['r']:.5f} −{abs(r['r_lo']):.4f}/+{r['r_hi']:.4f}"
 lines.append(f"| {a['era']} | {fmt(a)} | {fmt(b)} | {a['covQual_b']}/{a['covQual_s']} | {b['covQual_b']}/{b['covQual_s']} |")
lines+=['','Curve con normalizzazioni post-fit, non riscalate per farle coincidere. Componenti negative: si mostra il valore assoluto, con legenda esplicita. Il pre-fit usa seed numerici positivi comuni per inizializzare i polinomi, non vincoli aggiuntivi. Il confronto di modelli non valida la scelta del grado.']
(p/'overlay_degree5_degree7_allFloat_report.md').write_text('\n'.join(lines)+'\n')
print('Completed overlays')
