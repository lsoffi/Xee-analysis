import pathlib,json,numpy as np,csv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import norm
from matplotlib.ticker import FuncFormatter
labels={'original':'Envelope originale','degree5':'Chebyshev grado 5','cheb7':'Chebyshev grado 7','bernstein':'Bernstein grado 7'}
P=pathlib.Path('/Users/imb/Documents/Codex/Analisi-statistica-Xee/runs/20261029_background_validation');summary=json.loads((P/'toy_summary.json').read_text());rows=list(csv.DictReader((P/'toys_all.csv').open()))
fig,axs=plt.subplots(2,2,figsize=(12,8),layout='constrained')
for ax,(m,g) in zip(axs.flat,[(m,g) for m in [7.5,9.5] for g in ['cheb7','bernstein']]):
 z=np.linspace(.01,3.5,100);ax.plot(z,norm.sf(z),'k--',label='Atteso teorico asintotico')
 for model,color in [('original','#246eae'),('degree5','#d35226')]:
  v=np.array([float(r['Z']) for r in rows if r.get('valid_point')=='True' and float(r['mass'])==m and r['generator']==g and float(r['injection'])==0 and r['model']==model]);n=len(v)
  if n:
   frac=np.array([(v>=zz).mean() for zz in z]);ax.step(z,np.where(frac>0,frac,np.nan),where='post',color=color,label=f'{labels[model]}, n={n}');ax.axhline(1/n,color=color,alpha=.2,lw=.7)
 ax.set_yscale('log');ax.set_yticks([.01,.1,1]);ax.yaxis.set_major_formatter(FuncFormatter(lambda v,p:f'{v:g}'.replace('.',',')));ax.set_ylim(.003,1);ax.set_xlim(0,3.5);ax.set_title(f'm={m} GeV — fondo generatore: {labels[g]}');ax.set_xlabel('Soglia Z locale');ax.set_ylabel('Frazione di toy B-only con Z ≥ soglia');ax.grid(alpha=.15,which='both');ax.legend(fontsize=8)
fig.suptitle('Toy senza segnale: quanto spesso compare un eccesso apparente?\nn minimi utilizzabili in legenda; Z da q0 di FitDiagnostics, distinta dalla qualità degli intervalli.',fontsize=13)
for ext in ['png','pdf']:fig.savefig(P/f'null_Z_calibration.{ext}',dpi=170)
plt.close(fig)
fig,axs=plt.subplots(1,2,figsize=(14,6),layout='constrained')
for ax,injname in zip(axs,['zero','signal']):
 ticklabels=[]
 for j,(m,g) in enumerate([(m,g) for m in [7.5,9.5] for g in ['cheb7','bernstein']]):
  ticklabels.append(f'{m} GeV\n{labels[g]}')
  for model,col,offset in [('original','#246eae',-.1),('degree5','#d35226',.1)]:
   v=[r for r in summary if r['mass']==m and r['generator']==g and r['model']==model and r['injection']==(0 if injname=='zero' else (.06 if m==7.5 else 3.))]
   if not v or not v[0].get('n_valid'):continue
   r=v[0];y=r['coverage68'];ax.errorbar(j+offset,y,yerr=[[y-r['coverage68_lo']],[r['coverage68_hi']-y]],fmt='o',color=col,label=labels[model] if j==0 else None,capsize=3)
   ax.annotate(f"{r['n_valid']}/{r['n_total']}",(j+offset,y),xytext=(0,10 if offset<0 else -18),textcoords='offset points',ha='center',fontsize=8,color=col)
 ax.axhline(.68,color='black',ls='--',label='Nominale 68%');ax.set_ylim(0,1.05);ax.set_xticks(range(4),ticklabels);ax.set_ylabel('Copertura intervallo 68% tra fit utilizzabili');ax.set_title('Senza segnale' if injname=='zero' else 'Segnale iniettato: r=0,06 / r=3');ax.grid(alpha=.15);ax.legend(fontsize=8)
fig.suptitle('Copertura nella campagna principale — MINOS standard\nEtichette: fit utilizzabili / toy totali; barre: intervalli binomiali Wilson ≈68%',fontsize=13)
for ext in ['png','pdf']:fig.savefig(P/f'coverage68.{ext}',dpi=170)
