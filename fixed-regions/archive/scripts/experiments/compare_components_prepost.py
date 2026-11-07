import numpy as np,json,csv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path
p=Path('outputs');components=[('dy','Chebyshev'),('upsilon1s','Υ(1S)'),('upsilon2s','Υ(2S)'),('signal','Segnale a 9 GeV')];sets=[];summary=[]
for suffix,title in [('', 'Fit 6–12 GeV'),('_to11','Fit 6–11 GeV')]:
 ds={k:np.genfromtxt(p/f'combine_2022{suffix}_{k}_shapes.csv',delimiter=',',names=True) for k in ['prefit','fit_b','fit_s']};sets.append((title,ds))
 for k,label in components:
  row={'range':title,'component':label}
  for stage in ['prefit','fit_b','fit_s']:row[stage+'_histogram_integral']=float(ds[stage][k].sum())
  for stage in ['fit_b','fit_s']:row[stage+'_post_over_pre_integral']=row[stage+'_histogram_integral']/row['prefit_histogram_integral']
  summary.append(row)
fig,axs=plt.subplots(4,2,figsize=(13,13),layout='constrained',sharex='col')
figr,axr=plt.subplots(4,2,figsize=(13,11),layout='constrained',sharex='col')
for j,(title,ds) in enumerate(sets):
 x=ds['prefit']['mass'];edges=np.r_[x-1/150,x[-1]+1/150]
 for i,(k,label) in enumerate(components):
  ax=axs[i,j];ar=axr[i,j]
  for stage,color,style,lab in [('prefit','#424242','--','Pre-fit'),('fit_b','#1766ba',':','Post-fit B'),('fit_s','#cc2d2d','-','Post-fit S+B')]:
   y=ds[stage][k];ax.stairs(np.where(y>0,y,np.nan),edges,color=color,ls=style,lw=1.8,label=lab,baseline=None)
   if y.sum()>0:
    yn=y/y.sum();pre=ds['prefit'][k];pre=pre/pre.sum();ar.plot(x,yn/pre,color=color,ls=style,lw=1.6,label=lab)
  ax.set_yscale('log');positive=np.concatenate([ds[z][k][ds[z][k]>0] for z in ['prefit','fit_b','fit_s']]);ax.set_ylim(positive.min()*.6,positive.max()*1.7);ax.set_ylabel('Eventi / bin');ax.set_title(label+' — '+title);ax.grid(alpha=.15,which='both');ax.legend(fontsize=8,loc='best')
  ar.set_ylabel('Shape post / pre\n(normalizzate a 1)');ar.set_title(label+' — '+title);ar.grid(alpha=.2);ar.legend(fontsize=8)
  if k in ['upsilon1s','upsilon2s']:ar.set_ylim(.98,1.02)
  if k=='signal':
   ax.text(.02,.04,'Post-fit B: segnale nullo',transform=ax.transAxes,fontsize=9,color='#1766ba');ar.text(.02,.04,'B: rapporto non definito (segnale nullo)',transform=ar.transAxes,fontsize=8,color='#1766ba')
 for ax in [axs[-1,j],axr[-1,j]]:ax.set_xlabel('Massa [GeV]')
 axs[-1,j].set_xlim(6,12 if j==0 else 11);axr[-1,j].set_xlim(6,12 if j==0 else 11)
fig.suptitle('2022 — ogni componente pre-fit e post-fit\nCurve con le normalizzazioni salvate da Combine; asse verticale logaritmico',fontsize=14)
figr.suptitle('2022 — variazione della sola forma\nOgni componente normalizzata a integrale 1 nel proprio range di fit',fontsize=14)
for ext in ['png','pdf']:
 fig.savefig(p/f'componenti_pre_post_due_range.{ext}',dpi=160)
 figr.savefig(p/f'componenti_pre_post_shape_ratio.{ext}',dpi=160)
with (p/'componenti_pre_post_integrali.csv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=summary[0].keys());w.writeheader();w.writerows(summary)
print(json.dumps(summary,indent=2))
