import numpy as np,csv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path
p=Path('outputs');eras=[('2022','ch1'),('2022EE','ch2'),('2023','ch3'),('2023BPix','ch4')];comps=[('dy','Chebyshev grado 5'),('upsilon1s','Υ(1S)'),('upsilon2s','Υ(2S)'),('signal','Segnale 9 GeV')];summary=[]
for mode in ['individual','combined']:
 fig,axs=plt.subplots(4,4,figsize=(19,12),layout='constrained',sharex=True);fr,ars=plt.subplots(4,4,figsize=(19,11),layout='constrained',sharex=True)
 for j,(era,ch) in enumerate(eras):
  case=era if mode=='individual' else 'combined';ds={st:np.genfromtxt(p/f'degree5_allFloat_{case}_{ch}_{st}_shapes.csv',delimiter=',',names=True) for st in ['prefit','fit_b','fit_s']};x=ds['prefit']['mass'];ed=np.r_[x-1/150,x[-1]+1/150]
  for i,(key,title) in enumerate(comps):
   ax=axs[i,j];ar=ars[i,j];positive=[];signed=False
   summary.append(dict(fit=mode,era=era,component=key,**{st:float(ds[st][key].sum()) for st in ds}))
   for st,c,ls,lab in [('prefit','#444444','--','Pre-fit'),('fit_b','#1766ba',':','Post-fit B'),('fit_s','#cc2d2d','-','Post-fit S+B')]:
    y=ds[st][key];display=y.copy()
    if np.any(y<0):display=np.abs(y);ls='-.';lab='|Post-fit S+B|, negativo';signed=True
    positive.extend(display[display>0]);ax.stairs(np.where(display>0,display,np.nan),ed,color=c,ls=ls,lw=1.5,label=lab,baseline=None)
    if abs(y.sum())>0:
     base=ds['prefit'][key];ratio=(y/y.sum())/(base/base.sum());ar.plot(x,ratio,color=c,ls=ls,lw=1.5,label='Post-fit S+B' if st=='fit_s' else lab)
   ax.set_yscale('log');ax.set_ylim(min(positive)*.6,max(positive)*1.8);ax.set_title(era+' — '+title,fontsize=10);ax.grid(which='both',alpha=.13);ax.legend(fontsize=6.5,loc='best')
   ar.set_title(era+' — '+title,fontsize=10);ar.grid(alpha=.2);ar.legend(fontsize=6.5)
   if key in ['upsilon1s','upsilon2s']:ar.set_ylim(.98,1.02)
   if key=='signal':
    ax.text(.02,.03,'B: segnale nullo',transform=ax.transAxes,color='#1766ba',fontsize=7)
    if signed:ax.text(.02,.11,'S+B: mostrato valore assoluto',transform=ax.transAxes,color='#cc2d2d',fontsize=7)
   if j==0:ax.set_ylabel('Eventi / 0.01333 GeV');ar.set_ylabel('Shape post/pre\n(integrali a 1)')
   if i==3:ax.set_xlabel('Massa [GeV]');ar.set_xlabel('Massa [GeV]')
   ax.set_xlim(6,12);ar.set_xlim(6,12)
 label='Fit separati per era' if mode=='individual' else 'Fit combinato, r condiviso'
 fig.suptitle(label+' — componenti pre/post, grado5, range6–12, binning originale\nNormalizzazioni del fit; curve negative mostrate in valore assoluto e identificate in legenda',fontsize=13)
 fr.suptitle(label+' — variazione della sola forma\nRapporti post/pre dopo normalizzazione a integrale1; per il segnale negativo si confronta il profilo in valore assoluto',fontsize=13)
 for ext in ['png','pdf']:
  fig.savefig(p/f'degree5_allFloat_all_eras_{mode}_pre_post_components.{ext}',dpi=160);fr.savefig(p/f'degree5_allFloat_all_eras_{mode}_pre_post_shape_ratios.{ext}',dpi=160)
with open(p/'degree5_allFloat_all_eras_pre_post_component_yields.csv','w') as f:
 wr=csv.DictWriter(f,fieldnames=summary[0].keys());wr.writeheader();wr.writerows(summary)
print('Completed component comparisons for individual and combined fits')
