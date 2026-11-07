import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path
p=Path('outputs');fig,axs=plt.subplots(4,4,figsize=(19,12),layout='constrained',sharex=True)
for j,(era,ch) in enumerate([('2022','ch1'),('2022EE','ch2'),('2023','ch3'),('2023BPix','ch4')]):
 ds={s:np.genfromtxt(p/f'degree5_sequence_{era}_{ch}_{s}_shapes.csv',delimiter=',',names=True) for s in ['fit_b','fit_s']}
 for i,(key,label) in enumerate([('dy','Chebyshev'),('upsilon1s','Υ(1S)'),('upsilon2s','Υ(2S)'),('signal','Segnale9GeV')]):
  ax=axs[i,j];pos=[]
  for s,c,ls,lab in [('fit_b','#1766ba','--','Fit B = inizio S+B'),('fit_s','#cc2d2d','-','Post S+B')]:
   a=ds[s];y=a[key]
   if np.any(y<0):y=abs(y);lab='|Post S+B|, negativo';ls='-.'
   if np.any(y>0):pos.extend(y[y>0]);x=a['mass'];ax.stairs(np.where(y>0,y,np.nan),np.r_[x-1/150,x[-1]+1/150],color=c,lw=1.5,ls=ls,label=lab,baseline=None)
  ax.set_yscale('log');ax.set_ylim(min(pos)*.6,max(pos)*1.7);ax.set_xlim(6,12);ax.set_title(era+' — '+label,fontsize=10);ax.grid(which='both',alpha=.15);ax.legend(fontsize=7)
  if j==0:ax.set_ylabel('Eventi / 0.01333 GeV')
  if i==3:ax.set_xlabel('Massa [GeV]')
fig.suptitle('Grado5 — inizializzazione dai fit solo fondo della stessa era\nConfronto funzioni del fit B e post S+B; tutti i parametri del fondo restano liberi',fontsize=13)
for ext in ['png','pdf']:fig.savefig(p/f'degree5_sequence_B_to_SB_functions.{ext}',dpi=160)
