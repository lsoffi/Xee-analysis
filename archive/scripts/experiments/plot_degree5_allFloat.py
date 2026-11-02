import numpy as np,json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path
p=Path('outputs');meta={r['era']:r for r in json.load(open(p/'degree5_allFloat_all_eras_summary.json'))};eras=[('2022','ch1'),('2022EE','ch2'),('2023','ch3'),('2023BPix','ch4')]
for which in ['individual','combined']:
 fig,axs=plt.subplots(4,2,figsize=(13,14),layout='constrained',sharex=True)
 for i,(era,ch) in enumerate(eras):
  case=era if which=='individual' else 'combined';m=meta[case]
  for j,(stage,label) in enumerate([('fit_b','Solo fondo'),('fit_s','Segnale + fondo')]):
   ax=axs[i,j];a=np.genfromtxt(p/f'degree5_allFloat_{case}_{ch}_{stage}_shapes.csv',delimiter=',',names=True);x=a['mass'];edges=np.r_[x-1/150,x[-1]+1/150];n=a['data']
   ax.errorbar(x,n,yerr=np.sqrt(n),fmt='.',ms=2.6,color='black',lw=.5,label='Dati')
   for key,color,lab,ls in [('total','#c62828','Totale','-'),('dy','#1565c0','Chebyshev grado5','--'),('upsilon1s','#ef8a00','Υ(1S), norma libera','-.'),('upsilon2s','#9c27b0','Υ(2S)',':'),('signal','#188341','Segnale 9 GeV','-')]:
    y=a[key]
    if key=='signal' and np.any(y<0):y=abs(y);lab='|Segnale|, r negativo';ls='--'
    if np.any(y>0):ax.stairs(np.where(y>0,y,np.nan),edges,color=color,ls=ls,lw=1.4,label=lab,baseline=None)
   ax.set_yscale('log');ax.set_xlim(8,10.5);ax.set_ylim(.001,max(n[(x>=8)&(x<=10.5)])*3);ax.set_title(era+' — '+label,fontsize=11);ax.grid(alpha=.12,which='both');ax.legend(fontsize=7,loc='lower left',ncol=2)
   if j==0:ax.set_ylabel('Eventi / 0.01333 GeV')
   if j==1:ax.text(.03,.96,f'r = {m["r"]:.3f} −{abs(m["r_lo"]):.3f}/+{m["r_hi"]:.3f}',transform=ax.transAxes,va='top',fontsize=8,bbox=dict(facecolor='white',edgecolor='none',alpha=.9))
   if i==3:ax.set_xlabel('Massa [GeV]')
 title='Fit separati per era' if which=='individual' else 'Fit combinato, r condiviso'
 fig.suptitle(title+' — Chebyshev grado5, range6–12GeV, binning originale\nTutti i fit: status0, covQual3; zoom8–10.5GeV',fontsize=13)
 for ext in ['png','pdf']:fig.savefig(p/f'degree5_allFloat_all_eras_{which}_fits.{ext}',dpi=170)
