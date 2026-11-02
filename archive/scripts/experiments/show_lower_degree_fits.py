import numpy as np,json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path
p=Path('outputs');refwidth=8/75
for upper,stem in [(12,'rebin8'),(11,'to11_rebin8')]:
 fig,axs=plt.subplots(3,2,figsize=(13,12),layout='constrained',sharex=True,sharey=True)
 for i,degree in enumerate([6,5,4]):
  tag=f'{stem}_degree{degree}';m=json.load(open(p/f'combine_2022_{tag}_pre_post.json'));r=m['parameters'][0]
  for j,(stage,title) in enumerate([('fit_b','Solo fondo'),('fit_s','Segnale + fondo')]):
   ax=axs[i,j];a=np.genfromtxt(p/f'combine_2022_{tag}_{stage}_shapes.csv',delimiter=',',names=True);ed=np.r_[a['low'],a['high'][-1]];f=refwidth/a['width'];n=a['data']
   ax.errorbar(a['mass'],n*f,yerr=np.sqrt(n)*f,fmt='.',ms=4,color='black',lw=.7,label='Dati 2022',zorder=4)
   for key,color,label,ls in [('total','#c62828','Totale','-'),('dy','#1565c0','Chebyshev','--'),('upsilon1s','#ef8a00','Υ(1S), norma fissa','-.'),('upsilon2s','#9c27b0','Υ(2S)',':'),('signal','#188341','Segnale 9 GeV','-')]:
    y=a[key]*f
    if key=='signal':
     if np.all(y==0):continue
     if np.any(y<0):y=abs(y);label='|Segnale|, r negativo';ls='--'
    ax.stairs(np.where(y>0,y,np.nan),ed,color=color,ls=ls,lw=1.5,label=label,baseline=None)
   q=m['b' if stage=='fit_b' else 's'];ax.set_title(f'Grado {degree} — {title}; covQual {q["covQual"]}',fontsize=11)
   if stage=='fit_s':ax.text(.03,.92,f'r = {r["s"]:.3f}  −{abs(r["s_lo"]):.3f}/+{r["s_hi"]:.3f}',transform=ax.transAxes,fontsize=9,va='top',bbox=dict(facecolor='white',edgecolor='none',alpha=.9))
   else:ax.text(.03,.92,'r = 0 (fisso)',transform=ax.transAxes,fontsize=9,va='top',bbox=dict(facecolor='white',edgecolor='none',alpha=.9))
   ax.set_yscale('log');ax.set_ylim(.003,2200);ax.set_xlim(8,10.5);ax.grid(which='both',alpha=.13);ax.legend(fontsize=7,loc='lower left',ncol=2)
   if j==0:ax.set_ylabel('Eventi / 0.10667 GeV')
   if i==2:ax.set_xlabel('Massa [GeV]')
 fig.suptitle(f'Fit 2022 — range 6–{upper} GeV, rebin ×8\nGradi ridotti del Chebyshev; zoom 8–10.5 GeV, tutte le componenti post-fit',fontsize=14)
 for ext in ['png','pdf']:fig.savefig(p/f'fit_gradi_ridotti_6_{upper}_rebin8.{ext}',dpi=170)
