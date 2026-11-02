from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
root=Path('/Users/imb/Documents/Codex/Analisi-statistica-Xee/runs/20261029_original_mass_scan')
out=root/'plots/data_zoom'
fig=plt.figure(figsize=(15,10))
outer=fig.add_gridspec(2,2,hspace=.3,wspace=.2)
for cell,era in zip(outer,['2022','2022EE','2023','2023BPix']):
    gs=cell.subgridspec(2,1,height_ratios=[3,1],hspace=.06)
    ax=fig.add_subplot(gs[0]);res=fig.add_subplot(gs[1],sharex=ax)
    src=root/f'fits/mass_7.5/fitdiagnostics/{era}/Xee_ee_4_{era}_fit_b_shapes.csv'
    a=np.genfromtxt(src,delimiter=',',names=True)
    a=a[(a['low']>=6-1e-9)&(a['high']<=8+1e-9)]
    x=a['mass_GeV'];y=a['data'];b=a['total'];edges=np.r_[a['low'],a['high'][-1]]
    ax.errorbar(x,y,yerr=np.sqrt(y),fmt='.',ms=2.8,lw=.6,color='black',label='Dati')
    ax.stairs(b,edges,color='#d14d24',lw=1.5,label='Fondo totale, fit B-only')
    ax.stairs(a['dy'],edges,color='#2879ad',lw=1,ls='--',label='Continuo')
    ax.set_title(era);ax.set_ylabel('Eventi / 13,33 MeV');ax.set_ylim(bottom=0);ax.grid(alpha=.15)
    ax.legend(fontsize=8,loc='upper right');ax.tick_params(labelbottom=False)
    pull=(y-b)/np.sqrt(b)
    res.axhspan(-2,2,color='gray',alpha=.12)
    res.axhline(0,color='#d14d24',lw=1)
    res.plot(x,pull,'.',ms=3,color='black');res.set_ylim(-5,5);res.set_yticks([-4,-2,0,2,4])
    res.set_ylabel('(D−B)/√B',fontsize=10);res.set_xlabel('Massa dielettronica [GeV]');res.set_xlim(6,8);res.grid(alpha=.15)
    np.savetxt(out/f'residuals_{era}_6_8.csv',np.c_[x,y,b,pull],delimiter=',',header='mass_GeV,data,background,pearson_residual',comments='')
fig.suptitle('Zoom 6–8 GeV: dati e fondo post-fit per era',fontsize=17,y=.98)
fig.text(.5,.943,'Fit B-only originale su 6–12 GeV, ipotesi 7,5 GeV; Υ(1S) congelata. Residui senza incertezza del modello.',ha='center',fontsize=10)
fig.subplots_adjust(top=.90,bottom=.07,left=.07,right=.98)
for ext in ['png','pdf']:fig.savefig(out/f'data_background_residuals_four_eras_6_8.{ext}',dpi=180)
