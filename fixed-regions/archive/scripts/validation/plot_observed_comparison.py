from pathlib import Path
import numpy as np,uproot
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
p=Path('/Users/imb/Documents/Codex/Analisi-statistica-Xee/runs/20261029_background_validation')
for mass,limits in [('7.5',(6,8)),('9.5',(8.7,10.3))]:
 files={m:uproot.open(p/f'observed/m{mass}_{m}/fitDiagnostics_validation.root') for m in ['original','degree5']}
 fig=plt.figure(figsize=(14,9));outer=fig.add_gridspec(2,2,hspace=.3,wspace=.2)
 for i,(cell,era) in enumerate(zip(outer,['2022','2022EE','2023','2023BPix']),1):
  gs=cell.subgridspec(2,1,height_ratios=[3,1],hspace=.04);ax=fig.add_subplot(gs[0]);res=fig.add_subplot(gs[1],sharex=ax)
  for model,col in [('original','#246eae'),('degree5','#d35226')]:
   dr=files[model][f'shapes_fit_b/ch{i}'];v,e=dr['total'].to_numpy();w=np.diff(e);x=(e[:-1]+e[1:])/2;m=(x>=limits[0])&(x<limits[1]);ed=np.r_[e[:-1][m],e[1:][m][-1]];data=dr['data'].values('y')*w;B=v*w
   if model=='original':ax.errorbar(x[m],data[m],yerr=np.sqrt(data[m]),fmt='.',ms=2,lw=.5,color='black',label='Dati')
   ax.stairs(B[m],ed,color=col,lw=1.5,label='Fondo '+model,baseline=None);res.plot(x[m],((data-B)/np.sqrt(B))[m],'.',ms=2.5,color=col,label=model)
  ax.set_title(era);ax.set_ylabel('Eventi / 13,33 MeV');ax.tick_params(labelbottom=False);ax.grid(alpha=.15);ax.legend(fontsize=8);ax.set_xlim(*limits)
  res.axhline(0,color='gray',lw=.8);res.axhspan(-2,2,color='gray',alpha=.12);res.set_ylim(-5,5);res.set_ylabel('(D−B)/√B');res.set_xlabel('Massa [GeV]');res.grid(alpha=.15)
 fig.suptitle(f'Confronto B-only: envelope originale e Chebyshev grado 5 — ipotesi {mass} GeV\nFit combinato 6–12 GeV, Υ libere, binning originale',fontsize=14,y=.98);fig.subplots_adjust(top=.9,bottom=.07,left=.07,right=.98)
 for ext in ['png','pdf']:fig.savefig(p/f'observed_comparison_m{mass}.{ext}',dpi=170)
 plt.close(fig)
