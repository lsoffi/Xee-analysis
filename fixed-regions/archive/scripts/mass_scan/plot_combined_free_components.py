from pathlib import Path
import numpy as np,uproot
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
p=Path('/Users/imb/Documents/Codex/Analisi-statistica-Xee/runs/20261029_m7p5_combined_upsilon1s_free')
f=uproot.open(p/'fitDiagnostics_m7p5_upsilon1s_free.root')
fig=plt.figure(figsize=(16,11));outer=fig.add_gridspec(2,2,hspace=.32,wspace=.2)
for i,(cell,era) in enumerate(zip(outer,['2022','2022EE','2023','2023BPix']),1):
 gs=cell.subgridspec(2,1,height_ratios=[3,1],hspace=.06);ax=fig.add_subplot(gs[0]);res=fig.add_subplot(gs[1],sharex=ax)
 b=f[f'shapes_fit_b/ch{i}'];s=f[f'shapes_fit_s/ch{i}'];v,e=s['total'].to_numpy();w=np.diff(e);x=(e[:-1]+e[1:])/2;mask=(x>=6)&(x<8);edges=np.r_[e[:-1][mask],e[1:][mask][-1]]
 data=s['data'].values('y')*w;B=b['total'].values()*w;S=v*w
 ax.errorbar(x[mask],data[mask],yerr=np.sqrt(data[mask]),fmt='.',ms=2.5,lw=.5,color='black',label='Dati',zorder=6)
 ax.stairs(B[mask],edges,color='#6b6b6b',ls='--',lw=1.6,label='Totale B-only',baseline=None)
 for key,color,label,ls,lw in [('total','#d62728','Totale S+B','-',1.8),('dy','#246eae','Continuo (S+B)','--',1.2),('upsilon1s','#ef8a00','Υ(1S) (S+B)','-.',1.3),('upsilon2s','#9b40a4','Υ(2S) (S+B)',':',1.4),('Zd','#159447','Segnale 7,5 GeV','-',1.6)]:
  y=s[key].values()*w;ax.stairs(y[mask],edges,color=color,ls=ls,lw=lw,label=label,baseline=None)
 ax.set_title(era);ax.set_xlim(6,8);ax.set_yscale('log');ax.set_ylim(.01,max(data[mask])*2);ax.set_ylabel('Eventi / 13,33 MeV');ax.grid(alpha=.15,which='both');ax.tick_params(labelbottom=False);ax.legend(fontsize=8,ncol=2,loc='lower left')
 res.axhspan(-2,2,color='gray',alpha=.12);res.axhline(0,color='gray',lw=.7)
 res.plot(x[mask],((data-B)/np.sqrt(B))[mask],'.',color='#666666',ms=3,label='B-only')
 res.plot(x[mask],((data-S)/np.sqrt(S))[mask],'.',color='#d62728',ms=3,label='S+B')
 res.set_ylim(-5,5);res.set_ylabel('(D−M)/√M');res.set_xlabel('Massa dielettronica [GeV]');res.grid(alpha=.15);res.legend(fontsize=7,ncol=2,loc='lower left')
fig.suptitle('Fit combinato: B-only, S+B e componenti — zoom 6–8 GeV',fontsize=17,y=.98)
fig.text(.5,.948,'Fit su 6–12 GeV; Υ(1S) libera; binning originale. B e S+B: status 0, covQual 3. Residui senza incertezza del modello.',ha='center',fontsize=10)
fig.subplots_adjust(top=.91,bottom=.065,left=.065,right=.98)
for ext in ['png','pdf']:fig.savefig(p/f'combined_B_SB_components_6_8.{ext}',dpi=180)
