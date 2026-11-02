from pathlib import Path
import numpy as np,uproot
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
p=Path('/Users/imb/Documents/Codex/Analisi-statistica-Xee/runs/20261029_m7p5_combined_upsilon1s_free');f=uproot.open(p/'fitDiagnostics_m7p5_upsilon1s_free.root')
fig,axs=plt.subplots(2,2,figsize=(14,8),sharex=True,layout='constrained')
for i,(ax,era) in enumerate(zip(axs.flat,['2022','2022EE','2023','2023BPix']),1):
 b=f[f'shapes_fit_b/ch{i}'];s=f[f'shapes_fit_s/ch{i}'];v,e=s['total'].to_numpy();w=np.diff(e);x=(e[:-1]+e[1:])/2;m=(x>=6)&(x<8);edges=np.r_[e[:-1][m],e[1:][m][-1]]
 signal=s['Zd'].values()*w;delta=(s['total'].values()-s['Zd'].values()-b['total'].values())*w;dy=(s['dy'].values()-b['dy'].values())*w
 for y,col,label,ls in [(signal,'#159447','Segnale S','-'),(delta,'#246eae','ΔB = B(S+B) − B(B-only)','-'),(signal+delta,'#d62728','Δtotale = S + ΔB','-'),(dy,'#ac68a0','Δcontinuo','--')]:ax.stairs(y[m],edges,color=col,label=label,ls=ls,lw=1.5,baseline=None)
 ax.axhline(0,color='black',lw=.7);ax.set_title(era);ax.set_xlim(6,8);ax.set_xlabel('Massa dielettronica [GeV]');ax.set_ylabel('Eventi / 13,33 MeV');ax.grid(alpha=.2);ax.legend(fontsize=8)
 np.savetxt(p/f'compensation_{era}.csv',np.c_[x,signal,delta,signal+delta,dy],delimiter=',',header='mass_GeV,signal,delta_background,delta_total,delta_continuum',comments='')
fig.suptitle('Compensazione tra segnale e fondo — fit combinato, Υ(1S) libera\nFit originale su 6–12 GeV; ipotesi di segnale 7,5 GeV',fontsize=14)
for ext in ['png','pdf']:fig.savefig(p/f'signal_background_compensation_6_8.{ext}',dpi=180)
