from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
root=Path('/Users/imb/Documents/Codex/Analisi-statistica-Xee/runs/20261029_original_mass_scan')
out=root/'plots/data_zoom'
fig,axs=plt.subplots(2,2,figsize=(14,8),sharex=True)
for ax,era in zip(axs.flat,['2022','2022EE','2023','2023BPix']):
    src=root/f'fits/mass_7.5/fitdiagnostics/{era}/Xee_ee_4_{era}_fit_b_shapes.csv'
    a=np.genfromtxt(src,delimiter=',',names=True)
    a=a[(a['low']>=6-1e-9)&(a['high']<=8+1e-9)]
    y=a['data']; x=a['mass_GeV']
    ax.stairs(y,np.r_[a['low'],a['high'][-1]],color='#245e96',lw=.8)
    ax.errorbar(x,y,yerr=np.sqrt(y),fmt='.',ms=2.6,lw=.6,color='#245e96',capsize=0)
    ax.set_title(f'{era} — {y.sum():,.0f} eventi tra 6 e 8 GeV'.replace(',',' '))
    ax.set_xlim(6,8);ax.set_ylim(bottom=0);ax.grid(alpha=.2)
    ax.set_ylabel('Eventi / 13,33 MeV');ax.set_xlabel('Massa dielettronica [GeV]')
    np.savetxt(out/f'data_{era}_6_8.csv',np.c_[a['low'],a['high'],y],delimiter=',',header='low_GeV,high_GeV,data',comments='')
fig.suptitle('Spettri osservati — binning originale, nessun fit',fontsize=16)
fig.tight_layout(rect=[0,0,1,.95])
for ext in ['png','pdf']:fig.savefig(out/f'data_four_eras_6_8.{ext}',dpi=180)
