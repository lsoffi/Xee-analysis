import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import json
from pathlib import Path
p=Path('outputs');p.mkdir(exist_ok=True)
a=np.genfromtxt('work/xee_2022_nofit.csv',delimiter=',',names=True)
n=a['data'];sig=a['signal_fraction'];bg=a['background_fraction'];N=n.sum();ns=220.687;nb=29701.3
scale=N/(ns*sig.sum()+nb*bg.sum());s=scale*ns*sig;b=scale*nb*bg;y=s+b
x=(a['low']+a['high'])/2;edges=np.r_[a['low'],a['high'][-1]]
fig,axs=plt.subplots(1,2,figsize=(13,4.8),layout='constrained')
for ax in axs:
 ax.errorbar(x,n,yerr=np.sqrt(n),fmt='.',ms=3,color='black',elinewidth=.6,label='Dati 2022')
 ax.stairs(y,edges,color='#d33f36',lw=1.8,label='Somma normalizzata ai dati')
 ax.stairs(b,edges,color='#2171b5',lw=1.5,label='Fondo Chebyshev (grado 7)')
 ax.stairs(s,edges,color='#29934b',lw=1.5,label='Segnale DCB, massa 9, r = 1')
 ax.set_xlabel('Massa [GeV]');ax.set_ylabel('Eventi / 0.01333 GeV');ax.set_ylim(bottom=0);ax.grid(alpha=.15)
axs[0].set_xlim(6,12);axs[1].set_xlim(8,10.5);axs[1].set_ylim(0,max(n[(x>8)&(x<10.5)])*1.15)
axs[0].legend(fontsize=8);fig.suptitle('2022 — DCB + Chebyshev, parametri originali, nessun fit',fontsize=13)
fig.savefig(p/'xee_2022_chebyshev_somma_senza_fit.png',dpi=180)
fig.savefig(p/'xee_2022_chebyshev_somma_senza_fit.pdf')
np.savetxt(p/'xee_2022_chebyshev_somma_senza_fit.csv',np.c_[x,n,s,b,y],delimiter=',',header='mass,data,signal,chebyshev,sum',comments='')
meta=dict(data_total=float(N),signal_fraction_integral=float(sig.sum()),background_fraction_integral=float(bg.sum()),common_scale=float(scale),signal_after_scale=float(s.sum()),background_after_scale=float(b.sum()),total_after_scale=float(y.sum()),minimum_background_bin=float(b.min()),fit_performed=False,era='2022',signal_mass_GeV=9,r=1,components=['DCB signal','Chebyshev continuum degree 7'],excluded_components=['upsilon1s','upsilon2s'])
(p/'xee_2022_chebyshev_somma_senza_fit.json').write_text(json.dumps(meta,indent=2))
print(json.dumps(meta,indent=2))
