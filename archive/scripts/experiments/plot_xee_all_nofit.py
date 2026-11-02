import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import json
from pathlib import Path
p=Path('outputs');a=np.genfromtxt('work/xee_2022_nofit.csv',delimiter=',',names=True)
meta={}
for line in Path('work/xee_2022_nofit_meta.txt').read_text().splitlines():
 k,v=line.rsplit(' ',1);meta[k]=float(v)
n=a['data'];names=['Zd','dy','upsilon1s','upsilon2s'];cols=['signal_fraction','background_fraction','upsilon1s_fraction','upsilon2s_fraction'];raw=np.array([meta[k]*a[c] for k,c in zip(names,cols)])
scale=n.sum()/raw.sum();s,b,u1,u2=raw*scale;y=s+b+u1+u2;btotal=b+u1+u2
x=(a['low']+a['high'])/2;edges=np.r_[a['low'],a['high'][-1]]
fig,axs=plt.subplots(1,2,figsize=(13,5.1),layout='constrained')
for ax in axs:
 ax.errorbar(x,n,yerr=np.sqrt(n),fmt='.',ms=3,color='black',elinewidth=.6,label='Dati 2022')
 for arr,color,label,ls in [(y,'#d33f36','Somma S + fondi','-'),(btotal,'#6a51a3','Fondo totale: Chebyshev + Υ(1S) + Υ(2S)','--'),(b,'#2171b5','Chebyshev (grado 7)','-'),(u1,'#f39c12','Υ(1S)','-'),(u2,'#ad579f','Υ(2S)','-'),(s,'#29934b','Segnale DCB, massa 9, r = 1','-')]:
  ax.stairs(arr,edges,color=color,lw=1.5,label=label,linestyle=ls)
 ax.set_xlabel('Massa [GeV]');ax.set_ylabel('Eventi / 0.01333 GeV');ax.set_ylim(bottom=0);ax.grid(alpha=.15)
axs[0].set_xlim(6,12);axs[1].set_xlim(8,10.5);axs[1].set_ylim(0,max(n[(x>8)&(x<10.5)])*1.15)
axs[0].legend(fontsize=7.5);fig.suptitle('2022 — Chebyshev + Υ(1S) + Υ(2S) + segnale a 9 GeV\nParametri originali, somma normalizzata ai dati, nessun fit',fontsize=12)
base=p/'xee_2022_chebyshev_upsilon_senza_fit'
fig.savefig(str(base)+'.png',dpi=180);fig.savefig(str(base)+'.pdf')
np.savetxt(str(base)+'.csv',np.c_[x,n,s,b,u1,u2,btotal,y],delimiter=',',header='mass,data,signal,chebyshev,upsilon1s,upsilon2s,total_background,total',comments='')
r=dict(data_total=float(n.sum()),common_scale=float(scale),nominal_yields={k:meta[k] for k in names},scaled_yields={k:float(z.sum()) for k,z in zip(names,raw*scale)},total=float(y.sum()),fit_performed=False,era='2022',r=1,signal_mass_GeV=9)
Path(str(base)+'.json').write_text(json.dumps(r,indent=2));print(json.dumps(r,indent=2));assert abs(y.sum()-n.sum())<1e-6
