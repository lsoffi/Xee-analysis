from pathlib import Path
import numpy as np,uproot
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
P=Path('/Users/imb/Documents/Codex/Analisi-statistica-Xee');R=P/'runs/20261029_upsilon_control_cheb3'
fig,axs=plt.subplots(2,2,figsize=(12,8),sharex=True)
for ax,i,era in zip(axs.flat,range(1,5),['2022','2022EE','2023','2023BPix']):
 for degree,folder,color,style in [(2,'20261029_upsilon_control','#0072B2','-'),(3,'20261029_upsilon_control_cheb3','#D55E00','--')]:
  with uproot.open(P/'runs'/folder/'joint/fitdiagnostics/fitDiagnostics_local.root') as f:
   y,e=f[f'shapes_fit_s/SR_ch{i}/total'].to_numpy();x=(e[:-1]+e[1:])/2;y=y*np.diff(e)
   ax.plot(x,y,color=color,ls=style,lw=2,label=f'Chebyshev grado {degree}')
 ax.set_title(era);ax.set_ylabel('Eventi / 13,33 MeV');ax.set_xlabel('Massa [GeV]');ax.set_xlim(e[0],e[-1]);ax.grid(alpha=.15);ax.legend()
fig.suptitle('Curve totali post-fit S+B — ipotesi di massa 7,5 GeV\nFit simultaneo con regione di controllo Υ',fontsize=14)
fig.tight_layout(rect=[0,0,1,.93])
for ext in ['png','pdf']:fig.savefig(R/f'overlay_total_degree2_degree3.{ext}',dpi=170)
