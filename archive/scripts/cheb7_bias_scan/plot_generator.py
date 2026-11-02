from pathlib import Path
import json,numpy as np,uproot
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
R=Path('runs/20261030_cheb7_bias_scan');f=uproot.open(R/'generator_fit/fitDiagnostics_generator.root');f3=uproot.open('runs/20261030_fixed_sr_scan/m7.5/joint/fitdiagnostics/fitDiagnostics_local.root')
for region in ['SR','CR']:
 fig=plt.figure(figsize=(13,10));outer=fig.add_gridspec(2,2,hspace=.32,wspace=.2)
 for i,(era,cell) in enumerate(zip(['2022','2022EE','2023','2023BPix'],outer),1):
  d=f[f'shapes_fit_b/{region}_ch{i}'];v,e=d['total'].to_numpy();bw=np.diff(e);x=(e[:-1]+e[1:])/2;data=d['data'].values('y')*bw;pred=v*bw;gs=cell.subgridspec(2,1,height_ratios=[3,1],hspace=.06);ax=fig.add_subplot(gs[0]);rx=fig.add_subplot(gs[1],sharex=ax);ax.errorbar(x,data,yerr=np.sqrt(data),fmt='k.',ms=2,label='Dati');ax.stairs(pred,e,color='tab:red',label='Totale generatore grado7');ax.stairs(f3[f'shapes_fit_b/{region}_ch{i}/total'].values()*bw,e,color='gray',ls='--',label='Totale grado3')
  for proc,color in [('dy','tab:blue'),('upsilon1s','orange'),('upsilon2s','purple')]:ax.stairs(d[proc].values()*bw,e,color=color,label=proc)
  ax.set_title(era);ax.set_ylabel('Eventi / bin');ax.tick_params(labelbottom=False);ax.legend(fontsize=7);ax.grid(alpha=.2);rx.axhline(0,color='gray');rx.plot(x,(data-pred)/np.sqrt(pred),'.',ms=2,color='tab:red');rx.set_xlabel('Massa [GeV]');rx.set_ylabel('(D-F)/sqrt(F)',fontsize=8);rx.grid(alpha=.2)
 fig.suptitle(f'Generatore: Chebyshev7 + Upsilon1S + Upsilon2S, solo fondo — {region}\nFit ai dati per definire uno scenario alternativo, non una verita fisica');fig.subplots_adjust(top=.90,bottom=.07);fig.savefig(R/f'generator_{region}.png',dpi=140);fig.savefig(R/f'generator_{region}.pdf');plt.close(fig)
