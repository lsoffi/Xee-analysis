from pathlib import Path
import json,csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
R=Path('runs/20261029_upsilon_mu_sigma');A=json.loads((R/'results.json').read_text());eras=['2022','2022EE','2023','2023BPix'];plt.rcParams.update({'font.size':10})
fig,axes=plt.subplots(2,4,figsize=(15,7),sharex=True)
for j,era in enumerate(eras):
 for row,param in enumerate(['upsilon_mu','upsilon_sigma']):
  ax=axes[row,j]
  for deg,color in [(0,'tab:orange'),(1,'tab:blue')]:
   aa=[a for a in A if a['era']==era and a['background_degree']==deg]
   for a in aa:
    p=a['parameters'][param];xx=a['multiplier']+(-.025 if deg==0 else .025);ax.plot(xx,p['value'],marker='o' if a['reliable_errors'] else 'x',color=color)
    if a['reliable_errors']:ax.errorbar(xx,p['value'],yerr=[[-p['lo']],[p['hi']]],fmt='none',color=color,capsize=3)
   ax.plot([],[],color=color,label='Fondo costante' if deg==0 else 'Fondo lineare')
  ax.axhline(aa[0]['nominal_mu' if row==0 else 'nominal_sigma'],ls='--',color='gray',label='Nominale MC');ax.grid(alpha=.2);ax.set_title(era if row==0 else '');ax.set_ylabel(('mu' if row==0 else 'sigma core')+' [GeV]');ax.set_xticks([1.5,2,3]);ax.set_xlabel('Semifinestra / sigma MC' if row else '')
axes[0,0].legend(fontsize=8);fig.suptitle('Upsilon(1S): code DCB fisse; mu e sigma liberi\nIntervalli profilati 68%; croce = errore/qualita numerica non affidabile',fontsize=13);fig.tight_layout();fig.savefig(R/'parameters.png',dpi=150);fig.savefig(R/'parameters.pdf');plt.close(fig)
for k in [1.5,2,3]:
 for deg in [0,1]:
  fig,axes=plt.subplots(2,2,figsize=(12,8))
  for ax,era in zip(axes.flat,eras):
   a=next(a for a in A if a['era']==era and a['multiplier']==k and a['background_degree']==deg);xx=[p['mass'] for p in a['data']];yy=np.array([p['events'] for p in a['data']]);ax.errorbar(xx,yy,yerr=np.sqrt(yy),fmt='.',color='black',label='Dati',ms=4)
   for key,label,color in [('total','Totale','tab:red'),('Y1','Upsilon(1S)','tab:blue'),('Y2','Upsilon(2S)','tab:purple'),('continuum','Continuo','tab:green')]:ax.plot(xx,a['curves'][key],label=label,color=color)
   ax.set_title(f"{era}: status={a['status']}, covQual={a['covQual']}"+(' [intervalli problematici]' if not a['reliable_errors'] else ''));ax.set_xlabel('m [GeV]');ax.set_ylabel('Eventi / 13.33 MeV');ax.set_ylim(bottom=0);ax.legend(fontsize=8);ax.grid(alpha=.2)
  fig.suptitle(f"Finestra +/-{k:g} sigma MC: {a['range'][0]:.3f}–{a['range'][1]:.3f} GeV; fondo "+('costante' if deg==0 else 'lineare'));fig.tight_layout();fig.savefig(R/f'fits_k{k:g}_b{deg}.png',dpi=150);fig.savefig(R/f'fits_k{k:g}_b{deg}.pdf');plt.close(fig)
fields=['tag','era','multiplier','background_degree','status','covQual','profile_failed','boundaries','mu','mu_lo','mu_hi','sigma','sigma_lo','sigma_hi','Y2_scale']
with (R/'parameters.csv').open('w') as f:
 writer=csv.DictWriter(f,fieldnames=fields);writer.writeheader()
 for a in A:
  row={k:a[k] for k in fields if k in a};row['boundaries']=';'.join(a['boundaries'])
  for n,pn in [('mu','upsilon_mu'),('sigma','upsilon_sigma')]:
   p=a['parameters'][pn];row.update({n:p['value'],n+'_lo':p['lo'] if a['reliable_errors'] else '',n+'_hi':p['hi'] if a['reliable_errors'] else ''})
  row['Y2_scale']=a['parameters']['scale_upsilon2s']['value'];writer.writerow(row)
html='<html><meta charset="utf-8"><title>Misura Upsilon mu e sigma</title><body><h1>Misura locale Upsilon(1S)</h1><p>24 fit Combine. Code MC fisse; mu e sigma liberi. Normalizzazioni libere, scala Y2 nel dominio originale [0,10]. Tre risultati con problemi numerici negli intervalli; nessuna finestra minima validata da questo solo studio.</p><a href="REPORT.md">Report</a> | <a href="parameters.csv">Parametri CSV</a><br><img width="1100" src="parameters.png">'
for k in [1.5,2,3]:
 for deg in [0,1]:html+=f'<h2>±{k} sigma, fondo grado {deg}</h2><img width="1100" src="fits_k{k:g}_b{deg}.png">'
(R/'index.html').write_text(html+'</body></html>')
