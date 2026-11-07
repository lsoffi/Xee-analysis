import numpy as np,json,csv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path
p=Path('outputs');rows=[]
for degree in [7,6,5,4]:
 for stem,upper in [('rebin8',12),('to11_rebin8',11)]:
  tag=stem+('_degree'+str(degree) if degree<7 else '');meta=json.load(open(p/f'combine_2022_{tag}_pre_post.json'));r=meta['parameters'][0]
  row=dict(degree=degree,range=f'6-{upper}',r=r['s'],r_lo=r['s_lo'],r_hi=r['s_hi'],covQual_b=meta['b']['covQual'],covQual_s=meta['s']['covQual'])
  for stage in ['fit_b','fit_s']:
   a=np.genfromtxt(p/f'combine_2022_{tag}_{stage}_shapes.csv',delimiter=',',names=True);n=a['data'];mu=a['total'];row[stage+'_deviance']=float(2*np.sum(mu-n+n*np.log(n/mu)))
  rows.append(row)
with (p/'confronto_gradi_chebyshev_rebin8.csv').open('w') as f:
 wr=csv.DictWriter(f,fieldnames=rows[0].keys());wr.writeheader();wr.writerows(rows)
lines=['# Confronto gradi Chebyshev — 2022, rebin ×8','','Stessi dati, range, componenti, vincoli e opzioni di Combine. Υ(1S) congelata; tutti i gradi ridotti hanno status0/covQual3 in entrambi i fit e range. Grado7/range6–11: covQual2 nel fit S+B.','', '| Grado | Range | Devianza B | Devianza S+B | r (68% robustFit) | covQual B / S+B |','|---|---|---:|---:|---:|---|']
for r in rows:lines.append(f"| {r['degree']} | {r['range']} | {r['fit_b_deviance']:.2f} | {r['fit_s_deviance']:.2f} | {r['r']:.4f} −{abs(r['r_lo']):.4f}/+{r['r_hi']:.4f} | {r['covQual_b']} / {r['covQual_s']} |")
lines+=['','Devianze Poisson calcolate sui bin completi includendo il residuo; non includono le penalità dei vincoli. Non sono significatività del segnale. Il confronto di r non è un criterio per scegliere il grado. Grado5 è un candidato ridotto, non una selezione validata: la scelta finale richiede i controlli di bias/spurious signal e copertura.','','## Parametri del candidato grado5','']
for upper,stem in [(12,'rebin8'),(11,'to11_rebin8')]:
 m=json.load(open(p/f'combine_2022_{stem}_degree5_pre_post.json'));lines += [f'### Range6–{upper}','','| Parametro | Pre | Post B | Post S+B |','|---|---:|---:|---:|']
 for a in m['parameters']:
  def val(stage):
   if a[stage+'_err'] is None:return f"{a[stage]:.6g} fisso"
   if a['parameter']=='r' and stage=='s':return f"{a[stage]:.6g} −{abs(a[stage+'_lo']):.4g}/+{a[stage+'_hi']:.4g}"
   return f"{a[stage]:.6g} ± {a[stage+'_err']:.3g}"
  lines.append(f"| `{a['parameter']}` | {a['prefit']:.6g} | {val('b')} | {val('s')} |")
 lines+=['']
p.joinpath('confronto_gradi_chebyshev_rebin8.md').write_text('\n'.join(lines)+'\n')
fig,axs=plt.subplots(2,2,figsize=(13,7),gridspec_kw={'height_ratios':[3,1]},layout='constrained',sharex='col')
for j,(stem,upper) in enumerate([('rebin8',12),('to11_rebin8',11)]):
 for degree,color in [(7,'#222222'),(6,'#1f77b4'),(5,'#d62728'),(4,'#e69f00')]:
  tag=stem+('_degree'+str(degree) if degree<7 else '');a=np.genfromtxt(p/f'combine_2022_{tag}_fit_b_shapes.csv',delimiter=',',names=True);x=a['mass'];ed=np.r_[a['low'],a['high'][-1]];scale=(8/75)/a['width'];n=a['data'];mu=a['total']
  if degree==7:axs[0,j].errorbar(x,n*scale,yerr=np.sqrt(n)*scale,fmt='.',color='black',ms=5,label='Dati 2022')
  axs[0,j].stairs(mu*scale,ed,color=color,lw=1.7,label=f'Grado {degree}',baseline=None)
  axs[1,j].plot(x,(n-mu)/np.sqrt(mu),'.-',color=color,lw=.6,ms=3,label=f'Grado {degree}')
 axs[0,j].set_yscale('log');axs[0,j].set_title(f'Fit solo fondo, range6–{upper} GeV');axs[0,j].legend(fontsize=8);axs[0,j].set_ylabel('Eventi / 0.10667 GeV');axs[0,j].grid(alpha=.2);axs[1,j].axhline(0,color='gray',lw=.7);axs[1,j].set_ylabel('(Dati−fit)/√fit');axs[1,j].set_xlabel('Massa [GeV]');axs[1,j].set_xlim(6,upper);axs[1,j].grid(alpha=.2)
fig.suptitle('Rebin ×8 — confronto dei gradi 7, 6, 5 e 4 del Chebyshev',fontsize=14)
for ext in ['png','pdf']:fig.savefig(p/f'confronto_gradi_chebyshev_rebin8.{ext}',dpi=180)
