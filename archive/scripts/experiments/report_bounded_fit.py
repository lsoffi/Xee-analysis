import json,csv,numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path
p=Path('outputs');bounds=json.load(open(p/'chebyshev_degree5_bounds.json'));tags=['rebin8_degree5','rebin8_degree5_limit1','rebin8_degree5_limit2'];labels=['Coefficienti liberi','Limiti rigidi ±1σ pre-fit','Limiti rigidi ±2σ pre-fit'];summary=[];params=[]
fig,axs=plt.subplots(3,2,figsize=(13,12),layout='constrained',sharex=True,sharey=True)
for i,(tag,title) in enumerate(zip(tags,labels)):
 m=json.load(open(p/f'combine_2022_{tag}_pre_post.json'));r=m['parameters'][0]
 row=dict(case=title,r=r['s'],r_lo=r['s_lo'],r_hi=r['s_hi'],status_b=m['b']['status'],covQual_b=m['b']['covQual'],status_s=m['s']['status'],covQual_s=m['s']['covQual'])
 for a in m['parameters']:
  a=dict(a);a['case']=title
  if a['parameter'].startswith('t'):
   a['shift_b_in_saved_sigma']=(a['b']-a['prefit'])/a['prefit_saved_error'];a['shift_s_in_saved_sigma']=(a['s']-a['prefit'])/a['prefit_saved_error']
  params.append(a)
 for j,(stage,stitle) in enumerate([('fit_b','Solo fondo'),('fit_s','Segnale + fondo')]):
  ax=axs[i,j];a=np.genfromtxt(p/f'combine_2022_{tag}_{stage}_shapes.csv',delimiter=',',names=True);x=a['mass'];n=a['data'];mu=a['total'];ed=np.r_[a['low'],a['high'][-1]];scale=(8/75)/a['width'];row[stage+'_deviance']=float(2*np.sum(mu-n+n*np.log(n/mu)))
  ax.errorbar(x,n*scale,yerr=np.sqrt(n)*scale,fmt='.',color='black',ms=4,lw=.7,label='Dati 2022')
  for key,c,lab,ls in [('total','#c62828','Totale','-'),('dy','#1565c0','Chebyshev grado 5','--'),('upsilon1s','#ef8a00','Υ(1S), norma fissa','-.'),('upsilon2s','#9c27b0','Υ(2S)',':'),('signal','#188341','Segnale 9 GeV','-')]:
   y=a[key]*scale
   if key=='signal' and np.any(y<0):y=np.abs(y);lab='|Segnale|, r negativo';ls='--'
   if np.any(y>0):ax.stairs(np.where(y>0,y,np.nan),ed,color=c,lw=1.5,ls=ls,label=lab,baseline=None)
  ax.set_yscale('log');ax.set_xlim(6,12);ax.set_ylim(.001,10000);ax.grid(alpha=.12,which='both');ax.set_title(title+' — '+stitle,fontsize=10)
  ax.legend(fontsize=7,loc='lower left',ncol=2)
  if stage=='fit_s':ax.text(.04,.94,f'r = {r["s"]:.3f} −{abs(r["s_lo"]):.3f}/+{r["s_hi"]:.3f}',transform=ax.transAxes,va='top',fontsize=9,bbox=dict(facecolor='white',alpha=.9,edgecolor='none'))
  if j==0:ax.set_ylabel('Eventi / 0.10667 GeV')
  if i==2:ax.set_xlabel('Massa [GeV]')
 summary.append(row)
fig.suptitle('2022 — fit 6–12 GeV, rebin ×8, Chebyshev grado 5\nConfronto libero / limiti pre-fit ±1σ / ±2σ; stessi vincoli e opzioni Combine',fontsize=13)
for ext in ['png','pdf']:fig.savefig(p/f'chebyshev_degree5_bounded_fits.{ext}',dpi=170)
with open(p/'chebyshev_degree5_bounded_summary.csv','w') as f:
 wr=csv.DictWriter(f,fieldnames=summary[0].keys());wr.writeheader();wr.writerows(summary)
fields=list(dict.fromkeys(k for a in params for k in a))
with open(p/'chebyshev_degree5_bounded_parameters.csv','w') as f:
 wr=csv.DictWriter(f,fieldnames=fields);wr.writeheader();wr.writerows(params)
lines=['# Fit 6–12 con intervalli rigidi sui coefficienti Chebyshev','','Canale2022, grado5, rebin×8, 57 bin/33341 eventi. Υ(1S) congelata alla normalizzazione originale; r tra−30 e30; strategia0 e robustFit1. Errori dei coefficienti pre-fit: RooRealVar.getError() salvati nel workspace originale. Intervalli cartesiani indipendenti, non vincoli gaussiani e non ellissoide di covarianza. Non costituiscono una propagazione dell’incertezza del fondo.','', '| Caso | r (68% robustFit) | Devianza B | Devianza S+B | covQual B/S+B |','|---|---:|---:|---:|---|']
for a in summary:lines.append(f"| {a['case']} | {a['r']:.5f} −{abs(a['r_lo']):.4f}/+{a['r_hi']:.4f} | {a['fit_b_deviance']:.3f} | {a['fit_s_deviance']:.3f} | {a['covQual_b']}/{a['covQual_s']} |")
lines+=['','Il caso±1σ forza t4 sul bordo superiore in entrambi i fit. Gli errori HESSE di un parametro al bordo non vanno interpretati come intervalli simmetrici validi. Il caso±2σ è sostanzialmente identico al fit libero: tutti i coefficienti liberi erano già entro±2σ.','', '## Intervalli applicati','','| Parametro | Centro pre-fit | σ salvata | Limiti ±1σ | Limiti ±2σ |','|---|---:|---:|---:|---:|']
for a in bounds[:5]:
 b=next(v for v in bounds if v['parameter']==a['parameter'] and v['nsigma']==2);lines.append(f"| `{a['parameter']}` | {a['prefit']:.8g} | {a['saved_sigma']:.7g} | [{a['low']:.8g}, {a['high']:.8g}] | [{b['low']:.8g}, {b['high']:.8g}] |")
for n,title in [(1,'Limiti±1σ'),(2,'Limiti±2σ')]:
 m=json.load(open(p/f'combine_2022_rebin8_degree5_limit{n}_pre_post.json'));lines+=['','## '+title,'','| Parametro | Pre-fit | Post B | Post S+B | Spostamento B / σ | Spostamento S+B / σ |','|---|---:|---:|---:|---:|---:|']
 for a in m['parameters']:
  if not a['parameter'].startswith('t'):continue
  sb=(a['b']-a['prefit'])/a['prefit_saved_error'];ss=(a['s']-a['prefit'])/a['prefit_saved_error'];lines.append(f"| `{a['parameter']}` | {a['prefit']:.7g} | {a['b']:.7g} ± {a['b_err']:.3g} | {a['s']:.7g} ± {a['s_err']:.3g} | {sb:+.4f} | {ss:+.4f} |")
p.joinpath('chebyshev_degree5_bounded_report.md').write_text('\n'.join(lines)+'\n')
for n in [1,2]:
 tag=f'rebin8_degree5_limit{n}';s=Path('outputs/run_fitdiagnostics_2022_rebin8_degree5.sh').read_text().replace('rebin8_degree5','rebin8_degree5_limit'+str(n));out=p/f'run_fitdiagnostics_2022_{tag}.sh';out.write_text(s);out.chmod(0o755)
