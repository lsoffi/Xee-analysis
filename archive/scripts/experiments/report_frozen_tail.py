import json,csv,numpy as np,re,shutil
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path
p=Path('outputs');base={a['era']:a for a in json.load(open(p/'degree5_native_all_eras_summary.json'))};new=json.load(open(p/'degree5_frozenTail_all_eras_summary.json'));params=json.load(open(p/'degree5_frozenTail_all_eras_parameters.json'));choices=json.load(open(p/'degree5_frozen_tail_choices.json'));decomp=json.load(open(p/'degree5_tail_curvature_decomposition.json'));eras=[('2022','ch1'),('2022EE','ch2'),('2023','ch3'),('2023BPix','ch4')];quality=[]
for case in [e for e,_ in eras]+['combined']:
 for version,prefix in [('free','degree5_native'),('frozen','degree5_frozenTail')]:
  sums={s:0 for s in ['fit_b','fit_s']}
  for e,ch in (eras if case=='combined' else [next(t for t in eras if t[0]==case)]):
   for st in sums:
    a=np.genfromtxt(p/f'{prefix}_{case}_{ch}_{st}_shapes.csv',delimiter=',',names=True);n=a['data'];mu=a['total'];sums[st]+=float(2*np.sum(mu-n+np.where(n>0,n*np.log(np.maximum(n,1)/mu),0)))
  quality.append(dict(era=case,model=version,**sums))
with open(p/'degree5_frozenTail_fit_quality.csv','w') as f:
 wr=csv.DictWriter(f,fieldnames=quality[0].keys());wr.writeheader();wr.writerows(quality)
fig,axs=plt.subplots(2,4,figsize=(17,7),layout='constrained')
for j,(era,ch) in enumerate(eras):
 for i,(lo,hi) in enumerate([(6,12),(10.8,12)]):
  ax=axs[i,j]
  for prefix,c,ls,label in [('degree5_native','#1766ba','--','Post S+B libero'),('degree5_frozenTail','#cc2d2d','-','Post S+B congelato')]:
   a=np.genfromtxt(p/f'{prefix}_combined_{ch}_fit_s_shapes.csv',delimiter=',',names=True);ax.plot(a['mass'],a['dy'],color=c,ls=ls,label=label,lw=1.7)
  a=np.genfromtxt(p/f'degree5_frozenTail_combined_{ch}_prefit_shapes.csv',delimiter=',',names=True);ax.plot(a['mass'],a['dy'],color='#444444',ls=':',label='Pre-fit',lw=1.3);ax.set_yscale('log');ax.set_xlim(lo,hi);visible=np.concatenate([line.get_ydata()[(line.get_xdata()>=lo)&(line.get_xdata()<=hi)] for line in ax.lines]);ax.set_ylim(visible.min()*.8,visible.max()*1.2);ax.set_title(era+(' — range completo' if i==0 else ' — coda'));ax.set_xlabel('Massa [GeV]');ax.set_ylabel('Eventi / bin');ax.grid(alpha=.15,which='both');ax.legend(fontsize=7)
fig.suptitle('Chebyshev, fit combinato — prima e dopo il congelamento selettivo\nt4 congelato nel2022; t3 nellealtreere; altri coefficienti profilati',fontsize=13)
for ext in ['png','pdf']:fig.savefig(p/f'degree5_frozenTail_chebyshev_comparison.{ext}',dpi=170)
lines=['# Report — congelamento selettivo della curvatura del Chebyshev','','## Configurazione','','Input Xee_ee_4_allYears.root, segnale DCB a massa9GeV, tutte le ere. Range6–12GeV,450 bin per era, Chebyshev grado5. Fit individuali e combinato. r tra−30 e30, strategia0, robustFit1; saveNormalizations e saveShapes. Restano congelate le normalizzazioni Υ(1S); normalizzazioni Chebyshev e Υ(2S), altri coefficienti e nuisance restano profilati.','', '## Scelta dei coefficienti','','Diagnostica sul fit combinato precedente: variazione pre→post del contributo dei singoli termini alla curvatura discreta P(11.5)−[P(11)+P(12)]/2 del polinomio non normalizzato. Questa decomposizione identifica il contributo dominante, non una causa fisica indipendente né un parametro che controlla soltanto la coda. Per evitare valutazioni non fisiche durante il test, si usa la decomposizione algebrica del polinomio: riportare isolatamente alcuni coefficienti ai valori iniziali può rendere il PDF negativo.','', '| Era | Coefficiente congelato | Valore | Coefficienti ancora liberi |','|---|---|---:|---|']
for a in choices:
 if a['fit']!='combined':continue
 fixed=a['parameter'].split('_')[0];others=', '.join(f't{i}' for i in range(5) if f't{i}'!=fixed);lines.append(f"| {a['era']} | `{a['parameter']}` | {a['frozen_value']:.9g} | {others} |")
lines+=['','**Riferimento dei valori congelati:** sono i valori iniziali del nuovo modello di grado5, che usava i primi5 coefficienti2022 come seed per tutte le ere. Non sono una misura indipendente della shape del fondo né vincoli esterni.','', '## Risultati','','| Era | r prima | r dopo (68% robustFit) | Status B/S+B | covQual B/S+B |','|---|---:|---:|---|---|']
for a in new:lines.append(f"| {a['era']} | {base[a['era']]['r']:.5f} | {a['r']:.5f} −{abs(a['r_lo']):.5f}/+{a['r_hi']:.5f} | {a['status_b']}/{a['status_s']} | {a['covQual_b']}/{a['covQual_s']} |")
lines+=['','## Qualità della descrizione dei dati','','Devianze Poisson dai bin, senza penalità dei vincoli. Non sono valori di significanza.','', '| Era | B libero | B congelato | S+B libero | S+B congelato |','|---|---:|---:|---:|---:|']
for case in [e for e,_ in eras]+['combined']:
 old=next(q for q in quality if q['era']==case and q['model']=='free');now=next(q for q in quality if q['era']==case and q['model']=='frozen');lines.append(f"| {case} | {old['fit_b']:.3f} | {now['fit_b']:.3f} | {old['fit_s']:.3f} | {now['fit_s']:.3f} |")
lines+=['','## Significanza combinata','','Significance asintotica locale ottenuta con lo stesso congelamento: **Z=3.75619, p0=8.62595e−5**. Con minimizzazione di controllo strategia1 e tolleranza0.001: **Z=3.75616**. Il rapporto delle likelihood di FitDiagnostics dà −2ΔNLL=14.10887, coerente con Z. Prima del congelamento il risultato era Z=0.','', '**Interpretazione:** è la significanza condizionata a questo congelamento empirico. Non va presentata come significanza validata dell’analisi: i coefficienti sono stati scelti osservando gli stessi dati e i valori fissati sono seed numerici senza propagazione dell’incertezza. Il forte aumento dimostra la sensibilità della significanza alla rigidità imposta sul fondo. Prima di un risultato fisico servono validazione del bias/spurious signal e della procedura completa. Il congelamento non garantisce l’eliminazione della gobba: gli altri coefficienti possono compensare.','', '## Parametri pre/post','','Valori HESSE per i parametri liberi, intervallo robustFit per r; i parametri congelati sono elencati sopra.']
for case in [e for e,_ in eras]+['combined']:
 lines+=['','### '+case,'','| Parametro | Pre-fit | Post B | Post S+B |','|---|---:|---:|---:|'];group=[a for a in params if a['era']==case];names=list(dict.fromkeys(a['parameter'] for a in group))
 for name in names:
  b=next((a for a in group if a['fit']=='fit_b' and a['parameter']==name),None);s=next((a for a in group if a['fit']=='fit_s' and a['parameter']==name),None);a=s or b
  def fmt(v):
   if v is None:return '—'
   if v['parameter']=='r':return f"{v['postfit']:.6g} −{abs(v['error_lo']):.4g}/+{v['error_hi']:.4g}"
   return f"{v['postfit']:.7g} ± {v['error']:.3g}"
  lines.append(f"| `{name}` | {a['prefit']:.7g} | {fmt(b)} | {fmt(s)} |")
p.joinpath('degree5_frozenTail_fit_report.md').write_text('\n'.join(lines)+'\n')
print('\n'.join(lines[:45]))
