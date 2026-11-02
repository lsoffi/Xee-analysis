import pathlib,json,csv,hashlib,math
import numpy as np
import uproot,re
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import norm
P=pathlib.Path('/Users/imb/Documents/Codex/Analisi-statistica-Xee');R=P/'runs/20261029_local_gaussian_scan_7_8';old=json.loads((P/'runs/20261029_gaussian_signal_scan/summary/selected.json').read_text());rows=[];flags=[]
for i in range(70,81):
 m=i/10;d=R/f'm{m:.1f}';r=json.loads((d/'results.json').read_text());a=json.loads((d/'audit.json').read_text());q=json.loads((d/'quality.json').read_text());assert hashlib.sha256(pathlib.Path(a['source']).read_bytes()).hexdigest()==a['source_sha256'];assert r['limit_crossing']['verified']
 for st in ['b','s']:
  if q[st]['status']!=0 or q[st]['covQual']!=3:flags.append(dict(mass=m,stage=st,status=q[st]['status'],covQual=q[st]['covQual']))
  for v in q[st]['parameters']:
   if v['name'].startswith('local_') or v['name']=='r':
    if min(v['value']-v['minimum'],v['maximum']-v['value'])<1e-4*(v['maximum']-v['minimum']):flags.append(dict(mass=m,stage=st,bound=v['name']))
 for role in ['fitdiagnostics','significance','limits']:
  assert json.loads((d/role/'receipt.json').read_text())['exit_code']==0
  log=(d/role/'combine.log').read_text()
  if 'ERROR:' in log:flags.append(dict(mass=m,role=role,error_lines=[x for x in log.splitlines() if 'ERROR:' in x]))
 sig=next(x for x in old if x['model']=='gaussian' and x['mass']==m and x['role']=='significance');lim=next(x for x in old if x['model']=='gaussian' and x['mass']==m and x['role']=='limits');assert lim['valid_observed']
 original_limit=lim['observed']
 if m==7.4:
  check=R/'global_check_m7.4'
  with uproot.open(next(check.glob('higgsCombine*.root'))) as f:
   vs=f['limit']['limit'].array(library='np');qs=f['limit']['quantileExpected'].array(library='np')
   lim=dict(lim)
   for key,quantile in [('observed',-1),('expected_50',.5),('expected_16',.16),('expected_84',.84)]:lim[key]=float(vs[np.argmin(abs(qs-quantile))])
  matches=re.findall(r'At r =\s*([\d.eE+-]+):[^\n]*CLs\s*=\s*([\d.eE+-]+)',(check/'combine.log').read_text());near=min(matches,key=lambda a:abs(float(a[0])-lim['observed']));assert abs(float(near[1])-.05)<.01 and abs(float(near[0])/lim['observed']-1)<.02
 if r['r']['r']+r['r']['rHiErr']>=4.999 or r['r']['r']-r['r']['rLoErr']<=-4.999:flags.append(dict(mass=m,interval_bound=True))
 row=dict(global_limit_original_scan=original_limit,mass=m,local_Z=r['Z'],global_Z=sig['Z'],local_p0=r['p0'],global_p0=sig['p0'],local_r=r['r']['r'],local_limit=r['limits']['observed'],global_limit=lim['observed'],local_expected=r['limits']['expected_50'],global_expected=lim['expected_50'],local_expected16=r['limits']['expected_16'],local_expected84=r['limits']['expected_84'],global_expected16=lim['expected_16'],global_expected84=lim['expected_84'],local_acceptance_min=min(x['signal_fraction'] for x in a['eras']),local_acceptance_max=max(x['signal_fraction'] for x in a['eras']),local_B_deviance=r['B_deviance'],local_B_gof_p=r['B_approx_p']);rows.append(row)
with (R/'comparison.csv').open('w') as f:w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
(R/'comparison.json').write_text(json.dumps(rows,indent=2));(R/'quality_flags.json').write_text(json.dumps(flags,indent=2));x=np.array([r['mass'] for r in rows]);val=lambda k:np.array([r[k] for r in rows]);colors={'local':'#d66b17','global':'#256ca5'};labels={'local':'Esponenziale + Gaussiana, fit 7–8','global':'Chebyshev7 + Gaussiana, fit 6–12'}
for role in ['limits','pvalues']:
 fig,ax=plt.subplots(figsize=(10.5,6.2))
 for model in ['global','local']:
  if role=='limits':
   ax.plot(x,val(model+'_limit'),'o-',color=colors[model],label=labels[model]+' — osservato',lw=2,ms=4)
   ax.plot(x,val(model+'_expected'),'--',color=colors[model],label=labels[model]+' — atteso mediano',lw=1.5)
  else:ax.plot(x,val(model+'_p0'),'o-',color=colors[model],label=labels[model],lw=2,ms=4)
 ax.set_yscale('log');ax.set_xlim(7,8);ax.set_xticks(x);ax.grid(alpha=.2,which='both');ax.set_xlabel('Massa ipotizzata del segnale [GeV]');ax.legend(fontsize=9)
 if role=='limits':ax.set_ylabel('Limite superiore su r, 95% CLs');ax.set_title('Limiti: confronto scan locale e fit globale\nQuattro ere; stessa Gaussiana MC e convenzione di r')
 else:
  ax.set_ylabel('p-value locale p₀');ax.set_ylim(min(.001,val('local_p0').min()/2,val('global_p0').min()/2),1);ax.set_title('p-value: confronto scan locale e fit globale\nFinestra locale fissa 7–8 GeV per tutte le ipotesi di massa')
  for z in [1,2,3]:
   p=norm.sf(z)
   if p>ax.get_ylim()[0]:ax.axhline(p,color='gray',ls=':',lw=.8);ax.text(7.985,p,f'{z}σ',ha='right',va='bottom',color='gray')
 
 if role=='limits':
  fig.text(.5,.012,'Punto Cheb7 a 7,4 GeV ricontrollato: strategia 1, tolleranza 0,001 (dettagli nel report).',ha='center',fontsize=8,color='#555555');fig.tight_layout(rect=[0,.035,1,1])
 else:fig.tight_layout()
 for ext in ['png','pdf']:fig.savefig(R/f'{role}_comparison.{ext}',dpi=160)
 plt.close(fig)
fig,ax=plt.subplots(figsize=(10,5.5))
for model in ['global','local']:
 ax.fill_between(x,val(model+'_expected16'),val(model+'_expected84'),color=colors[model],alpha=.15)
 ax.plot(x,val(model+'_expected'),'o--',color=colors[model],label=labels[model]+'; banda ±1σ',ms=4)
ax.set_yscale('log');ax.set_xlim(7,8);ax.set_xticks(x);ax.set_xlabel('Massa ipotizzata [GeV]');ax.set_ylabel('Limite atteso su r, 95% CLs');ax.set_title('Sensibilità attesa, mediana e bande 68%');ax.legend(fontsize=9);ax.grid(alpha=.2,which='both');fig.tight_layout()
for ext in ['png','pdf']:fig.savefig(R/f'expected_comparison.{ext}',dpi=160)
plt.close(fig)
tr=''.join(f"<tr><td>{r['mass']:.1f}</td><td>{r['local_limit']:.5f}</td><td>{r['global_limit']:.5f}</td><td>{r['local_p0']:.5g}</td><td>{r['global_p0']:.5g}</td><td>{r['local_Z']:.3f}</td><td>{r['global_Z']:.3f}</td><td><a href='m{r['mass']:.1f}/fits.png'>Fit</a> · <a href='m{r['mass']:.1f}/parameters.csv'>Parametri</a></td></tr>" for r in rows)
(R/'index.html').write_text(f'''<!doctype html><meta charset="utf-8"><title>Scan locale 7–8 contro Cheb7 globale</title><style>body{{font:17px system-ui;max-width:1100px;margin:40px auto;padding:20px;color:#263647}}img{{width:100%}}table{{border-collapse:collapse}}th,td{{padding:9px;border:1px solid #ddd}}</style><h1>Scan locale 7–8 GeV e confronto Chebyshev7 + Gaussiana</h1><p>11 masse da 7,0 a 8,0 GeV, passo 0,1. Per ogni massa: FitDiagnostics, AsymptoticLimits e Significance reali CMS Combine, quattro ere, binning originale (75 bin/era nel locale), r comune. Range locale fisso 7–8 per tutti i punti; il confronto blu proviene dallo scan globale 6–12 già eseguito. Non è un fit Chebyshev nella finestra locale.</p><p>Fondo locale: una esponenziale per era, pendenza e normalizzazione floating, descrizione del fondo totale comprese code Υ. Segnale: Gaussiana con stessi nodi di media/σ e nuisance MC, yield per r identico al test globale prima della selezione. Accettanza nella finestra calcolata analiticamente e profilata con le nuisance, senza rinormalizzare artificialmente il segnale completo nella finestra. Ai bordi 7 e 8 GeV si osserva solo circa metà Gaussiana, quindi la sensibilità e l'identificabilità cambiano.</p><p>33 esecuzioni principali concluse, più4controlli intervalliFD e1ricontrollo del limite globale; limiti osservati verificati all'attraversamento CLs≈0,05. Flag di qualità: {len(flags)} (vedi <a href='quality_flags.json'>dettagli</a>). Input immutati, checksum verificati su tutti gli11workspace. Fit strategia1/tolleranza0,001, robustFit1 perFD; range r[-5,5] perFD/significanza, rMax1 perlimiti. Nessun vincolo del fondo derivato dagli stessi dati.</p><p><b>Correzione numerica a7,4GeV nella curva globale:</b> il valore storico0,79300 è stato ricontrollato con strategia1/tolleranza0,001 e rMax1: limite0,07046 (CLs≈0,04992). La curva blu usa il punto ricalcolato; il valore storico rimane nella colonna global_limit_original_scan del CSV e negli output originali immutati. Questa variazione evidenzia sensibilità numerica: il solo attraversamentoCLs non garantisce un profilo corretto. Il p-value globale rimane quello dello scan precedente.</p><p>I4fitFD locali7,0–7,3 sono ricalcolati con stepSize0,001 e strategia1 anche nella ricerca robusta degli intervalli, dopo che il passo iniziale attraversava regioni con previsioni negative; i tentativi iniziali sono archiviati. Limiti e significanze locali non sono modificati da questo controllo.</p><h2>Limiti osservati e attesi</h2><a href='limits_comparison.pdf'>PDF</a><img src='limits_comparison.png'><h2>p-value locale</h2><a href='pvalues_comparison.pdf'>PDF</a><img src='pvalues_comparison.png'><h2>Sensibilità attesa</h2><a href='expected_comparison.pdf'>PDF</a><img src='expected_comparison.png'><h2>Tabella e fit</h2><p><a href='comparison.csv'>Scarica CSV</a></p><table><tr><th>Massa</th><th>Limite exp locale</th><th>Limite Cheb globale</th><th>p₀ exp</th><th>p₀ Cheb</th><th>Z exp</th><th>Z Cheb</th><th>Output locali</th></tr>{tr}</table><p>I p-value sono locali, non corretti per lo scan. Cambiano insieme modello, range e trattamento delle code Υ: il confronto non attribuisce le differenze a una sola causa. Convergenza e bontà dei residui non sostituiscono test bias/injection/coverage. Il modello non viene scelto minimizzando la significanza osservata.</p>''')
(R/'completion.json').write_text(json.dumps(dict(complete=True,masses=len(rows),new_combine_jobs=33,verified_limits=11,quality_flags=len(flags),source_checksums_unchanged=True),indent=2))
print(json.dumps(rows,indent=2));print('FLAGS',flags)
