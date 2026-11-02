import pathlib,json,csv,hashlib
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import norm
P=pathlib.Path('/Users/imb/Documents/Codex/Analisi-statistica-Xee');R=P/'runs/20261029_local_dcb_scan_7_8';prior=P/'runs/20261029_local_gaussian_scan_7_8';rows=json.loads((prior/'comparison.json').read_text());flags=[]
for row in rows:
 m=row['mass'];d=R/f'm{m:.1f}';r=json.loads((d/'results.json').read_text());q=json.loads((d/'quality.json').read_text());a=json.loads((d/'audit.json').read_text());assert hashlib.sha256(pathlib.Path(a['source']).read_bytes()).hexdigest()==a['source_sha256'];assert r['limit_crossing']['verified']
 for st in ['b','s']:
  if q[st]['status']!=0 or q[st]['covQual']!=3:flags.append(dict(mass=m,stage=st,status=q[st]['status'],covQual=q[st]['covQual']))
  for v in q[st]['parameters']:
   if v['name'].startswith('local_') or v['name']=='r':
    if min(v['value']-v['minimum'],v['maximum']-v['value'])<1e-4*(v['maximum']-v['minimum']):flags.append(dict(mass=m,stage=st,bound=v['name']))
 for role in ['fitdiagnostics','significance','limits']:
  assert json.loads((d/role/'receipt.json').read_text())['exit_code']==0
  log=(d/role/'combine.log').read_text()
  if 'ERROR:' in log:flags.append(dict(mass=m,role=role,error_lines=[x for x in log.splitlines() if 'ERROR:' in x]))
 if r['r']['r']+r['r']['rHiErr']>=4.999 or r['r']['r']-r['r']['rLoErr']<=-4.999:flags.append(dict(mass=m,interval_bound=True))
 row.update(dcb_Z=r['Z'],dcb_p0=r['p0'],dcb_r=r['r']['r'],dcb_limit=r['limits']['observed'],dcb_expected=r['limits']['expected_50'],dcb_expected16=r['limits']['expected_16'],dcb_expected84=r['limits']['expected_84'],dcb_acceptance_min=min(e['signal_fraction'] for e in a['eras']),dcb_acceptance_max=max(e['signal_fraction'] for e in a['eras']))
with (R/'comparison.csv').open('w') as f:w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
(R/'comparison.json').write_text(json.dumps(rows,indent=2));(R/'quality_flags.json').write_text(json.dumps(flags,indent=2));x=np.array([r['mass'] for r in rows]);val=lambda k:np.array([r[k] for r in rows]);colors={'local':'#d66b17','global':'#256ca5','dcb':'#188548'};labels={'local':'Exp + Gaussiana, fit 7–8','global':'Chebyshev7 + Gaussiana, fit 6–12','dcb':'Exp + DCB workspace, fit 7–8'}
for role in ['limits','pvalues','expected']:
 fig,ax=plt.subplots(figsize=(11,6.5))
 for model in ['global','local','dcb']:
  if role=='limits':
   ax.plot(x,val(model+'_limit'),'o-',color=colors[model],label=labels[model]+' — osservato',lw=2,ms=4)
   ax.plot(x,val(model+'_expected'),'--',color=colors[model],label=labels[model]+' — atteso',lw=1.5)
  elif role=='expected':
   ax.plot(x,val(model+'_expected'),'o--',color=colors[model],label=labels[model],lw=1.8,ms=4);ax.fill_between(x,val(model+'_expected16'),val(model+'_expected84'),color=colors[model],alpha=.12)
  else:ax.plot(x,val(model+'_p0'),'o-',color=colors[model],label=labels[model],lw=2,ms=4)
 ax.set_yscale('log');ax.set_xlim(7,8);ax.set_xticks(x);ax.grid(alpha=.2,which='both');ax.set_xlabel('Massa ipotizzata del segnale [GeV]');ax.legend(fontsize=8.5)
 if role=='limits':ax.set_ylabel('Limite superiore su r, 95% CLs');ax.set_title('Limiti: aggiunta del segnale DCB originale del workspace\nQuattro ere combinate; stessa convenzione di r')
 elif role=='expected':ax.set_ylabel('Limite atteso su r, 95% CLs');ax.set_title('Sensibilità attesa: mediana e bande 68%')
 else:
  ax.set_ylabel('p-value locale p₀');ax.set_ylim(min(.001,*[val(k+'_p0').min()/2 for k in labels]),1);ax.set_title('p-value: aggiunta del segnale DCB originale del workspace\nFinestra locale fissa 7–8 GeV per tutte le ipotesi di massa')
  for z in [1,2,3,4]:
   p=norm.sf(z)
   if p>ax.get_ylim()[0]:ax.axhline(p,color='gray',ls=':',lw=.8);ax.text(7.985,p,f'{z}σ',ha='right',va='bottom',color='gray')
 if role in ['limits','expected']:fig.text(.5,.012,'Punto Cheb7 a 7,4 GeV: controllo numerico già documentato nel confronto precedente.',ha='center',fontsize=8,color='#555555');fig.tight_layout(rect=[0,.035,1,1])
 else:fig.tight_layout()
 for ext in ['png','pdf']:fig.savefig(R/f'{role}_comparison.{ext}',dpi=155)
 plt.close(fig)
tr=''.join(f"<tr><td>{r['mass']:.1f}</td><td>{r['dcb_limit']:.5f}</td><td>{r['dcb_expected']:.5f}</td><td>{r['dcb_p0']:.5g}</td><td>{r['dcb_Z']:.3f}</td><td>{r['dcb_acceptance_min']:.3f}–{r['dcb_acceptance_max']:.3f}</td><td><a href='m{r['mass']:.1f}/fits.png'>Fit 4 ere</a> · <a href='m{r['mass']:.1f}/fit_sum.png'>Somma/componenti</a> · <a href='m{r['mass']:.1f}/parameters.csv'>Parametri</a></td></tr>" for r in rows)
gallery=''.join(f"<h3>Massa {r['mass']:.1f} GeV</h3><p><a href='m{r['mass']:.1f}/fits.pdf'>PDF</a> · <a href='m{r['mass']:.1f}/fit_sum.pdf'>Somma e segnale PDF</a></p><img loading='lazy' src='m{r['mass']:.1f}/fits.png'>" for r in rows)
(R/'index.html').write_text(f'''<!doctype html><meta charset="utf-8"><title>Esponenziale + DCB, scan 7–8</title><style>body{{font:17px system-ui;max-width:1150px;margin:35px auto;padding:20px;color:#263647}}img{{width:100%}}table{{border-collapse:collapse}}th,td{{padding:8px;border:1px solid #ddd}}a{{color:#176da5}}</style><h1>Esponenziale + signal shape DCB originale</h1><p>11 masse7,0–8,0, passo0,1. Fit locale fisso7–8GeV, quattro ere combinate, r comune,75bin/era originali. Esponenziale totale locale con pendenza e normalizzazione indipendenti perera; comprende le code Υ senza alterare le shape sorgenti. DCB originale RooCrystalBall del workspace di ciascuna massa/era: media, sigma, alphaL/R,nL/R e nuisance originali preservati. Nessun nuovo parametro di forma del segnale libero.</p><p>Yield originale perr1 moltiplicato per l'accettanza DCB7–8/6–12, dipendente dalle nuisance originali. Integrale analitico delle code e nucleo verificato contro ROOT in220configurazioni (44mass/era ×5configurazioni nominale/scale±1/sigma±1); input immutati con checksum verificati. La Gaussiana del confronto non ha le code, quindi la frazione di segnale mantenuta nella finestra è diversa: r mantiene la stessa convenzione di yield totale6–12.</p><p>33esecuzioni Combine: FitDiagnostics robustFit1, AsymptoticLimits95%CLs, Significance locale. Strategia1,tolleranza0,001; intervalli FD stepSize0,001 e strategia1; r[-5,5] inFD/significanza, rMax1 inlimiti. Tutti11limiti verificati con CLs≈0,05. Flag finali: {len(flags)} (<a href='quality_flags.json'>dettagli</a>). I fit B-only sono lo stesso modello locale e gli stessi dati del test Gaussiano; pGOF asintotico≈0,431. Non sono stati eseguiti nuovi test bias/injection/coverage.</p><h2>Curve sovrapposte ai grafici precedenti</h2><p>Verde: nuova exp+DCB7–8. Arancione: exp+Gauss7–8. Blu: Cheb7+Gauss6–12. Curva blu con la correzione numerica7,4 giàdocumentata nel report precedente; nessuna modifica agli output storici.</p><a href='pvalues_comparison.pdf'>PDF p-value</a><img src='pvalues_comparison.png'><a href='limits_comparison.pdf'>PDF limiti</a><img src='limits_comparison.png'><a href='expected_comparison.pdf'>PDF sensibilità</a><img src='expected_comparison.png'><p>Confrontare verde/arancione isola il cambio della shape di segnale a fondo/range uguali. Confrontare con blu cambia anche modello e range. I p-value sono locali, senza correzione per lo scan. Non scegliere il modello minimizzando l'eccesso osservato.</p><h2>Risultati DCB e tutti i fit</h2><a href='comparison.csv'>CSV confronto completo</a><table><tr><th>Massa</th><th>Limite obs.</th><th>Limite atteso</th><th>p₀</th><th>Z</th><th>Accettanza fra ere</th><th>File</th></tr>{tr}</table><p>Nei fit: dati neri, fondo solo blu tratteggiato, S+B rosso e fondo S+B giallo; residui sotto. La normalizzazione del segnale può essere negativa nel fit diagnostico: non corrisponde a un segnale fisico negativo.</p>{gallery}''')
(R/'completion.json').write_text(json.dumps(dict(complete=True,masses=11,combine_jobs=33,verified_observed_limits=11,quality_flags=len(flags),source_unchanged=True,acceptance_checks=220),indent=2));print('FLAGS',flags)
for r in rows:print(r['mass'],r['dcb_Z'],r['dcb_p0'],r['dcb_limit'],r['dcb_expected'])
