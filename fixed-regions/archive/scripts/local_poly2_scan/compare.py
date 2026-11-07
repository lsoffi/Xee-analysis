import pathlib,json,csv,hashlib
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import norm
P=pathlib.Path('/Users/imb/Documents/Codex/Analisi-statistica-Xee');R=P/'runs/20261029_local_poly2_scan_7_8';old=P/'runs/20261029_local_dcb_scan_7_8';rows=json.loads((old/'comparison.json').read_text());flags=[];auditrows=[];summaries=[]
for row in rows:
 m=row['mass']
 for model in ['bernstein2','chebyshev2']:
  d=R/model/f'm{m:.1f}';r=json.loads((d/'results.json').read_text());q=json.loads((d/'quality.json').read_text());a=json.loads((d/'audit.json').read_text());assert hashlib.sha256(pathlib.Path(a['source']).read_bytes()).hexdigest()==a['source_sha256'];assert r['limit_crossing']['verified'];assert abs(r['limit_crossing']['CLs']-.05)<.002
  for st in ['b','s']:
   if q[st]['status']!=0 or q[st]['covQual']!=3:flags.append(dict(model=model,mass=m,stage=st,status=q[st]['status'],covQual=q[st]['covQual']))
   pars={v['name']:v for v in q[st]['parameters']}
   for v in pars.values():
    if v['name'].startswith('local_') or v['name']=='r':
     if min(v['value']-v['minimum'],v['maximum']-v['value'])<1e-4*(v['maximum']-v['minimum']):flags.append(dict(model=model,mass=m,stage=st,bound=v['name']))
   for era in ['2022','2022EE','2023','2023BPix']:
    if model=='bernstein2':
     b1=pars['local_b1_'+era]['value'];b2=pars['local_b2_'+era]['value'];A=1-2*b1+b2;B=2*(b1-1);C=1;lo=0;hi=1
    else:
     c1=pars['local_c1_'+era]['value'];c2=pars['local_c2_'+era]['value'];A=2*c2;B=c1;C=1-c2;lo=-1;hi=1
    points=[lo,hi]
    if A!=0 and lo<-B/(2*A)<hi:points.append(-B/(2*A))
    minimum=min(A*x*x+B*x+C for x in points);assert minimum>0,(model,m,st,era,minimum)
    auditrows.append(dict(model=model,mass=m,stage=st,era=era,polynomial_minimum_unnormalized=minimum))
  for role in ['fitdiagnostics','significance','limits']:
   assert json.loads((d/role/'receipt.json').read_text())['exit_code']==0
   log=(d/role/'combine.log').read_text()
   if 'ERROR:' in log:flags.append(dict(model=model,mass=m,role=role,error_lines=[x for x in log.splitlines() if 'ERROR:' in x]))
  if r['r']['r']+r['r']['rHiErr']>=4.999 or r['r']['r']-r['r']['rLoErr']<=-4.999:flags.append(dict(model=model,mass=m,interval_bound=True))
  quant=[r['limits'][k] for k in ['expected_2p5','expected_16','expected_50','expected_84','expected_97p5']];assert min(quant)>0 and quant==sorted(quant)
  row.update({model+'_Z':r['Z'],model+'_p0':r['p0'],model+'_r':r['r']['r'],model+'_limit':r['limits']['observed'],model+'_expected':r['limits']['expected_50'],model+'_expected16':r['limits']['expected_16'],model+'_expected84':r['limits']['expected_84'],model+'_B_deviance':r['B_deviance'],model+'_B_gof_p':r['B_approx_p']})
  summaries.append(dict(model=model,mass=m,**r))
with (R/'comparison.csv').open('w') as f:w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
(R/'comparison.json').write_text(json.dumps(rows,indent=2));(R/'quality_flags.json').write_text(json.dumps(flags,indent=2));(R/'positivity_audit.json').write_text(json.dumps(auditrows,indent=2));x=np.array([r['mass'] for r in rows]);val=lambda k:np.array([r[k] for r in rows]);colors={'local':'#d66b17','global':'#256ca5','dcb':'#188548','bernstein2':'#9846b3','chebyshev2':'#c73535'};labels={'local':'Exp + Gaussiana, fit 7–8','global':'Cheb7 + Gaussiana, fit 6–12','dcb':'Exp + DCB, fit 7–8','bernstein2':'Bernstein2 + DCB, fit 7–8','chebyshev2':'Chebyshev2 + DCB, fit 7–8'}
for focused in [False,True]:
 models=['dcb','bernstein2','chebyshev2'] if focused else ['global','local','dcb','bernstein2','chebyshev2']
 for role in ['limits','pvalues','expected']:
  fig,ax=plt.subplots(figsize=(11,6.6))
  for model in models:
   style=dict(color=colors[model],label=labels[model],lw=1.8,marker='s' if model=='chebyshev2' else 'o',ms=6 if model=='chebyshev2' else 4,ls='--' if model=='chebyshev2' else '-',mfc='none' if model=='chebyshev2' else colors[model])
   if role=='limits':
    ax.plot(x,val(model+'_limit'),**style);ax.plot(x,val(model+'_expected'),':',color=colors[model],lw=1.1,alpha=.8)
   elif role=='expected':
    ax.plot(x,val(model+'_expected'),**style);ax.fill_between(x,val(model+'_expected16'),val(model+'_expected84'),color=colors[model],alpha=.10)
   else:ax.plot(x,val(model+'_p0'),**style)
  ax.set_yscale('log');ax.set_xlim(7,8);ax.set_xticks(x);ax.grid(alpha=.2,which='both');ax.set_xlabel('Massa ipotizzata del segnale [GeV]');ax.legend(fontsize=8.5)
  if role=='limits':ax.set_ylabel('Limite superiore su r, 95% CLs');ax.set_title('Aggiunta Bernstein2 e Chebyshev2 con DCB originale\nCurve con marcatori: osservato; puntinate: atteso mediano')
  elif role=='expected':ax.set_ylabel('Limite atteso su r, 95% CLs');ax.set_title('Sensibilità attesa: mediana e bande 68%')
  else:
   ax.set_ylabel('p-value locale p₀');ax.set_ylim(min(.001,*[val(k+'_p0').min()/2 for k in models]),1);ax.set_title('p-value: Bernstein2 e Chebyshev2 con DCB originale\nModelli locali: finestra fissa 7–8 GeV, quattro ere combinate')
   for z in [1,2,3,4]:
    p=norm.sf(z)
    if p>ax.get_ylim()[0]:ax.axhline(p,color='gray',ls=':',lw=.8);ax.text(7.985,p,f'{z}σ',ha='right',va='bottom',color='gray')
  fig.text(.5,.012,'Bernstein2 e Chebyshev2 descrivono entrambi un polinomio quadratico: le curve possono sovrapporsi.',ha='center',fontsize=8,color='#555555');fig.tight_layout(rect=[0,.035,1,1])
  for ext in ['png','pdf']:fig.savefig(R/f'{role}_{"local_DCB" if focused else "comparison"}.{ext}',dpi=150)
  plt.close(fig)
tr=''.join(f"<tr><td>{r['mass']:.1f}</td>"+''.join(f"<td>{r[model+'_Z']:.3f}</td><td>{r[model+'_limit']:.5f}</td><td>{r[model+'_expected']:.5f}</td>" for model in ['bernstein2','chebyshev2'])+'</tr>' for r in rows)
gallery=''.join(f"<h3>Massa {m/10:.1f} GeV</h3>"+''.join(f"<p><b>{model}</b> · <a href='{model}/m{m/10:.1f}/fits.pdf'>PDF fit</a> · <a href='{model}/m{m/10:.1f}/parameters.csv'>Parametri</a> · <a href='{model}/m{m/10:.1f}/fit_sum.png'>Somma/componenti</a></p><img loading='lazy' src='{model}/m{m/10:.1f}/fits.png'>" for model in ['bernstein2','chebyshev2']) for m in range(70,81))
(R/'index.html').write_text(f'''<!doctype html><meta charset="utf-8"><title>Polinomi grado2 + DCB, scan7–8</title><style>body{{font:17px system-ui;max-width:1150px;margin:35px auto;padding:20px;color:#263647}}img{{width:100%}}table{{border-collapse:collapse}}th,td{{padding:8px;border:1px solid #ddd}}</style><h1>Bernstein grado2 e Chebyshev grado2 + DCB</h1><p>11 masse7–8GeV,passo0,1, finestra di fit fissa7–8,4ere combinate,75bin originali perera, rcomune. 22workspace e66nuove esecuzioni Combine FitDiagnostics/AsymptoticLimits/Significance. ShapeDCBoriginale, code e nuisance MC preservate; stessa accettanza DCB7–8/6–12 del precedente scan. Tutti22limiti osservati verificati aCLs≈0,05; quantili attese positive e ordinate. Flag finali:{len(flags)} (<a href='quality_flags.json'>dettagli</a>); positività analitica controllata nei minimi dei polinomi per176configurazioni massa/era/fitB-S. Inputimmutati, hashverificati.</p><h2>Definizione delle due forme</h2><p>Bernstein2: b(u)∝(1−u)²+2b₁u(1−u)+b₂u², u=(m−7)/(8−7). b₀=1 elimina la scala comune ridondante; b₁,b₂∈[10⁻⁶,20] liberi, yieldseparato. Chebyshev2: b(t)∝1+c₁t+c₂(2t²−1), t=2(m−7)−1; c₁∈[−2,2],c₂∈[−1,1], liberi, yieldseparato. Due parametri di forma e una normalizzazione perera, nessuna prior dagli stessi dati. Sono due basi dello stesso spazio di polinomi quadratici, con domini diversi: se i minimi sono interni alla regione comune, ci aspettiamo forme e risultati coincidenti. Non sono due prove indipendenti di validità del fondo.</p><p>Fondo totale locale, incluse codeΥ senza alterare le shapeΥ sorgenti; confrontabile con l'esponenziale precedente. DCB e rate r1 identici per tutti imodelli; integraliaccettanza verificati controROOT nominale/scale±1/sigma±1 (440check). Minimiz.Strategia1,tol.001,FDrobustFit1/stepSize.001/robuststrategy1,r[-5,5]; limiti rMax1 fino a7,8 e rMax20 a7,9–8,0; accuratezza relativa0,001 e assoluta0,00001. Prima esecuzione archiviata in initial_attempt: ricalcolo dopo inizializzazione dell’osservabile dentro la finestra e precisione migliorata sui limiti. Nessun toyGOF/bias eseguito: in accordo con la richiesta sono rinviati.</p><h2>Curve aggiunte al confronto precedente</h2><a href='pvalues_comparison.pdf'>PDF p-value</a><img src='pvalues_comparison.png'><a href='limits_comparison.pdf'>PDF limiti</a><img src='limits_comparison.png'><h2>Confronto a parità di DCB e range: exp, Bernstein2, Chebyshev2</h2><img src='pvalues_local_DCB.png'><img src='limits_local_DCB.png'><img src='expected_local_DCB.png'><p>Il confronto bluCheb7Gauss usa il fit6–12 e il punto7,4 giàricontrollato nel precedente studio. I p-value sono locali e il modello non viene scelto minimizzandoli. Le bande attese riflettono l'ipotesi di ogni modello, non una validazione del bias.</p><h2>Tabella</h2><a href='comparison.csv'>CSV completo</a><table><tr><th>Massa</th><th>Z Bern2</th><th>Limite obs Bern2</th><th>Atteso Bern2</th><th>Z Cheb2</th><th>Limite obs Cheb2</th><th>Atteso Cheb2</th></tr>{tr}</table><h2>Tutti i fit, componenti e residui</h2>{gallery}''')
(R/'completion.json').write_text(json.dumps(dict(complete=True,masses=11,models=['bernstein2','chebyshev2'],combine_jobs=66,archived_initial_jobs=66,verified_observed_limits=22,quality_flags=len(flags),source_unchanged=True,acceptance_checks=440,positivity_checks=len(auditrows)),indent=2));print('FLAGS',flags)
for r in rows:print(r['mass'],[(k,r[k+'_Z'],r[k+'_limit'],r[k+'_expected']) for k in ['bernstein2','chebyshev2']])
print('BGOF',[(k,rows[0][k+'_B_deviance'],rows[0][k+'_B_gof_p']) for k in ['bernstein2','chebyshev2']])
