import pathlib,json,csv,hashlib,re
import numpy as np,uproot
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import norm
P=pathlib.Path('/Users/imb/Documents/Codex/Analisi-statistica-Xee');R=P/'runs/20261029_local_center78_scan';LO=7.293333333333333;HI=8.306666666666667;models=['exp','bernstein2','chebyshev2'];colors={'exp':'#188548','bernstein2':'#9846b3','chebyshev2':'#c73535'};labels={'exp':'Esponenziale + DCB','bernstein2':'Bernstein2 + DCB','chebyshev2':'Chebyshev2 + DCB'};rows=[];flags=[];positive=[];old=[]
def olddir(k,m):return P/('runs/20261029_local_dcb_scan_7_8' if k=='exp' else 'runs/20261029_local_poly2_scan_7_8/'+k)/f'm{m:.1f}'
for k in models:
 for i in range(70,81):
  m=i/10;r=json.loads((olddir(k,m)/'results.json').read_text());old.append(dict(model=k,mass=m,Z=r['Z'],p0=r['p0'],limit=r['limits']['observed'],expected=r['limits']['expected_50'],r=r['r']['r'],rLo=r['r']['rLoErr'],rHi=r['r']['rHiErr']))
 for i in range(73,84):
  m=i/10;d=R/k/f'm{m:.1f}';r=json.loads((d/'results.json').read_text());a=json.loads((d/'audit.json').read_text());q=json.loads((d/'quality.json').read_text());assert hashlib.sha256(pathlib.Path(a['source']).read_bytes()).hexdigest()==a['source_sha256'];assert r['limit_crossing']['verified'];assert abs(r['limit_crossing']['CLs']-.05)<.002
  ex=[r['limits'][s] for s in ['expected_2p5','expected_16','expected_50','expected_84','expected_97p5']];assert min(ex)>0 and ex==sorted(ex)
  for st in ['b','s']:
   if q[st]['status']!=0 or q[st]['covQual']!=3:flags.append(dict(model=k,mass=m,stage=st,status=q[st]['status'],covQual=q[st]['covQual']))
   par={v['name']:v for v in q[st]['parameters']}
   for v in par.values():
    if v['name'].startswith('local_') or v['name']=='r':
     if min(v['value']-v['minimum'],v['maximum']-v['value'])<1e-4*(v['maximum']-v['minimum']):flags.append(dict(model=k,mass=m,stage=st,bound=v['name']))
   for era in ['2022','2022EE','2023','2023BPix']:
    if k=='exp':minimum=np.exp(min(LO*par['local_slope_'+era]['value'],HI*par['local_slope_'+era]['value']))
    else:
     if k=='bernstein2':b1=par['local_b1_'+era]['value'];b2=par['local_b2_'+era]['value'];A=1-2*b1+b2;B=2*(b1-1);C=1;lo=0;hi=1
     else:c1=par['local_c1_'+era]['value'];c2=par['local_c2_'+era]['value'];A=2*c2;B=c1;C=1-c2;lo=-1;hi=1
     pts=[lo,hi]
     if A!=0 and lo<-B/(2*A)<hi:pts.append(-B/(2*A))
     minimum=min(A*t*t+B*t+C for t in pts)
    assert minimum>0;positive.append(dict(model=k,mass=m,era=era,stage=st,minimum=float(minimum)))
  for role in ['fitdiagnostics','significance','limits']:
   assert json.loads((d/role/'receipt.json').read_text())['exit_code']==0
   log=(d/role/'combine.log').read_text()
   if 'ERROR:' in log:flags.append(dict(model=k,mass=m,role=role,error_lines=[x for x in log.splitlines() if 'ERROR:' in x]))
  if r['r']['r']+r['r']['rHiErr']>=29.999 or r['r']['r']-r['r']['rLoErr']<=-29.999:flags.append(dict(model=k,mass=m,interval_bound=True))
  rows.append(dict(model=k,mass=m,Z=r['Z'],p0=r['p0'],limit=r['limits']['observed'],expected=r['limits']['expected_50'],expected16=r['limits']['expected_16'],expected84=r['limits']['expected_84'],r=r['r']['r'],rLo=r['r']['rLoErr'],rHi=r['r']['rHiErr'],B_deviance=r['B_deviance'],B_df=r['B_df'],B_gof_approx_p=r['B_approx_p'],data_events=a['total_events']))
for name,items in [('results_scan',rows),('previous_window',old)]:
 (R/(name+'.json')).write_text(json.dumps(items,indent=2))
 with (R/(name+'.csv')).open('w') as f:w=csv.DictWriter(f,fieldnames=list(items[0]));w.writeheader();w.writerows(items)
(R/'quality_flags.json').write_text(json.dumps(flags,indent=2));(R/'positivity_audit.json').write_text(json.dumps(positive,indent=2))
for role in ['pvalues','limits','expected']:
 fig,ax=plt.subplots(figsize=(11.5,6.8))
 for k in models:
  for group,data in [('nuova',rows),('precedente',old)]:
   a=[r for r in data if r['model']==k];key={'pvalues':'p0','limits':'limit','expected':'expected'}[role];color=colors[k];ax.plot([r['mass'] for r in a],[r[key] for r in a],color=color,lw=1.8 if group=='nuova' else 1.1,ls=('-.' if k=='chebyshev2' else '-') if group=='nuova' else ':',marker=('s' if k=='chebyshev2' else 'o') if group=='nuova' else None,mfc='none' if k=='chebyshev2' else color,ms=5,label=labels[k]+(' — 7,293–8,307' if group=='nuova' else ' — 7–8'))
 ax.axvline(7.8,color='gray',lw=.9,ls='--');ax.set_yscale('log');ax.set_xlim(7,8.3);ax.set_xticks(np.arange(7,8.31,.1));ax.grid(alpha=.18,which='both');ax.set_xlabel('Massa ipotizzata [GeV]');ax.legend(fontsize=8.5,loc='best')
 if role=='pvalues':
  ax.set_ylim(.001,1);ax.set_ylabel('p-value locale p₀');ax.set_title('Stabilità rispetto alla finestra: DCB originale in tutti i modelli')
  for z in [1,2,3]:p=norm.sf(z);ax.axhline(p,color='gray',ls=':',lw=.7);ax.text(8.29,p,f'{z}σ',ha='right',va='bottom',color='gray')
 else:ax.set_ylabel('Limite '+('osservato' if role=='limits' else 'atteso mediano')+' su r, 95% CLs');ax.set_title('Confronto finestra centrata a 7,8 e finestra precedente')
 fig.text(.5,.012,'76 bin/era, finestra 7,2933–8,3067 GeV. A 8,3 GeV: covQual 2 nel FitDiagnostics Bernstein2.',ha='center',fontsize=8);fig.tight_layout(rect=[0,.035,1,1])
 for ext in ['png','pdf']:fig.savefig(R/f'{role}_windows.{ext}',dpi=150)
 plt.close(fig)
eras=json.loads((R/'era_quality.json').read_text())
for a in eras:
 d=pathlib.Path(a['path']).parent
 for role in ['fitdiagnostics','significance']:assert json.loads((d/role/'receipt.json').read_text())['exit_code']==0
 with uproot.open(next((d/'significance').glob('higgsCombine*.root'))) as f:a['Z']=float(f['limit']['limit'].array(library='np')[0]);a['p0']=float(norm.sf(a['Z']))
 with uproot.open(d/'fitdiagnostics/fitDiagnostics_era.root') as f:t=f['tree_fit_sb'];a['rLo']=float(t['rLoErr'].array(library='np')[0]);a['rHi']=float(t['rHiErr'].array(library='np')[0])
 if a['s_status']!=0 or a['s_covQual']!=3 or a['b_status']!=0 or a['b_covQual']!=3:flags.append(dict(era=a['era'],model=a['model'],s_status=a['s_status'],s_covQual=a['s_covQual']))
(R/'individual_era_results.json').write_text(json.dumps(eras,indent=2))
with (R/'individual_era_results.csv').open('w') as f:w=csv.DictWriter(f,fieldnames=list(eras[0]));w.writeheader();w.writerows(eras)
fig,ax=plt.subplots(figsize=(9,5.5));en=['2022','2022EE','2023','2023BPix']
for j,k in enumerate(models):
 a=[next(a for a in eras if a['model']==k and a['era']==e) for e in en];y=np.arange(4)+(j-1)*.16;ax.errorbar([v['r'] for v in a],y,xerr=[[v['rLo'] for v in a],[v['rHi'] for v in a]],fmt='o',ms=4,color=colors[k],label=labels[k],capsize=3)
ax.axvline(0,color='gray',ls='--');ax.set_yticks(np.arange(4),en);ax.invert_yaxis();ax.set_xlabel('r best fit e intervallo 68%');ax.set_title('Massa 7,8 GeV: fit indipendenti delle singole ere\nFinestra 7,2933–8,3067 GeV');ax.grid(alpha=.2);ax.legend();fig.tight_layout()
for ext in ['png','pdf']:fig.savefig(R/f'individual_eras_m7p8.{ext}',dpi=150)
plt.close(fig)
tr=''
for k in models:
 n=next(r for r in rows if r['model']==k and r['mass']==7.8);o=next(r for r in old if r['model']==k and r['mass']==7.8);tr+=f"<tr><td>{labels[k]}</td><td>{o['Z']:.3f} → {n['Z']:.3f}</td><td>{o['r']:.5f} → {n['r']:.5f}</td><td>{o['limit']:.5f} → {n['limit']:.5f}</td><td>{o['expected']:.5f} → {n['expected']:.5f}</td></tr>"
etr=''.join(f"<tr><td>{a['model']}</td><td>{a['era']}</td><td>{a['r']:.5f} −{a['rLo']:.5f}/+{a['rHi']:.5f}</td><td>{a['Z']:.3f}</td><td>{a['p0']:.4f}</td><td>{a['b_status']}/{a['b_covQual']}; {a['s_status']}/{a['s_covQual']}</td></tr>" for a in eras)
gallery=''.join(f"<h3>{k}, massa {i/10:.1f}</h3><p><a href='{k}/m{i/10:.1f}/fits.pdf'>PDF</a> · <a href='{k}/m{i/10:.1f}/parameters.csv'>Parametri</a> · <a href='{k}/m{i/10:.1f}/fit_sum.png'>Somma e componenti</a></p><img loading='lazy' src='{k}/m{i/10:.1f}/fits.png'>" for i in range(73,84) for k in models)
(R/'index.html').write_text(f'''<!doctype html><meta charset="utf-8"><title>Stabilità della struttura a7,8</title><style>body{{font:17px system-ui;max-width:1150px;margin:35px auto;padding:20px;color:#263647}}img{{width:100%}}table{{border-collapse:collapse}}th,td{{padding:8px;border:1px solid #ddd}}</style><h1>Stabilità cambiando finestra: centro7,8GeV</h1><p>Nuova finestra7,293333–8,306667GeV centrataesattamentea7,8: 76bin/era originali da13,333MeV. La finestra7,3–8,3 tagliava bin a metà: sono inclusi ibininteri che coprono il range, senza interpolare gli eventi. 11ipotesi7,3–8,3passo0,1;4ere combinate. Modelli exp, Bernstein2 eChebyshev2 con DCBoriginale. Medesime shapeDCB e nuisance, rate r1 riferito a6–12, accettanza integrata sulla nuova finestra e profilata. Fondo totale locale include le codeΥ, le shapeΥoriginali nonsono alterate. Le basi quadratiche sono definite sulla nuova finestra.</p><p>99job di scan +24job FitDiagnostics/Significance perereindipendenti a7,8. Tutti i33limiti verificatiCLs≈0,05; inputchecksumverificati; positività dei fondi verificata analiticamente. Flag qualità:{len(flags)}. Strategia1,tol.001;FDrobust1,stepSize.0001,robuststrategy1,r[-30,30];limiti rMax5 fino8,0 e50 da8,1,accuratezza relativa.001/assoluta1e-5. Il range numerico di r è ampliato rispetto allo studio precedente, senza nuovi vincoli fisici del modello.</p><p><b>Avvertenza numerica a8,3GeV:</b> il fitS+B Bernstein2 ha covQual2 anche dopo strategia2/tolleranza1e-5. Gli altri fit combinati hanno status0/covQual3; tutti i fit separati perera a7,8 hanno status0/covQual3. Limite e p-value a8,3 sono riportati come diagnostica nominale, con il flag conservato. I fit a7,8 non hanno questo problema. I primi tentativi tecnici di estrazione delle ere sono archiviati: i risultati finali usano una singola categoria e gli stessi pesi originali.</p><h2>A7,8: precedente7–8 → nuova finestra centrata</h2><table><tr><th>Modello</th><th>Z</th><th>r</th><th>Limite osservato</th><th>Atteso mediano</th></tr>{tr}</table><p>La preferenza positiva non rimane comune a tutti imodelli: con l'esponenziale il bestfit cambia segno; con iquadratici resta debole. Sono gli stessi eventi inlargaparte, quindi le finestre sono correlate: questo confronto non è una confermaindipendente e nondimostra da solo che l'eccesso fosse spurio. Nessun nuovo testbias oGOFtoy è stato eseguito.</p><img src='pvalues_windows.png'><img src='limits_windows.png'><img src='expected_windows.png'><h2>Singole ere a7,8</h2><p>Fitindipendenti di ciascuna era, stessi modelli e vincoli originali; qui r viene stimato separatamente. Le Z non si sommano inquadratura per ricostruire automaticamente il combinato: nuisance condivise e rcomune cambiano il problema. Gli intervalli mostrano se ciascuna era è compatibile con zero e con gli altri valori.</p><img src='individual_eras_m7p8.png'><table><tr><th>Modello</th><th>Era</th><th>r e68%</th><th>Z</th><th>p₀</th><th>B;S status/cov</th></tr>{etr}</table><p><a href='results_scan.csv'>CSV nuovo scan</a> · <a href='previous_window.csv'>CSV precedente</a> · <a href='individual_era_results.csv'>CSV ere</a> · <a href='quality_flags.json'>Flag</a></p><h2>Tutti i fit combinati, componenti e residui</h2>{gallery}''')
(R/'quality_flags.json').write_text(json.dumps(flags,indent=2));(R/'completion.json').write_text(json.dumps(dict(complete=True,scan_jobs=99,era_jobs=24,masses=11,models=models,verified_observed_limits=33,flags=len(flags),source_unchanged=True,fit_range=[LO,HI],bins_per_era=76),indent=2));print('FLAGS',flags)
print('AT7.8',[r for r in rows if r['mass']==7.8]);print('ERAS',[(a['model'],a['era'],a['r'],a['Z'],a['s_covQual']) for a in eras])
