import json,numpy as np,csv
from pathlib import Path
p=Path('outputs');rows=[];metas=[]
for c,label in [('free','Liberi'),('limit1','±1σ pre-fit'),('limit2','±2σ pre-fit')]:
 tag='degree7_to10_'+c;m=json.load(open(p/f'combine_2022_{tag}_pre_post.json'));metas.append((label,m));r=m['parameters'][0]
 row=dict(case=label,r=r['s'],r_lo=r['s_lo'],r_hi=r['s_hi'],status_b=m['b']['status'],covQual_b=m['b']['covQual'],status_s=m['s']['status'],covQual_s=m['s']['covQual'])
 for st in ['fit_b','fit_s']:
  a=np.genfromtxt(p/f'combine_2022_{tag}_{st}_shapes.csv',delimiter=',',names=True);n=a['data'];mu=a['total'];row[st+'_deviance']=float(2*np.sum(mu-n+n*np.log(n/mu)))
 rows.append(row)
 s=Path('outputs/run_fitdiagnostics_2022.sh').read_text().replace('Xee_ee_4_2022_chebyshev.root','Xee_ee_4_2022_chebyshev_'+tag+'.root').replace('_inclusive_2022_chebyshev -v','_inclusive_2022_chebyshev_'+tag+' -v');out=p/f'run_fitdiagnostics_2022_{tag}.sh';out.write_text(s);out.chmod(0o755)
with open(p/'degree7_to10_bounded_summary.csv','w') as f:
 wr=csv.DictWriter(f,fieldnames=rows[0].keys());wr.writeheader();wr.writerows(rows)
lines=['# Fit originale grado7 senza rebin, range6–10: liberi, ±1σ, ±2σ','','2022, range6–10GeV,300 bin originali,32014 eventi. Chebyshev originale grado7, Υ(1S) congelata; altre opzioni originali Combine, strategia0/robustFit1/r tra−30 e30. Stessa likelihood del fit originale senza rebin.','', 'Gli intervalli dei sette coefficienti sono rigidi e centrati ai valori originali del workspace, usando RooRealVar.getError() salvati. Non sono vincoli gaussiani indipendenti, né un ellissoide che incorpora le correlazioni.','', '| Caso | r (68% robustFit) | Devianza B | Devianza S+B | Status B/S+B | covQual B/S+B |','|---|---:|---:|---:|---|---|']
for a in rows:lines.append(f"| {a['case']} | {a['r']:.6f} −{abs(a['r_lo']):.5f}/+{a['r_hi']:.5f} | {a['fit_b_deviance']:.3f} | {a['fit_s_deviance']:.3f} | {a['status_b']}/{a['status_s']} | {a['covQual_b']}/{a['covQual_s']} |")
lines+=['','Il fit è numericamente instabile: libero covQual2 in B e S+B; ±1σ covQual2 in B e3 in S+B; ±2σ covQual3 in B e2 in S+B. Le matrici covQual2 sono forzate positive definite e gli errori HESSE non sono affidabili. I valori e gli intervalli di r sono risultati della procedura, non una validazione del fit.','', 'Nel fit S+B ±1σ raggiungono i bordi t0(+1), t2(−1), t5(+1), t6(−1). Con ±2σ t5 è prossimo al bordo superiore.','', 'La base Chebyshev e la definizione dei yield restano sul range originale6–12. Le normalizzazioni nel range6–10 sono corrette tramite integrali. Il taglio10 interrompe la regione del picco Υ(2S).']
for stage,title in [('b','Solo fondo'),('s','Segnale + fondo')]:
 lines+=['','## '+title,'','| Parametro | Pre-fit | Post libero | Post ±1σ | Post ±2σ |','|---|---:|---:|---:|---:|']
 maps=[{a['parameter']:a for a in m['parameters']} for _,m in metas]
 for name,a in maps[0].items():
  def fmt(v):
   if v[stage+'_err'] is None:return f"{v[stage]:.6g} fisso"
   if name=='r' and stage=='s':return f"{v[stage]:.6g} −{abs(v[stage+'_lo']):.4g}/+{v[stage+'_hi']:.4g}"
   return f"{v[stage]:.7g} ± {v[stage+'_err']:.3g}"
  lines.append(f"| `{name}` | {a['prefit']:.7g} | {fmt(maps[0][name])} | {fmt(maps[1][name])} | {fmt(maps[2][name])} |")
lines+=['','## Spostamenti dei coefficienti nel fit S+B','','| Parametro | Libero / σ salvata | ±1σ / σ salvata | ±2σ / σ salvata |','|---|---:|---:|---:|']
for name,a in maps[0].items():
 if not name.startswith('t'):continue
 shifts=[(mp[name]['s']-a['prefit'])/a['prefit_saved_error'] for mp in maps];lines.append(f"| `{name}` | {shifts[0]:+.4f} | {shifts[1]:+.4f} | {shifts[2]:+.4f} |")
p.joinpath('degree7_to10_bounded_report.md').write_text('\n'.join(lines)+'\n')
print('\n'.join(lines[:12]))
