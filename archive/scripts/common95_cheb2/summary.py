from pathlib import Path
import json,csv,time,html
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
P=Path('/Users/imb/Documents/Codex/Analisi-statistica-Xee');R=P/'runs/20261029_common95_cheb2';rows=sorted([json.loads(f.read_text()) for f in R.glob('m*/results.json')],key=lambda r:r['mass']);R.mkdir(exist_ok=True);(R/'summary.json').write_text(json.dumps(rows,indent=2))
fig,ax=plt.subplots(figsize=(10,5.5))
band_rows=[r for r in rows if 'limits' in r and r.get('limit_crossing',{}).get('verified')]
for r in band_rows:
 vals=[r['limits'][k] for k in ['expected_2p5','expected_16','expected_50','expected_84','expected_97p5']]
 assert all(np.isfinite(v) and v>0 for v in vals) and all(a<=b for a,b in zip(vals,vals[1:])),r['mass']
band_x=[r['mass'] for r in band_rows]
ax.fill_between(band_x,[r['limits']['expected_2p5'] for r in band_rows],[r['limits']['expected_97p5'] for r in band_rows],color='#F6D743',label='Atteso ±2σ',zorder=1)
ax.fill_between(band_x,[r['limits']['expected_16'] for r in band_rows],[r['limits']['expected_84'] for r in band_rows],color='#69BE46',label='Atteso ±1σ',zorder=2)
for key,label,ls in [('observed','Osservato','-'),('expected_50','Atteso mediano','--')]:
 good=[r for r in rows if 'limits' in r and r.get('limit_crossing',{}).get('verified')];xx=[r['mass'] for r in good];yy=[r['limits'][key] for r in good];ax.plot(xx,yy,ls,color='black' if key=='expected_50' else '#0072B2',label=label,zorder=3)
 if key=='observed':
  for ok,marker,col,lab in [(True,'o','#0072B2','Controlli numerici di base superati'),(False,'x','#D55E00','Covarianza / minimi non stabili')]:
   g=[r for r in good if r.get('numerical_ok',r['fit_ok'])==ok];ax.scatter([r['mass'] for r in g],[r['limits'][key] for r in g],marker=marker,color=col,label=lab,zorder=4)
ax.set_yscale('log');ax.set_xlabel('Massa ipotizzata [GeV]');ax.set_ylabel('Limite 95% CLs su r');ax.set_title(f'Chebyshev 2 + DCB + Υ(1S) + Υ(2S) — {len(rows)}/36 punti\nFinestre comuni al 95%, allineate ai bin originali');ax.grid(alpha=.2);ax.legend(fontsize=9);fig.tight_layout()
for ext in ['png','pdf']:fig.savefig(R/f'limits.{ext}',dpi=160)
plt.close(fig)
fig,ax=plt.subplots(figsize=(10,5.5));g=[r for r in rows if 'p0' in r];ax.plot([r['mass'] for r in g],[max(r['p0'],1e-300) for r in g],color='#0072B2',lw=1)
for ok,marker,col,lab in [(True,'o','#0072B2','Controlli numerici di base superati'),(False,'x','#D55E00','Covarianza / minimi non stabili')]:
 a=[r for r in g if r.get('numerical_ok',r['fit_ok'])==ok];ax.scatter([r['mass'] for r in a],[max(r['p0'],1e-300) for r in a],marker=marker,color=col,label=lab,zorder=4)
from scipy.stats import norm
for z in [1,2,3,4,5]:ax.axhline(norm.sf(z),ls=':',color='gray',alpha=.5);ax.text(10.02,norm.sf(z),f'{z}σ',fontsize=8)
ax.set_yscale('log');ax.set_xlabel('Massa ipotizzata [GeV]');ax.set_ylabel('p-value locale');ax.set_title('Chebyshev 2 — risultati preliminari, non validati con toy');ax.set_xlim(6.45,10.15);ax.grid(alpha=.2);ax.legend(fontsize=9);fig.tight_layout()
for ext in ['png','pdf']:fig.savefig(R/f'pvalues.{ext}',dpi=160)
plt.close(fig)
with (R/'summary.csv').open('w') as f:
 w=csv.writer(f);w.writerow(['mass','low','high','r','Z','p0','limit_observed','limit_expected','fit_ok','CLs_crossing_verified'])
 for r in rows:w.writerow([r['mass'],*r['range'],r['r']['r'],r.get('Z'),r.get('p0'),r.get('limits',{}).get('observed'),r.get('limits',{}).get('expected_50'),r.get('numerical_ok',r['fit_ok']),r.get('limit_crossing',{}).get('verified')])
body='<!doctype html><meta charset="utf-8"><meta http-equiv="refresh" content="60"><title>Chebyshev2 finestre 95%</title><style>body{font:17px system-ui;max-width:1400px;margin:30px auto}img{width:100%}td,th{padding:7px;border:1px solid #ddd}table{border-collapse:collapse}</style>'
body+=f'<h1>Chebyshev2 + DCB + Υ(1S) + Υ(2S)</h1><p>{len(rows)}/36 risultati; {sum(r["fit_ok"] for r in rows)} con status/covarianza B/S soddisfacenti. I simboli arancioni includono anche discrepanze fra minimi ottenuti con inizializzazioni diverse. Nessuna validazione toy. Ultimo aggiornamento {time.strftime("%H:%M:%S")}.</p><p>Finestre nominali comuni al 95% nel dominio 6–12 GeV, estese verso l’esterno ai bin originali. Υ libere; continuo positivo. FitDiagnostics: Hesse, senza MINOS/robustFit. Simboli arancioni: risultati con problemi numerici. Collegamenti fra punti solo guida visiva.</p><img src="limits.png"><img src="pvalues.png"><p><a href="summary.csv">Tabella CSV</a></p>'
body+='<table><tr><th>Massa</th><th>Range</th><th>Status/cov B; S</th><th>Fit</th></tr>'
for r in rows:
 m=f'{r["mass"]:.1f}';qs='; '.join(f'{st}: {v["status"]}/{v["covQual"]}' for st,v in r['quality'].items());body+=f'<tr><td>{m}</td><td>{r["range"][0]:.3f}–{r["range"][1]:.3f}</td><td>{qs}</td><td><a href="#m{m}">vedi</a> · <a href="m{m}/fit.pdf">PDF</a></td></tr>'
body+='</table>'
for r in rows:
 m=f'{r["mass"]:.1f}';body+=f'<h2 id="m{m}">Massa {m} GeV</h2><img loading="lazy" src="m{m}/fit.png">'
(R/'index.html').write_text(body)
print('RESULTS',len(rows),'FIT_GOOD',sum(r['fit_ok'] for r in rows))
