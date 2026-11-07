from pathlib import Path
import json,csv,time,os
import numpy as np
from scipy.stats import norm
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
P=Path.cwd();R=P/'runs/20261030_sr6_cr4_scan';rows=sorted([json.loads(p.read_text()) for p in R.glob('m*/results.json')],key=lambda x:x['mass']);(R/'summary.json').write_text(json.dumps(rows,indent=2));good=[r for r in rows if r['limit_crossing']['verified']];x=[k/10 for k in range(60,91)];by={r['mass']:r for r in good}
def values(key):return [by[m]['limits'][key] if m in by else np.nan for m in x]
fig,ax=plt.subplots(figsize=(10,5.5))
for lo,hi,col,lab in [('expected_2p5','expected_97p5','#F6D743','Atteso ±2σ (provvisorio)'),('expected_16','expected_84','#69BE46','Atteso ±1σ (provvisorio)')]:ax.fill_between(x,values(lo),values(hi),color=col,label=lab)
for key,col,ls,lab in [('expected_50','black','--','Atteso mediano'),('observed','#0072B2','-','Osservato')]:ax.plot(x,values(key),color=col,ls=ls,label=lab)
for r in good:
 ax.vlines(r['mass'],r['limits']['expected_2p5'],r['limits']['expected_97p5'],color='#F6D743',lw=8,zorder=0)
 ax.vlines(r['mass'],r['limits']['expected_16'],r['limits']['expected_84'],color='#69BE46',lw=8,zorder=1)
 ax.plot(r['mass'],r['limits']['expected_50'],'_',color='black')
for r in good:ax.plot(r['mass'],r['limits']['observed'],'o' if r['numerical_ok'] else 'x',color='#0072B2' if r['numerical_ok'] else '#D55E00')

if not all(r.get('expected_bands_verified',False) for r in good):ax.text(.02,.97,'Bande attese provvisorie: estremi superiori non verificati',transform=ax.transAxes,va='top',color='firebrick')
ax.set_yscale('log');ax.set_xlabel('Massa ipotizzata [GeV]');ax.set_ylabel('Limite 95% CLs su r');ax.set_title(f'SR envelope Cheb6/Bern6 / CR Chebyshev4 + DCB + Υ — {len(rows)}/31 punti\nSR fissa: 6–9,093 GeV; CR: 9,093–10,507 GeV');ax.legend(loc='lower right');ax.grid(alpha=.2);fig.tight_layout()
for ext in ['png','pdf']:
 tmp=R/f'.limits_{os.getpid()}.{ext}';fig.savefig(tmp,dpi=160);tmp.replace(R/f'limits.{ext}')
plt.close(fig)
fig,ax=plt.subplots(figsize=(10,5.5));pd={r['mass']:r for r in rows};ax.plot(x,[max(pd[m]['p0'],1e-300) if m in pd else np.nan for m in x],color='#0072B2')
for r in rows:ax.plot(r['mass'],max(r['p0'],1e-300),'o' if r['numerical_ok'] else 'x',color='#0072B2' if r['numerical_ok'] else '#D55E00')
for z in [1,2,3,4,5]:
 if norm.sf(z)>=min([r['p0']*.3 for r in rows]+[norm.sf(3)*.7]):ax.axhline(norm.sf(z),color='gray',ls=':',alpha=.5);ax.text(9.03,norm.sf(z),f'{z}σ',fontsize=9)
ax.set_yscale('log');ax.set_xlim(5.95,9.15);ax.set_ylim(min([r['p0']*.3 for r in rows]+[norm.sf(3)*.7]),.65);ax.set_xlabel('Massa ipotizzata [GeV]');ax.set_ylabel('p-value locale');ax.set_title('SR envelope Cheb6/Bern6 / CR Chebyshev4 + DCB + Υ — SR e CR simultanee\nCode DCB incluse anche nella CR');ax.grid(alpha=.2);fig.tight_layout()
for ext in ['png','pdf']:
 tmp=R/f'.pvalues_{os.getpid()}.{ext}';fig.savefig(tmp,dpi=160);tmp.replace(R/f'pvalues.{ext}')
plt.close(fig)
with (R/'summary.csv').open('w') as f:
 keys=['mass','SR_low','SR_high','r','Z','p0','limit_observed','limit_expected','numerical_ok','max_CR_signal_fraction'];w=csv.DictWriter(f,keys);w.writeheader()
 for r in rows:w.writerow(dict(mass=r['mass'],SR_low=r['range'][0],SR_high=r['range'][1],r=r['r'],Z=r['Z'],p0=r['p0'],limit_observed=r['limits']['observed'],limit_expected=r['limits']['expected_50'],numerical_ok=r['numerical_ok'],max_CR_signal_fraction=max(r['CR_signal_fraction'].values())))
body='<meta charset="utf-8"><meta http-equiv="refresh" content="60"><style>body{font:17px system-ui;max-width:1400px;margin:30px auto}img{width:100%}td,th{padding:8px;border-bottom:1px solid #ddd}</style><h1>Scan6–9: SR envelope Cheb6/Bern6 / CR Chebyshev4 + DCB + Υ</h1>'
body+=f'<p>{len(rows)}/31 punti; {sum(r["numerical_ok"] for r in rows)} con controlli numerici superati. Aggiornato {time.strftime("%H:%M:%S")}. FitDiagnostics, AsymptoticLimits e Significance reali Combine. SegnaleDCBoriginale, nuisanceoriginali; Envelope Cheb6/Bern6 inSR,Chebyshev4 inCR,indipendenti perregione/era; normalizzazioniΥcondivise/libere. SR eCRdisgiunte, CR9,093333–10,506667; binningoriginale. SegnaleDCBincluso ancheinCR. SR fissa6–9,093333GeV per tutte le masse. Bias Asimov e GoF con 100 toy completati; copertura non verificata; pvaluelocali, selezione delmodelloesplorativa.</p><p>Croci arancioni: problemi numerici non risolti; non trattarle come risultati validati. <a href="summary.csv">CSV</a> · <a href="summary.json">JSON</a></p><p style="color:firebrick"><strong>Bande attese provvisorie:</strong> gli estremi superiori sono troncati dalla ricerca numerica del crossing. Conservate su richiesta; ricalcolo interrotto. Non usarle per conclusioni quantitative.</p><img src="limits.png"><img src="pvalues.png"><table><tr><th>Massa</th><th>SR</th><th>Z</th><th>p0</th><th>Status/cov B;S</th><th>Fit</th></tr>'
for r in rows:
 m=f'{r["mass"]:.1f}';q='; '.join(f'{k}:{v["status"]}/{v["covQual"]}' for k,v in r['quality'].items());body+=f'<tr><td>{m}</td><td>{r["range"][0]:.3f}–{r["range"][1]:.3f}</td><td>{r["Z"]:.3f}</td><td>{r["p0"]:.4g}</td><td>{q}</td><td><a href="#m{m}">SR</a> · <a href="m{m}/fit_CR.png">CR</a> · <a href="m{m}/joint/parameters.csv">parametri</a></td></tr>'
body+='</table>'
if (R/'bias.png').exists():body+='<h2>Bias Asimov (parziale finché la campagna è in corso)</h2><img src="bias.png">'
if (R/'gof.png').exists():body+='<h2>Goodness of fit B-only: 100 toy</h2><img src="gof.png"><a href="gof/summary.json">Risultati GoF</a>'
for r in rows:
 m=f'{r["mass"]:.1f}';body+=f'<h2 id="m{m}">Massa {m} GeV</h2><img loading="lazy" src="m{m}/fit_SR.png"><p><a href="m{m}/fit_SR.pdf">PDF SR</a> · <a href="m{m}/fit_CR.pdf">PDF CR</a></p>'
(R/'index.html').write_text(body);print('points',len(rows),'flags',[r['mass'] for r in rows if not r['numerical_ok']]);print('maximum',max(rows,key=lambda r:r['Z']) if rows else None)
