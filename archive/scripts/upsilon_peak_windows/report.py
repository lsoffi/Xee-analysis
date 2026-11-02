from pathlib import Path
import json,csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
P=Path.cwd();R=P/'runs/20261029_upsilon_peak_windows';common=json.loads((R/'common_windows.json').read_text());cores=json.loads((R/'core_windows.json').read_text());sf=json.loads((R/'signal_fractions.json').read_text());sr=json.loads((R/'SR_overlap.json').read_text());eras=['2022','2022EE','2023','2023BPix'];current=[9.093333333333334,10.506666666666668];allc=[('CR attuale',[current],'#777777')]+[(f'Contenimento {round(100*r["target_fraction"])}%',r['union'],'#D55E00') for r in common]+[(f'Nucleo ±{r["sigma_core_multiplier"]:g}σ',r['union'],'#0072B2') for r in reversed(cores)]
fig,ax=plt.subplots(figsize=(11,5))
for i,(label,regions,color) in enumerate(allc):
 for a,b in regions:ax.plot([a,b],[i,i],color=color,lw=15,solid_capstyle='butt');ax.text(a,i+.19,f'{a:.3f}',ha='left',fontsize=9);ax.text(b,i+.19,f'{b:.3f}',ha='right',fontsize=9)
ax.set_yticks(range(len(allc)),[x[0] for x in allc]);ax.invert_yaxis();ax.set_xlabel('Massa [GeV]');ax.set_title('CR comuni alle quattro ere, allineate ai bin originali\nGli intervalli sovrapposti sono uniti: nessun doppio conteggio');ax.grid(axis='x',alpha=.2);ax.set_xlim(7.95,10.7);fig.tight_layout();fig.savefig(R/'windows.png',dpi=170);fig.savefig(R/'windows.pdf');plt.close(fig)
fig,ax=plt.subplots(figsize=(11,5.5))
for key,label,col in [('current','CR attuale','#777777'),('68','Finestre al 68% (unione)','#D55E00'),('core1.5','Due nuclei ±1,5σ','#009E73'),('core1.0','Due nuclei ±1σ','#0072B2')]:
 masses=sorted({r['mass'] for r in sf});lo=[];hi=[]
 for m in masses:
  a=[100*r['signal_fraction_in_CR'] for r in sf if r['mass']==m and r['candidate']==key];lo.append(min(a));hi.append(max(a))
 ax.plot(masses,np.mean([lo,hi],axis=0),color=col,label=label);ax.fill_between(masses,lo,hi,color=col,alpha=.15)
ax.set_xlim(8,11);ax.set_ylim(0,100);ax.set_xlabel('Massa ipotizzata del segnale DCB [GeV]');ax.set_ylabel('Frazione di segnale nella CR [%]');ax.set_title('Segnale nella CR: da includere nel fit simultaneo\nFasce: minimo–massimo fra le quattro ere, non incertezze');ax.legend();ax.grid(alpha=.2);fig.tight_layout();fig.savefig(R/'signal_in_CR.png',dpi=170);fig.savefig(R/'signal_in_CR.pdf');plt.close(fig)
fig,ax=plt.subplots(figsize=(11,5.5))
for key,label,col in [('current','CR attuale','#777777'),('core1.5','Due nuclei ±1,5σ','#009E73'),('core1.0','Due nuclei ±1σ','#0072B2')]:
 a=[r for r in sr if r['candidate']==key];m=[r['mass'] for r in a];lo=[100*r['retained_signal_SR_min'] for r in a];hi=[100*r['retained_signal_SR_max'] for r in a];ax.plot(m,np.mean([lo,hi],axis=0),color=col,label=label);ax.fill_between(m,lo,hi,color=col,alpha=.15)
ax.axhline(95,color='gray',ls=':');ax.set_xlim(8,10);ax.set_ylim(0,100);ax.set_xlabel('Massa ipotizzata [GeV]');ax.set_ylabel('Segnale contenuto nella sola SR privata dei bin CR [%]');ax.set_title('Finestre SR precedenti al 95%, con esclusione delle CR\nIl segnale nelle CR resta utilizzabile nel fit simultaneo');ax.legend();ax.grid(alpha=.2);fig.tight_layout();fig.savefig(R/'signal_in_masked_SR.png',dpi=170);fig.savefig(R/'signal_in_masked_SR.pdf');plt.close(fig)
summary=[]
for r in cores:
 key=r['label'];a=[x for x in sr if x['candidate']==key];small=dict(label=key,regions=r['union'],width=r['union_width'],gap=r['gap'],recovered_vs_current_CR=r['recovered_vs_current_CR'],largest_existing_unmodified_disjoint_mass=max(x['mass'] for x in a if x['disjoint']),coverage_own_peak={p['process']:[min(p['coverage'].values()),max(p['coverage'].values())] for p in r['pieces']},examples={})
 for m in [8.5,8.8,9.,9.2,9.5,9.8,10.]:
  x=next(x for x in a if x['mass']==m);f=[x['signal_fraction_in_CR'] for x in sf if x['candidate']==key and x['mass']==m];small['examples'][str(m)]=dict(signal_in_CR=[min(f),max(f)],signal_in_masked_SR=[x['retained_signal_SR_min'],x['retained_signal_SR_max']])
 summary.append(small)
(R/'summary.json').write_text(json.dumps(dict(definition='Nominal original DCB shapes conditioned on6–12GeV; no fitting, no new model selection. Common envelopes across eras, outward bin alignment. Overlapping CR candidates merged. Core sigma is the Gaussian-core parameter, not sigma_eff.',core_candidates=summary,containment_candidates=common,signal_mass_range=[6.5,11],existing_SR_mass_range=[6.5,10]),indent=2))
body='<meta charset="utf-8"><style>body{font:17px system-ui;max-width:1300px;margin:30px auto}img{width:100%}td,th{padding:8px;border-bottom:1px solid #ddd}</style><h1>Finestre attorno ai picchi Υ e contaminazione DCB</h1><p>Calcolo geometrico sulle shape originali, senza nuovi fit. Tutte le frazioni sono normalizzate sul dominio6–12GeV. Input disponibili per46masse6,5–11; finestreSR precedenti disponibili finoa10. Due CR al68% non sono separate: le code basse causano una sovrapposizione. Per isolare i nuclei occorre accettare un contenimento inferiore.</p><img src="windows.png"><table><tr><th>Candidato</th><th>CR1</th><th>CR2</th><th>Contenimento Υ1 nel proprio picco</th><th>Υ2</th><th>Spettro recuperato</th></tr>'
for r in cores:
 p1,p2=r['pieces'];cov=lambda p:f'{100*min(p["coverage"].values()):.1f}–{100*max(p["coverage"].values()):.1f}%';body+=f'<tr><td>±{r["sigma_core_multiplier"]:g}σ</td><td>{p1["low"]:.3f}–{p1["high"]:.3f}</td><td>{p2["low"]:.3f}–{p2["high"]:.3f}</td><td>{cov(p1)}</td><td>{cov(p2)}</td><td>{r["recovered_vs_current_CR"]:.3f}GeV</td></tr>'
body+='</table><p>Entrambe le Υ devono essere presenti in entrambe le CR: le code della2S contribuiscono anche sotto la1S. Le quote di segnale nellaCR non sono eventi persi se la DCB è inclusa nella likelihood con lo stesso r. Questo calcolo non determina la precisione delle normalizzazioni né la separabilità segnale/Υ: richiedono fit e test di bias/iniezione. La baselineCR9,1–10,5 non è stata cambiata.</p><img src="signal_in_CR.png"><img src="signal_in_masked_SR.png"><p><a href="signal_fractions.csv">Frazioni per massa/era</a> · <a href="per_era_windows.csv">Finestre Υ per era</a> · <a href="SR_overlap.csv">Intersezioni con SR precedenti</a> · <a href="validation.json">Validazione ROOT</a> · <a href="summary.json">Riepilogo numerico</a></p>';(R/'index.html').write_text(body);print(json.dumps(summary,indent=2))
