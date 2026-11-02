from pathlib import Path
import json,numpy as np,matplotlib,hashlib,csv
matplotlib.use('Agg')
import matplotlib.pyplot as plt
P=Path.cwd();R=P/'runs/20261030_sr_upper_mass_bias';SIG=json.loads((R/'signals.json').read_text());COL={3:'#9467bd',4:'#0072B2',5:'#E69F00',6:'#009E73',7:'#D55E00'};rows={};summ=[]
for direction in ['forward','reverse']:
 a=[json.loads(l) for l in (R/(direction+'.jsonl')).read_text().splitlines()];assert len(a)==60;rows[direction]=a
 for x in a:
  sig=next(s for s in SIG if s['mass']==x['mass']);assert max(abs(x['signal_CR_fraction'][ei]-q['fraction_CR']) for ei,q in enumerate(sig['eras']))<1e-8
 for d in range(3,8):
  aa=[x for x in a if x['degree']==d];valid=[x for x in aa if x['valid'] and not x['r_at_boundary']];summ.append(dict(direction=direction,degree=d,total=len(aa),valid=len(valid),maximum_bias_sigma=max(abs(x['bias_sigma']) for x in valid),passes=len(valid)==12 and max(abs(x['bias_sigma']) for x in valid)<.2))
for s in SIG:assert hashlib.sha256(Path(s['source']).read_bytes()).hexdigest()==s['sha256']
(R/'summary.json').write_text(json.dumps(summ,indent=2))
fig,axes=plt.subplots(2,2,figsize=(12,8),sharex=True,sharey=True)
for col,direction in enumerate(['forward','reverse']):
 for row,inj in enumerate(['zero','one']):
  ax=axes[row,col];model='Chebyshev' if direction=='forward' else 'Bernstein';generator='Bernstein 7' if direction=='forward' else 'Chebyshev 7'
  for d in range(3,8):
   aa=sorted([x for x in rows[direction] if x['degree']==d and x['injection']==inj],key=lambda x:x['mass']);ax.plot([x['mass'] for x in aa],[x['bias_sigma'] if x['valid'] and not x['r_at_boundary'] else np.nan for x in aa],'o-',ms=4,color=COL[d],label=f'{model} {d}')
  ax.axhspan(-.2,.2,color='grey',alpha=.2);ax.axhline(0,color='k',lw=.6);ax.set_ylim(-.8,.8);ax.set_xticks(np.arange(8.5,9.01,.1));ax.set_xlabel('Massa ipotizzata [GeV]');ax.set_ylabel('(r fit − r vero) / σ del fit');ax.set_title(f'Generatore {generator}; r vero = '+('0' if inj=='zero' else '1'));ax.grid(alpha=.2);ax.legend(fontsize=8,ncol=2)
fig.suptitle('Bias Asimov 8,5–9 GeV — zoom ±0,8σ\nDCB del segnale in SR e CR, unico r; CR definitiva invariata; valori fuori scala non visibili');fig.tight_layout();fig.savefig(R/'bias_zoom.png',dpi=160);fig.savefig(R/'bias_zoom.pdf');plt.close(fig)
# All degrees without clipping.
fig,axes=plt.subplots(1,2,figsize=(12,4.5),sharex=True,sharey=True)
for col,direction in enumerate(['forward','reverse']):
 ax=axes[col];model='Chebyshev' if direction=='forward' else 'Bernstein'
 for d in range(3,8):
  aa=sorted([x for x in rows[direction] if x['degree']==d and x['injection']=='zero'],key=lambda x:x['mass']);ax.plot([x['mass'] for x in aa],[x['bias_sigma'] if x['valid'] and not x['r_at_boundary'] else np.nan for x in aa],'o-',ms=4,color=COL[d],label=f'{model} {d}')
 ax.axhspan(-.2,.2,color='grey',alpha=.2);ax.axhline(0,color='k',lw=.6);ax.set_xlabel('Massa [GeV]');ax.set_ylabel('Bias / σ del fit');ax.set_title('Generatore '+('Bernstein 7' if direction=='forward' else 'Chebyshev 7'));ax.grid(alpha=.2);ax.legend(fontsize=8,ncol=2)
fig.suptitle('8,5–9 GeV: tutti i gradi senza tagli verticali — Asimov r=0');fig.tight_layout();fig.savefig(R/'bias_all_degrees.png',dpi=160);plt.close(fig)
# Whole scanned mass range, only zero injection (identical definition at every mass).
oldforward=[json.loads(l) for l in (P/'runs/20261030_sr_bernstein7_bias/stable_asimov.jsonl').read_text().splitlines()];oldreverse=json.loads((P/'runs/20261030_sr_cheb7_bern_bias/final_asimov.json').read_text());fig,axes=plt.subplots(1,2,figsize=(13,4.8),sharex=True,sharey=True)
for col,(direction,old) in enumerate([('forward',oldforward),('reverse',oldreverse)]):
 ax=axes[col];model='Chebyshev' if direction=='forward' else 'Bernstein';a=[x for x in old if x['mass']<8.5]+rows[direction]
 for d in range(3,8):
  aa=sorted([x for x in a if x['degree']==d and x['injection']=='zero'],key=lambda x:x['mass']);ax.plot([x['mass'] for x in aa],[x['bias_sigma'] if x['valid'] and not x.get('r_at_boundary',False) else np.nan for x in aa],'o-',ms=3,color=COL[d],label=f'{model} {d}')
 ax.axhspan(-.2,.2,color='grey',alpha=.2);ax.axvline(8.5,color='grey',ls=':',lw=.8);ax.axhline(0,color='k',lw=.6);ax.set_ylim(-.8,.8);ax.set_xlabel('Massa [GeV]');ax.set_ylabel('Bias / σ del fit');ax.set_title('Generatore '+('Bernstein 7' if direction=='forward' else 'Chebyshev 7'));ax.grid(alpha=.2);ax.legend(fontsize=8,ncol=2)
fig.suptitle('Scan completo 6,5–9 GeV — bias Asimov r=0, zoom ±0,8σ\nSegnale DCB incluso anche nella CR; fascia grigia ±0,2σ');fig.tight_layout();fig.savefig(R/'bias_full_scan.png',dpi=160);plt.close(fig)
# Fractions and precision to distinguish bias from loss of sensitivity.
fig,axes=plt.subplots(1,2,figsize=(11,4.2))
for ei,e in enumerate(['2022','2022EE','2023','2023BPix']):axes[0].plot([s['mass'] for s in SIG],[100*s['eras'][ei]['fraction_CR'] for s in SIG],'o-',label=e)
for direction,label in [('forward','Chebyshev 6'),('reverse','Bernstein 6')]:
 aa=sorted([x for x in rows[direction] if x['degree']==6 and x['injection']=='zero'],key=lambda x:x['mass']);axes[1].plot([x['mass'] for x in aa],[x['error'] for x in aa],'o-',label=label)
axes[0].set_ylabel('Segnale nella CR [% del totale 6–12]');axes[1].set_ylabel('Errore Hesse su r');
for ax in axes:ax.set_xlabel('Massa [GeV]');ax.grid(alpha=.2);ax.legend(fontsize=8)
fig.suptitle('Contenimento del segnale e precisione: due aspetti distinti dal bias');fig.tight_layout();fig.savefig(R/'signal_fraction_precision.png',dpi=150);plt.close(fig)
lines=['# Estensione bias Asimov a8,5–9GeV','','SR6–9,093333 eCRdefinitiva9,093333–10,506667GeV, Chebyshev4nellaCR,μ/normΥcondivisi liberi,σ/codeMCfisse. Quattroere. DCBsegnaleoriginale inENTRAMBEle regioni, unico r,accettanze ricalcolate. Era giàinclusa nellaCR nei test precedenti.','6masse8,5–9a passi0,1; generatoriBernstein7 eChebyshev7 fitted B-only ai dati, fit incrociatiChebyshev3–7/Bernstein3–7. Asimovr0 e r1:120fit. r1 è una ampiezza fissa, NONunlimiteattesocalcolato perquesti punti. Non sono nuovi limiti/p-value né una validazione della copertura.','NuisanceMCsegnale nominali inquestotestcondizionale. LikelihoodPoisson conintegralibin esatti, stabilelog1p/serie, iminuitprecision1e-12/strategy2/tol.001. FrazioniDCB confrontateconintegraliROOT a<1e-8, hash6inputimmutati.','', '|Generatore → fit|Grado|Fit validi|Massimo bias/σ su8,5–9 e due iniezioni|','|---|---:|---:|---:|']
for q in summ:lines.append(f"|{'Bernstein7 → Cheb' if q['direction']=='forward' else 'Cheb7 → Bernstein'}|{q['degree']}|{q['valid']}/{q['total']}|{q['maximum_bias_sigma']:.4f}|")
lines+=['','## Conclusione','Ilgrado6mantiene ilbiasAsimov sotto0,2σ nell’estensione fino9GeV perentrambi igeneratoririchiesti. Questa èla conclusionedeltestcondizionale, non unagaranziauniversale.','La frazioneDCB nellaCR cresce da1,24–1,52% a8,5GeV a15,87–17,17% a9GeV. L’erroreHesse sur delgrado6cresce da≈0,17 a≈1,07 neltest a r0. Biascontenuto nonimplica sensibilitàcostante. Non attribuireautomaticamentetutta la crescita alleΥ: lecorrelazionir–Υ sono riportate neiJSON, ancheilcontinuopartecipa alprofiling.','I profili vicino albordo, la copertura dei limiti e iltrattamentocompleto delle nuisanceMC rimangono verifiche separate. Nessun cambiamento dellaCRdefinitiva. Bernstein eChebyshev a parigrado restano parametrizzazioni dello stesso spaziopolinomiale con vincoli diversi, non duefamiglieindipendenti.']
(R/'REPORT.md').write_text('\n'.join(lines));(R/'completion.json').write_text(json.dumps(dict(masses=6,range=[8.5,9.0],fits=120,all_valid=sum(x['valid'] and not x['r_at_boundary'] for a in rows.values() for x in a),signal_in_SR_and_CR=True,injections=[0,1],input_unchanged=True,integral_check_tolerance=1e-8),indent=2));html='<meta charset="utf-8"><style>body{font:17px system-ui;max-width:1300px;margin:auto}img{max-width:100%}</style><h1>Bias fino a9GeV — DCB inSR+CR</h1><p><a href="REPORT.md">Report</a> · <a href="summary.json">Tabella</a></p><img src="bias_zoom.png"><img src="bias_all_degrees.png"><img src="bias_full_scan.png"><img src="signal_fraction_precision.png">';(R/'index.html').write_text(html);print(summ)
