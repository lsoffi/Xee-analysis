from pathlib import Path
import json,numpy as np,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
P=Path.cwd();R=P/'runs/20261030_sr_cheb7_bern_bias';initial=[json.loads(l) for l in (R/'asimov.jsonl').read_text().splitlines()];ref=[json.loads(l) for l in (R/'stable_refined.jsonl').read_text().splitlines()];A=[x for x in initial if x['degree']<6]+ref;FROOT=P/'runs/20261030_sr_bernstein7_bias';fn='stable_asimov.jsonl' if (FROOT/'stable_asimov.jsonl').exists() and len((FROOT/'stable_asimov.jsonl').read_text().splitlines())==210 else 'asimov.jsonl';F=[json.loads(l) for l in (FROOT/fn).read_text().splitlines()];COL={3:'#9467bd',4:'#0072B2',5:'#E69F00',6:'#009E73',7:'#D55E00'};summary=[]
for d in range(3,8):
 aa=[x for x in A if x['degree']==d];good=[x for x in aa if x['valid']];worst=max(good,key=lambda x:abs(x['bias_sigma']));summary.append(dict(degree=d,valid=len(good),total=len(aa),maximum_bias_sigma=abs(worst['bias_sigma']),worst_mass=worst['mass'],worst_injection=worst['injection'],passes_asimov=len(good)==42 and abs(worst['bias_sigma'])<.2))
(R/'summary.json').write_text(json.dumps(summary,indent=2));(R/'final_asimov.json').write_text(json.dumps(A))
fig,axes=plt.subplots(2,2,figsize=(13,8),sharex=True)
for row,inj in enumerate(['zero','expected']):
 for col,ax in enumerate(axes[row]):
  for d in range(3,8):
   aa=sorted([x for x in A if x['degree']==d and x['injection']==inj],key=lambda x:x['mass']);ax.plot([x['mass'] for x in aa],[x['bias_sigma'] if x['valid'] else np.nan for x in aa],'o-',ms=3,color=COL[d],label=f'Bernstein {d}');bad=[x for x in aa if not x['valid']];ax.plot([x['mass'] for x in bad],[x['bias_sigma'] for x in bad],'x',color=COL[d],ms=6)
  ax.axhspan(-.2,.2,color='grey',alpha=.2);ax.axhline(0,color='k',lw=.6);ax.set_title(('r vero = 0' if inj=='zero' else 'r vero = limite atteso di riferimento')+(' · zoom' if col else ''));ax.grid(alpha=.2);ax.set_xlabel('Massa ipotizzata [GeV]');ax.set_ylabel('(r fit − r vero) / σ del fit')
  if col:ax.set_ylim(-.8,.8)
axes[0,0].legend(fontsize=8,ncol=2);fig.suptitle('Generatore Chebyshev 7 → fit Bernstein 3–7, sempre SR+CR\nAsimov; μ Υ liberi, σ/code MC fisse; CR definitiva invariata; fascia grigia ±0,2σ');fig.tight_layout();fig.savefig(R/'asimov_bias.png',dpi=160);fig.savefig(R/'asimov_bias.pdf');plt.close(fig)
fig,axes=plt.subplots(1,2,figsize=(13,4.6),sharex=True,sharey=True)
for j,(data,title,model) in enumerate([(F,'Generatore Bernstein 7','Chebyshev'),(A,'Generatore Chebyshev 7','Bernstein')]):
 ax=axes[j]
 for d in [4,5,6,7]:
  aa=sorted([x for x in data if x['degree']==d and x['injection']=='zero'],key=lambda x:x['mass']);ax.plot([x['mass'] for x in aa],[x['bias_sigma'] if x['valid'] else np.nan for x in aa],'o-',ms=3,color=COL[d],label=f'{model} {d}')
 ax.axhspan(-.2,.2,color='grey',alpha=.2);ax.axhline(0,color='k',lw=.6);ax.set_ylim(-.8,.8);ax.grid(alpha=.2);ax.set_title(title);ax.set_xlabel('Massa [GeV]');ax.set_ylabel('Bias / σ del fit');ax.legend(fontsize=8)
fig.suptitle('Confronto inverso: grado minimo 6 nei due test Asimov richiesti\nZoom r vero=0; non una validazione della copertura dei limiti');fig.tight_layout();fig.savefig(R/'comparison.png',dpi=160);plt.close(fig)
# Generator GoF from the existing identical simultaneous Chebyshev7 model/toys.
gofrows=[json.loads(l) for l in (P/'runs/20261030_sr_degree_selection/toys.jsonl').read_text().splitlines()];G=json.loads((R/'generators.json').read_text());D=sum(g['D'] for g in G);gg=[]
for k in range(100):
 fs=[x for x in gofrows if x['generator_degree']==7 and x['degree']==7 and x['toy']==k]
 if len(fs)==4 and all(f['valid'] for f in fs):gg.append(sum(f['D'] for f in fs))
ex=sum(x>=D for x in gg);gof=dict(observed=D,total=100,valid=len(gg),exceedances=ex,failure_bounds=[ex/100,(ex+100-len(gg))/100]);(R/'gof_summary.json').write_text(json.dumps(gof))
lines=['# Generatore Chebyshev7 → fit Bernstein3–7','','## Configurazione','SR unica6–9,093333GeV, CRdefinitiva9,093333–10,506667GeV con Chebyshev4. Quattroere. μ/normΥ1S/2Sliberi e condivisi SR/CR;σecodeMCfisse. CRinvariata. Nessunfitgrado8.','GeneratoreCheb7 B-only fitted simultaneamente ai dati.21masse6,5–8,5GeV, passo0,1. SegnaleDCBMC nominale presente inSR eCR conrstesso. Per ogni massa fitBernstein3–7 suAsimovr0 e r=limiteattesomediano dello scanstoricoSR7/CR3 usato esclusivamente comeampiezza di iniezione.210fitS+B. Nonsono nuovi limiti; non è una simulazione di copertura; nuisanceMCsegnale nominali, nonprofilate inquestotestcondizionale.','Bernstein n+1 coefficienti positivi di intensità, senza un ulteriore yield ridondante;CRconcontinuoindipendente. Integrali esatti bin13,333MeV. Implementazioneiminuit, nonCombine.','Gradi6/7: ricontrollo uniforme84fit, stessa likelihood valutata conlog1p/serie per evitare cancellazioni numeriche vicino alminimo, precisioneMinuit1e-12/strategy2/tol.001. Tutti84validi. Variazione massima rrispettoalprimotentativo<0,005σ; iniziali e retryconservati. Gradi3–5 hanno alcuni fallimenti espliciti e biasgiàsoprasoglia nei puntivalidi.','', '|Grado Bernstein|FitAsimovvalidi|Massimo bias/σ nei validi|Passa tutti i punti a0,2σ?|','|---|---:|---:|---|']
for s in summary:lines.append(f"|{s['degree']}|{s['valid']}/{s['total']}|{s['maximum_bias_sigma']:.4f}|{'Sì' if s['passes_asimov'] else 'No'}|")
lines+=['','## Conclusione','Bernstein6 è il grado minimo che supera il criterio Asimov in questo confronto inverso. Nella direzioneBernstein7→Chebyshev il minimo èChebyshev6. Le due conclusioni sono condizionate ai generatori richiesti e alle21masse provate. Non annullano i risultati degli altri scenari precedenti.','Bernstein eChebyshev allo stesso grado descrivono lo stesso spazio polinomiale; la positività dei coefficientiBernstein restringe il dominio. PertantoBernstein6/Chebyshev6 sono principalmente un confronto di parametrizzazione e vincoli, NONdue famiglie di fondo indipendenti.','Il controllo con toy fluttuatiBernstein7→Cheb è separato nella directory runs/20261030_sr_bernstein7_bias. Non viene presentato come un nuovo test fluttuato delladirezioneinversa. I profili delle nuisanceΥ e la copertura dei limiti rimangono check separati.']
(R/'REPORT.md').write_text('\n'.join(lines));html='<meta charset="utf-8"><style>body{font:17px system-ui;max-width:1300px;margin:auto}img{max-width:100%}</style><h1>Chebyshev 7 → Bernstein 3–7</h1><p><a href="REPORT.md">Report</a> · <a href="summary.json">Tabella</a></p><img src="asimov_bias.png"><img src="comparison.png">';(R/'index.html').write_text(html);print(summary);print(gof)
