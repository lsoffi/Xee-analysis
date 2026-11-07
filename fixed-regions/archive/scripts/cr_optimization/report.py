from study import *
from scipy.stats import beta
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import csv
S=json.loads((R/'scan.json').read_text());T=[json.loads(s) for s in (R/'toys.jsonl').read_text().splitlines()];G=json.loads((R/'bernstein_generator.json').read_text());summ=[]
def interval(k,n):return [0 if k==0 else float(beta.ppf(.025,k,n-k+1)),1 if k==n else float(beta.ppf(.975,k+1,n-k))]
for r in S:
 li=round((r['lo']-6)*75);ui=round((r['hi']-6)*75);deg=r['degree'];go=[t for t in T if t['kind']=='gof' and t['li']==li and t['ui']==ui and t['degree']==deg];bi=[t for t in T if t['kind']=='bias' and t['li']==li and t['ui']==ui and t['degree']==deg]
 if not go:continue
 good=[t for t in go if t['valid']];k=sum(t['D']>=r['D'] for t in good);rec={k:v for k,v in r.items() if k!='fits'};rec.update(gof_p=k/len(good),gof_CP95=interval(k,len(good)),gof_valid=len(good),gof_total=len(go))
 if deg<3:
  nxt=next(z for z in S if z['lo']==r['lo'] and z['hi']==r['hi'] and z['degree']==deg+1);delta=r['D']-nxt['D'];paired=[t for t in go if t['valid'] and t['valid_next']];kk=sum(t['delta']>=delta for t in paired);rec.update(LR_delta=delta,LR_p=kk/len(paired),LR_CP95=interval(kk,len(paired)),LR_valid=len(paired),LR_negative=sum(t['delta']<-.001 for t in paired))
 if bi:
  valid=[t for t in bi if t['valid']];p=np.array([t['prediction'] for t in valid]);er=np.array([t['error'] for t in valid]);truth=np.array(valid[0]['truth']);dif=p-truth;bias=np.mean(dif,axis=0);spread=np.std(p,axis=0,ddof=1);sem=spread/np.sqrt(len(p));meanerr=np.mean(er,axis=0);coverage=np.mean(abs(dif)<=er,axis=0);rec.update(bias_valid=len(valid),bias_total=len(bi),truth=truth.tolist(),bias=bias.tolist(),bias_MCerror=sem.tolist(),bias_over_sigma=(bias/meanerr).tolist(),bias_over_sigma_MCerror=(sem/meanerr).tolist(),coverage68_Wald=coverage.tolist(),coverage_CP95=[interval(int(sum(abs(dif[:,j])<=er[:,j])),len(valid)) for j in range(2)],mean_error=meanerr.tolist(),toy_spread=spread.tolist())
 summ.append(rec)
(R/'toy_summary.json').write_text(json.dumps(summ,indent=2))
# precision/stability across the entire grid
fig,ax=plt.subplots(2,2,figsize=(12,8),sharex=True)
for deg,col in [(1,'#0072B2'),(2,'#E69F00'),(3,'#009E73')]:
 a=[s for s in S if s['degree']==deg];xx=[z['width'] for z in a]
 for j in range(2):
  ax[0,j].scatter(xx,[z['prediction'][j] for z in a],color=col,label=f'Chebyshev {deg}',s=30);ax[1,j].scatter(xx,[100*z['error'][j]/z['prediction'][j] for z in a],color=col,s=30)
for j,label in enumerate(['Υ previste in 6–9,093 GeV','Υ previste in 10,507–12 GeV']):
 ax[0,j].set_title(label);ax[0,j].set_ylabel('Eventi Υ(1S) + Υ(2S)');ax[1,j].set_ylabel('Errore Hesse relativo [%]');ax[1,j].set_xlabel('Ampiezza CR [GeV]')
for a in ax.flat:a.grid(alpha=.2)
ax[0,0].legend();fig.suptitle('48 fit CR: stessa shape Υ, diversi bordi e gradi\nLe finestre strette possono passare il GoF ma perdere identificabilità');fig.tight_layout();fig.savefig(R/'window_scan.png',dpi=150);fig.savefig(R/'window_scan.pdf');plt.close(fig)
if all('bias' in s for s in summ):
 fig,ax=plt.subplots(2,2,figsize=(13,8),sharex=True);ranges=sorted(set((s['lo'],s['hi']) for s in summ),key=lambda x:x[1]-x[0],reverse=True);xx=np.arange(len(ranges))
 for deg,col in [(1,'#0072B2'),(2,'#E69F00'),(3,'#009E73')]:
  aa=[next(s for s in summ if s['lo']==lo and s['hi']==hi and s['degree']==deg) for lo,hi in ranges]
  for j in range(2):
   ax[0,j].errorbar(xx+.12*(deg-2),[s['bias_over_sigma'][j] if s['bias_valid']/s['bias_total']>=.95 else np.nan for s in aa],yerr=[s['bias_over_sigma_MCerror'][j] for s in aa],fmt='o-',color=col,label=f'Grado {deg}');ax[1,j].plot(xx+.12*(deg-2),[s['coverage68_Wald'][j] if s['bias_valid']/s['bias_total']>=.95 else np.nan for s in aa],'o-',color=col)
 for j in range(2):
  ax[0,j].axhspan(-.2,.2,color='grey',alpha=.15,label='Obiettivo ±0,2σ' if j==0 else None);ax[0,j].axhline(0,color='k',lw=.7);ax[0,j].set_ylabel('Bias / errore medio');ax[1,j].axhline(.6827,color='k',ls='--');ax[1,j].set_ylim(0,1);ax[1,j].set_ylabel('Copertura intervallo Hesse ±1σ');ax[1,j].set_xticks(xx);ax[1,j].set_xticklabels([f'{a:.3f}\n{b:.3f}' for a,b in ranges]);ax[1,j].set_xlabel('Bordi CR [GeV]');ax[0,j].set_title(['Υ estrapolate sotto la CR','Υ estrapolate sopra la CR'][j])
 for a in ax.flat:a.grid(alpha=.2)
 ax[0,0].legend();fig.suptitle('200 toy per finestra/grado — generatore alternativo Bernstein 4\nBarre: errore Monte Carlo; intervalli Hesse; omesso grado 3 a 9,4–10,2 (83/200 fit invalidi)');fig.tight_layout();fig.savefig(R/'bias_coverage.png',dpi=150);fig.savefig(R/'bias_coverage.pdf');plt.close(fig)
# all observed fit plots
for r in S:
 fig,axes=plt.subplots(2,2,figsize=(11,7));deg=r['degree']
 for i,ax in enumerate(axes.flat):
  A,n,mask=matrices(r['lo'],r['hi'],deg,i);f=r['fits'][i];p=np.array(f['p']);xc=(edges[:-1]+edges[1:])[mask]/2;ax.errorbar(xc,n,yerr=np.sqrt(n),fmt='.',color='k',ms=3,label='Dati');ax.plot(xc,f['mu'],label='Totale');ax.plot(xc,A[:,:deg+1]@p[:deg+1],'--',label='Continuo');ax.plot(xc,A[:,-2]*p[-2],label='Υ(1S)');ax.plot(xc,A[:,-1]*p[-1],label='Υ(2S)');ax.set_title(ER[i]['era']);ax.set_xlabel('Massa [GeV]');ax.set_ylabel('Eventi / bin originale');ax.grid(alpha=.2)
 axes[0,0].legend(fontsize=8);fig.suptitle(f"CR {r['lo']:.3f}–{r['hi']:.3f} GeV, Chebyshev {deg}");fig.tight_layout();name=f"fit_{round((r['lo']-6)*75)}_{round((r['hi']-6)*75)}_d{deg}.png";fig.savefig(R/name,dpi=120);plt.close(fig)
lines=['# Ottimizzazione rapida della CR Υ','','48 fit: 16 finestre × Chebyshev gradi 1, 2, 3; quattro ere indipendenti; shape Υ originali fissate, normalizzazioni libere [0,10]. Nessun segnale nel test CR-only.','Likelihood Poisson binned con integrali esatti per bin (13,333 MeV), minimizzazione Minuit/iminuit, continuo parametrizzato linearmente. Positività controllata su griglia durante il fit e agli estremi e punti stazionari dopo. È un test dedicato della CR, non un nuovo scan Combine di limiti/p-value.','I nuovi fit integrano le PDF nei bin, mentre i precedenti Combine valutavano la densità ai centri. Cross-check CR baseline grado3: scale Υ1 differiscono meno di 0,0005; scale Υ2 meno di 0,005. Non confondere piccole differenze di implementazione con cambiamenti fisici.','Intervalli/errori Hesse: possibili problemi ai bordi, in particolare Υ2=0 in alcune ere. Le coperture qui misurate si riferiscono solo agli intervalli Hesse sulle somme delle rese estrapolate, non a intervalli profilati o a limiti sul segnale.','Ogni fit dei toy usa due seed, risultato minimo; tutti i fallimenti sono conteggiati e conservati. GoF e LR condizionali ai parametri fitted: 200 toy per grado/finestra; LR fra gradi consecutivi con gli stessi dati toy. Non si confrontano NLL di finestre differenti.','Generatore alternativo positivo Bernstein4 fittato alla CR larga, poi ritagliato: stessa verità e stessi toy appaiati tra finestre/gradi. 200 pseudoesperimenti indipendenti, ciascuno rifittato nelle 15 configurazioni. Un solo generatore alternativo non certifica robustezza generale. Bernstein e Chebyshev a pari grado non sono alternative indipendenti: qui il generatore ha grado4.','Le rese estrapolate sono riferite a intervalli FISSI 6–9,0933 e 10,5067–12: si confronta la stessa quantità senza modificare la SR durante la scelta CR. Non misura la sensibilità del successivo fit SR+CR.','Soglia bias esplorativa |bias|<0,2 volte errore medio. Riportato errore Monte Carlo: un risultato vicino alla soglia non certifica il passaggio. GoF>0,05 non sufficiente per scegliere la finestra.','','|CR|Grado|p GoF|p LR al grado dopo|Bias/σ sotto|Bias/σ sopra|Copertura sotto/sopra|Toy validi bias|','|---|---:|---:|---:|---:|---:|---|---|']
for s in summ:
 if 'bias' not in s:continue
 lines.append(f"|{s['lo']:.3f}–{s['hi']:.3f}|{s['degree']}|{s['gof_p']:.3f}|{s.get('LR_p','—')}|{s['bias_over_sigma'][0]:+.2f} ± {s['bias_over_sigma_MCerror'][0]:.2f}|{s['bias_over_sigma'][1]:+.2f} ± {s['bias_over_sigma_MCerror'][1]:.2f}|{s['coverage68_Wald'][0]:.3f} / {s['coverage68_Wald'][1]:.3f}|{s['bias_valid']}/{s['bias_total']}|")
(R/'REPORT.md').write_text('\n'.join(lines));html='<meta charset="utf-8"><style>body{font:17px system-ui;max-width:1200px;margin:auto}img{max-width:100%}td{padding:8px}</style><h1>Ottimizzazione CR Υ</h1><p><a href="REPORT.md">Report e limiti del test</a> · <a href="toy_summary.json">Risultati toy</a> · <a href="scan.json">Tutti i parametri</a></p><img src="window_scan.png"><img src="bias_coverage.png"><h2>Tutti i fit</h2>'
for r in S:
 name=f"fit_{round((r['lo']-6)*75)}_{round((r['hi']-6)*75)}_d{r['degree']}.png";html+=f'<details><summary>{r["lo"]:.3f}–{r["hi"]:.3f}, grado {r["degree"]}</summary><img src="{name}"></details>'
(R/'index.html').write_text(html)
print('report',len(T),'toy rows',flush=True)
