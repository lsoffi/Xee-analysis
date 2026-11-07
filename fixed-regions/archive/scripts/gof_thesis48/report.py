from pathlib import Path
import json,csv,re,hashlib
import numpy as np,uproot
from scipy.stats import beta
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
P=Path.cwd();R=P/'runs/20261030_gof_thesis48';eras=['2022','2022EE','2023','2023BPix'];rows=[]
for task in json.loads((R/'tasks.json').read_text()):
 d=Path(task['path'])
 if not (d/'toys/receipt.json').exists():continue
 if json.loads((d/'toys/receipt.json').read_text())['exit_code']!=0:continue
 def read(name):
  with uproot.open(next((d/name).glob('higgsCombine*.root'))) as f:return f['limit']['limit'].array(library='np')
 obs=read('observed');toys=read('toys');assert len(obs)==1 and np.all(np.isfinite(obs));finite=np.isfinite(toys);N=len(toys);ex=int(np.sum(toys>=obs[0]));p=ex/N if N else None;ci=[float(beta.ppf(.025,ex,N-ex+1)) if ex else 0,float(beta.ppf(.975,ex+1,N-ex)) if ex<N else 1]
 logs='\n'.join((d/x/'combine.log').read_text() for x in ['observed','generate','toys']);flags=re.findall(r'^.*(?:fit failed|failed to converge|invalid NLL|ERROR|nan|infinite).*$' ,logs,re.M|re.I)
 audit=json.loads((d/'audit.json').read_text());unchanged=hashlib.sha256(Path(audit['source']).read_bytes()).hexdigest()==audit['source_sha256'];assert unchanged
 row=dict(era=task['era'],mass=task['mass'],observed=float(obs[0]),N=N,N_finite=int(sum(finite)),exceedances=ex,p=p,CP95=ci,valid=N==100 and bool(np.all(finite)) and not flags,log_flags=flags,seed=task['seed'],source_unchanged=unchanged)
 (d/'result.json').write_text(json.dumps(row,indent=2));rows.append(row)
rows.sort(key=lambda a:(a['era'],a['mass']));(R/'summary.json').write_text(json.dumps(rows,indent=2))
with (R/'summary.csv').open('w') as f:
 keys=['era','mass','observed','N','N_finite','exceedances','p','CP95','valid','seed'];w=csv.DictWriter(f,keys,extrasaction='ignore');w.writeheader();w.writerows(rows)
colors=['#0072B2','#D55E00','#009E73','#CC79A7'];grid=np.arange(60,91)/10
for era in eras:
 a=[v for v in rows if v['era']==era];by={v['mass']:v for v in a};fig,ax=plt.subplots(figsize=(9,5.8));y=[by[m]['p'] if m in by else np.nan for m in grid];ax.plot(grid,y,'o-',ms=4,color=colors[eras.index(era)],label='SR+CR, solo fondo')
 if a:ax.errorbar([v['mass'] for v in a],[v['p'] for v in a],yerr=np.array([[v['p']-v['CP95'][0] for v in a],[v['CP95'][1]-v['p'] for v in a]]),fmt='none',alpha=.35,color=colors[eras.index(era)],label='Intervallo MC 95%')
 ax.axhline(.05,color='firebrick',ls='--',label='Soglia p = 0,05');ax.set_yscale('log');ax.set_ylim(.005,1.1);ax.set_xlim(5.95,9.05);ax.set_xlabel('Ipotesi di massa [GeV]');ax.set_ylabel('GoF p-value');ax.set_title(f'Era {era}: GoF saturated — {len(a)}/31 punti\n100 toy per punto; r = 0 nella generazione e nel fit');ax.text(.02,.04,'SR 6–9,093; CR 9,093–10,507 GeV\nStesso test B-only: differenze tra masse da campioni toy distinti',transform=ax.transAxes,fontsize=9);ax.legend(loc='upper right' if era=='2023BPix' else 'center left');ax.grid(alpha=.2);fig.tight_layout();fig.savefig(R/f'gof_{era}.png',dpi=160);plt.close(fig)
fig,ax=plt.subplots(figsize=(10,6))
for era,color in zip(eras,colors):
 by={v['mass']:v for v in rows if v['era']==era};ax.plot(grid,[by[m]['p'] if m in by else np.nan for m in grid],'o-',ms=4,color=color,label=era)
ax.axhline(.05,color='firebrick',ls='--',label='p = 0,05');ax.set_yscale('log');ax.set_ylim(.005,1.1);ax.set_xlabel('Ipotesi di massa [GeV]');ax.set_ylabel('GoF p-value');ax.set_title(f'GoF per era — SR envelope Cheb6/Bern6 + CR Cheb4\n{len(rows)}/124 punti; 100 toy/punto, r fissato a zero');ax.legend(loc='lower left',ncol=5,fontsize=10);ax.grid(alpha=.2);fig.tight_layout();fig.savefig(R/'gof_all_eras.png',dpi=160);plt.close(fig)
pooled={}
for era in eras:
 a=[v for v in rows if v['era']==era];n=sum(v['N'] for v in a);k=sum(v['exceedances'] for v in a)
 if n:pooled[era]=dict(points=len(a),toys=n,p=k/n,exceedances=k,observed_min=min(v['observed'] for v in a),observed_max=max(v['observed'] for v in a),point_p_min=min(v['p'] for v in a),point_p_max=max(v['p'] for v in a),CP95=[float(beta.ppf(.025,k,n-k+1)) if k else 0,float(beta.ppf(.975,k+1,n-k)) if k<n else 1])
(R/'pooled.json').write_text(json.dumps(pooled,indent=2));complete=len(rows)==124
report=['# GoF per era: confronto con la figura 4.8','',f'Esecuzioni: {len(rows)}/124; 100 toy per punto; r=0 sia nella generazione sia nella valutazione del GoF.','', '**Precisazione sul confronto:** la tesi, sezione 4.2.3 e figura 4.8 (pagina stampata109), specifica 100 toy generati dal modello di solo fondo per era/massa, statistica saturated e soglia p>0.05. Non esplicita se r sia fisso o profilato durante il GoF. Questa campagna testa esplicitamente il solo fondo e non costituisce una replica verificata dello script originale, che non è stato fornito.','', 'Il modello dell’analisi attuale resta invariato: fit simultaneo SR6–9.093333 e CR9.093333–10.506667 GeV; envelope SR Cheb6/Bern6 positivo, CR Cheb4; medie e rese Upsilon condivise/libere, sigma/code fisse MC. Ogni era è testata separatamente. Nella tesi le tre regioni sono invece0.3–2.4,2–7,6–12 GeV. Non sono state ricostruite queste altre regioni né modificato il modello per copiarne i valori.','', 'Dentro ogni era il test B-only è identico alle31masse: stesso dataset e background, segnale annullato. I campioni da100toy usano seed distinti per ciascun punto. La dispersione della curva è Monte Carlo, non una sensibilità locale del GoF a una risonanza. Il grafico usa una scala logaritmica e la soglia0.05 come la figura4.8, ma evita estendere l’asse a1e-9:100toy risolvono frequenze di1/100. Barre per era: intervalli binomiali Clopper-Pearson95%.','', 'Il precedente p=0.47 riguarda il test combinato delle quattro ere e non è la media dei p-value per era.','', '|Era|Punti|Toy|P aggregato dei campioni|Intervallo MC95%|Range p nei singoli punti|','|---|---|---|---|---|---|']
for era,a in pooled.items():report.append(f"|{era}|{a['points']}|{a['toys']}|{a['p']:.4f}|{a['CP95'][0]:.4f}–{a['CP95'][1]:.4f}|{a['point_p_min']:.2f}–{a['point_p_max']:.2f}|")
report+=['','Il risultato aggregato combina i campioni toy dello stesso test ripetuto per era, non combina le ere. Va usato soltanto dopo la verifica della stabilità della statistica osservata e del modello B-only attraverso le masse. Non valida bias o coverage; nessun multistart sistematico dei toy. Tutti comandi, log, ROOT, seed e ricevute sono conservati.']
(R/'REPORT.md').write_text('\n'.join(report)+'\n')
html='<meta charset="utf-8"><style>body{font:18px system-ui;max-width:1200px;margin:35px auto}img{width:100%}</style><h1>GoF per era: riferimento figura4.8</h1>'+f'<p>{len(rows)}/124 punti completati. r=0 sia generazione sia fit;100toy per punto. Stesso test di fondo a tutte le masse: differenze Monte Carlo. Il comando originale della tesi non è disponibile: il trattamento di r nel suo GoF resta da verificare.</p><p><a href="REPORT.md">Metodo e risultati</a> · <a href="summary.csv">CSV</a></p><img src="gof_all_eras.png">'
for era in eras:html+=f'<h2>{era}</h2><img src="gof_{era}.png">'
(R/'index.html').write_text(html);(R/'completion.json').write_text(json.dumps(dict(executions_complete=complete,points=len(rows),valid=sum(v['valid'] for v in rows),requested_points=124,toys=sum(v['N'] for v in rows),exact_thesis_r_treatment_verified=False,background_only_test=True),indent=2));print('points',len(rows),'valid',sum(v['valid'] for v in rows),'pooled',pooled)
