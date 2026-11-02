from run import *
import hashlib
from scipy.stats import beta
F=json.loads((R/'fits.json').read_text());rows=[json.loads(l) for l in (R/'pilot.jsonl').read_text().splitlines()];O=json.loads((R/'summary.json').read_text());out=[]
def cp(k,n):return [0 if k==0 else float(beta.ppf(.025,k,n-k+1)),1 if k==n else float(beta.ppf(.975,k+1,n-k))]
def truth(g):
 return [sum(np.array(I['eras'][ei]['templates'])[:,(edges[:-1]>=a-1e-8)&(edges[1:]<=b+1e-8)].sum(axis=1)@np.array(g['fits'][ei]['p'][-2:]) for ei in range(4)) for a,b in [(6,9.093333333333334),(10.506666666666668,12)]]
for wi in [0,5,10]:
 for kind in ['gof3','gof4','bias4','bias5','asimov4','asimov5']:
  for deg in ([4] if kind=='gof4' else [3,4]):
   a=[r for r in rows if r['idx']==2*wi+deg-3 and r['kind']==kind];toys=[]
   for k in sorted(set(x['toy'] for x in a)):
    fs=[x for x in a if x['toy']==k]
    if len(fs)!=4:continue
    toys.append(dict(toy=k,valid=all(f['valid'] for f in fs),mu_boundary=any(any(f['boundary'][2:]) for f in fs),D=sum(f['D'] for f in fs),prediction=np.sum([f['prediction'] for f in fs],axis=0).tolist(),error=np.sqrt(np.sum(np.array([f['prediction_error'] for f in fs])**2,axis=0)).tolist()))
   good=[t for t in toys if t['valid']];obs=O[2*wi+deg-3];rec=dict(window_id=wi,lo=obs['lo'],hi=obs['hi'],degree=deg,kind=kind,total=len(toys),valid=len(good),mu_boundary=sum(t['mu_boundary'] for t in toys),toys=toys)
   if kind==f'gof{deg}':
    n=len(toys);k=sum(t['D']>=obs['D'] for t in good);rec.update(observed_D=obs['D'],exceedances_valid=k,p_conditional=k/len(good) if good else None,p_possible_counting_failures=[k/n,(k+n-len(good))/n] if n else None,CP95_valid=cp(k,len(good)) if good else None)
   if kind.startswith(('bias','asimov')) and good:
    g=json.loads(((old/'bernstein_generator.json') if kind.endswith('4') else (R/'generator5.json')).read_text());tr=truth(g);p=np.array([t['prediction'] for t in good]);err=np.mean([t['error'] for t in good],axis=0);rec.update(truth=tr,conditional_bias_sigma=((p.mean(axis=0)-tr)/err).tolist(),conditional_coverage68=np.mean(abs(p-tr)<=np.array([t['error'] for t in good]),axis=0).tolist(),warning='Conditioned on successful fits; Wald intervals; not validated coverage')
   out.append(rec)
LR=[]
for wi in [0,5,10]:
 a=next(x for x in out if x['window_id']==wi and x['kind']=='gof3' and x['degree']==3);b=next(x for x in out if x['window_id']==wi and x['kind']=='gof3' and x['degree']==4);pairs=[(x,y) for x,y in zip(a['toys'],b['toys']) if x['valid'] and y['valid']];delta=O[2*wi]['D']-O[2*wi+1]['D'];ds=[x['D']-y['D'] for x,y in pairs];n=len(a['toys']);k=sum(d>=delta for d in ds);LR.append(dict(window_id=wi,delta_observed=delta,total=n,valid=len(ds),negative_deltas=sum(d<-.001 for d in ds),exceedances_valid=k,p_conditional=k/len(ds) if ds else None,p_possible_counting_failures=[k/n,(k+n-len(ds))/n],CP95_valid=cp(k,len(ds)) if ds else None))
(R/'pilot_summary.json').write_text(json.dumps(dict(configurations=out,LR=LR),indent=2));assert hashlib.sha256(Path(I['source']).read_bytes()).hexdigest()==I['sha256']
for f in F+rows:
 for j in range(2):assert abs(f['p'][-3+2*j]-f['nominal_shapes'][j][1])<1e-12
lines=['# Chebyshev grado 4 nella CR, soli μ delle Υ liberi','','## Configurazione','16 finestre, gradi3 e4 rifittati con procedura uniforme (128fit di era). Quattro ere indipendenti, μ1/μ2 nominali ±0,2GeV; σ e code originali MC fisse, normalizzazioni libere [0,10] volte nominale6–12. Accettanze ricalcolate. Input originale immutato, baselineSR+CR invariata.','Minuit/iminuit strategy2, tre seed, integrali esatti nei bin originali. Positività del continuo su griglia durante il fit e agli estremi/punti stazionari dopo. Si conserva il minimo tra i tentativi, anche se non supera i controlli: niente sostituzione silenziosa con un minimo peggiore.','Solo fitCR: non nuova significanza o limite del segnale. Le rese delle code sono confrontate sempre negli intervalli fissi6–9,0933 e10,5067–12. Errori Hesse solo diagnostici, non intervalli profilati.','', '## Pilot dei toy','50 pseudoesperimenti per finestra/generatore in tre finestre: attuale,9,2–10,4,9,293–10,307. GoF generato dal proprio fit per gradi3 e4. LR3→4 calibrato sugli stessi toy generati dal grado3. Fit3 e4 appaiati anche nei test delle rese.','Bernstein4 con shapeΥoriginali: medesimo generatore del test precedente, ma perCheb4 è un test di chiusura nella stessa famiglia polinomiale. Aggiunto Bernstein5positivo fitted alla CRlarga come scenario alternativo (GoF asintotico≈0,72, non calibrato contoy). Nessuno dei generatori è assunto verità dei dati. Asimov di entrambi in tutte3finestre e2gradi.','Fallimenti sempre conteggiati. P condizionati ai soli successi sono diagnostici: riportato anche intervallo[k/N,(k+fallimenti)/N] assegnando tutti i fallimenti ai due estremi. Questo NON è un intervallo di confidenza; non include l’errore Monte Carlo. CP95 separato nelJSON. Zero superamenti NON significa p=0.','', '|CR|Grado|Ere valide/4|Ere con μ al bordo|Devianza|Υ previste sotto CR ± Hesse|','|---|---:|---:|---:|---:|---|']
for wi in [0,5,10]:
 for deg in [3,4]:
  o=O[wi*2+deg-3];lines.append(f"|{o['lo']:.3f}–{o['hi']:.3f}|{deg}|{o['valid_eras']}|{o['mu_boundary_eras']}|{o['D']:.2f}|{o['prediction'][0]:.1f} ± {o['error'][0]:.1f}|")
lines+=['','## Test di grado','|CR|ΔD osservata|Toy LR validi|Superamenti validi|Intervallo contando fallimenti|','|---|---:|---:|---:|---|']
for q in LR:
 o=O[q['window_id']*2];lines.append(f"|{o['lo']:.3f}–{o['hi']:.3f}|{q['delta_observed']:.3f}|{q['valid']}/{q['total']}|{q['exceedances_valid']}|{q['p_possible_counting_failures']}|")
lines+=['','## Qualità dei toy e bias diagnostico','|CR|Generatore/test|Grado|Toy validi|Con μ al bordo|Bias/σ sotto / sopra (solo successi)|','|---|---|---:|---:|---:|---|']
for q in out:
 if q['kind']=='gof3' and q['degree']==4:continue
 lines.append(f"|{q['lo']:.3f}–{q['hi']:.3f}|{q['kind']}|{q['degree']}|{q['valid']}/{q['total']}|{q['mu_boundary']}|{q.get('conditional_bias_sigma','—')}|")
(R/'REPORT.md').write_text('\n'.join(lines));print('LR',LR)
for q in out:
 if q['kind'] in ['bias4','bias5','gof4']:print(q['window_id'],q['degree'],q['kind'],q['valid'],q['total'],q.get('conditional_bias_sigma'),q.get('p_possible_counting_failures'))
