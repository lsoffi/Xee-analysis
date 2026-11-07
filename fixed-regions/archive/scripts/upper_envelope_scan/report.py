from pathlib import Path
import json,math,time,sys,datetime
import numpy as np,uproot
from scipy.stats import beta,norm
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
P=Path.cwd();R=P/'runs/20261030_upper_envelope_scan';eras=['2022','2022EE','2023','2023BPix'];colors={1:'#D55E00',2:'#0072B2'}
def collect():
    a=[json.loads(p.read_text()) for p in sorted(R.glob('degree*/m*/results.json'))]
    b=[json.loads(p.read_text()) for p in sorted(R.glob('bias/degree*/m*/*/m*/result.json'))]
    g=[]
    for x in a:
        toys=[json.loads(p.read_text()) for p in (R/f'gof/degree{x["degree"]}').glob(f'toy*/m{x["mass"]:.1f}/result.json')];n=len(toys)
        if not n:continue
        k=sum(t['score']>x['gof_observed'] for t in toys);lo=0 if k==0 else beta.ppf(.025,k,n-k+1);hi=1 if k==n else beta.ppf(.975,k+1,n-k)
        failed=sum(not t['CR_valid'] for t in toys);valid_exceed=sum(t['CR_valid'] and t['score']>x['gof_observed'] for t in toys)
        g.append(dict(degree=x['degree'],mass=x['mass'],n=n,k=k,p=k/n,low=float(lo),high=float(hi),CR_invalid=failed,p_failure_low=valid_exceed/n,p_failure_high=(valid_exceed+failed)/n,r_boundaries=sum(t['r_boundary'] for t in toys),observed=x['gof_observed']))
    return a,b,g
def save(fig,name):
    fig.tight_layout();fig.savefig(R/(name+'.png'),dpi=160);fig.savefig(R/(name+'.svg'));plt.close(fig)
def plotfits(a):
    D=R/f'degree{a["degree"]}/m{a["mass"]:.1f}';src=D/'joint/fitdiagnostics/fitDiagnostics_upper.root'
    if (D/'fit_SR.png').exists() and (D/'fit_SR.png').stat().st_mtime>src.stat().st_mtime:return
    with uproot.open(src) as f:
        for region in ['SR','CR']:
            fig=plt.figure(figsize=(13,9));outer=fig.add_gridspec(2,2,hspace=.33,wspace=.2)
            for i,(era,cell) in enumerate(zip(eras,outer),1):
                dr=f[f'shapes_fit_s/{region}_ch{i}'];v,e=dr['total'].to_numpy();bw=np.diff(e);x=(e[:-1]+e[1:])/2;data=dr['data'].values('y')*bw;pred=v*bw
                gs=cell.subgridspec(2,1,height_ratios=[3,1],hspace=.06);ax=fig.add_subplot(gs[0]);rx=fig.add_subplot(gs[1],sharex=ax)
                ax.errorbar(x,data,yerr=np.sqrt(data),fmt='.',ms=2,color='black',lw=.5,label='Dati');ax.stairs(pred,e,color='#D55E00',lw=1.5,label='Totale S+B');ax.stairs(f[f'shapes_fit_b/{region}_ch{i}/total'].values()*bw,e,ls='--',color='gray',label='Solo fondo')
                for proc,col,lab in [('dy','#0072B2',f'Envelope grado {a["degree"]}' if region=='SR' else 'Chebyshev4'),('upsilon1s','#E69F00','Υ(1S)'),('upsilon2s','#CC79A7','Υ(2S)'),('Zd','#009E73','Segnale DCB')]:ax.stairs(dr[proc].values()*bw,e,color=col,label=lab)
                ax.set_title(era);ax.set_ylabel('Eventi / 13,33 MeV');ax.tick_params(labelbottom=False);ax.legend(fontsize=7,ncol=2);ax.grid(alpha=.15);rx.axhline(0,color='gray',lw=.7);rx.plot(x,(data-pred)/np.sqrt(np.maximum(pred,1e-10)),'.',color='#D55E00',ms=2.5);rx.set_ylabel('(D-M)/√M');rx.set_xlabel('Massa [GeV]');rx.set_xlim(e[0],e[-1]);rx.grid(alpha=.15)
            fig.suptitle(f'mX={a["mass"]:.1f} GeV | {region} | Envelope Cheb/Bern grado {a["degree"]}, CR4\nμ Υ dalla sola CR e fissate; segnale DCB in SR e CR'+(' | Resa estrapolata' if a['mass']>11 else ''),fontsize=13);fig.subplots_adjust(top=.89,bottom=.07)
            fig.text(.5,.015,f'r̂={a["r"]:.4g}; Z locale={a["Z"]:.3f}; FD B/S status={a["fit"]["b"]["status"]}/{a["fit"]["s"]["status"]}, cov={a["fit"]["b"]["covQual"]}/{a["fit"]["s"]["covQual"]}'+(' | FLAG NUMERICO' if not a['numerical_ok'] else ''),ha='center',fontsize=9)
            for ext in ['png','svg']:fig.savefig(D/f'fit_{region}.{ext}',dpi=140)
            plt.close(fig)
def report():
    a,b,g=collect();R.mkdir(exist_ok=True)
    for x in a:plotfits(x)
    for what in ['limits','pvalues']:
        fig,ax=plt.subplots(figsize=(11,5.6))
        for degree in [1,2]:
            rows=sorted([q for q in a if q['degree']==degree],key=lambda q:q['mass']);x=[q['mass'] for q in rows]
            if not rows:continue
            if what=='limits':
                ax.plot(x,[q['limits']['observed'] for q in rows],'o-',color=colors[degree],label=f'Grado {degree}: osservato');ax.plot(x,[q['limits']['expected_50'] for q in rows],'--',color=colors[degree],label=f'Grado {degree}: atteso mediano');ax.fill_between(x,[q['limits']['expected_2p5'] for q in rows],[q['limits']['expected_97p5'] for q in rows],color=colors[degree],alpha=.08);ax.fill_between(x,[q['limits']['expected_16'] for q in rows],[q['limits']['expected_84'] for q in rows],color=colors[degree],alpha=.16)
            else:ax.plot(x,[q['p0'] for q in rows],'o-',color=colors[degree],label=f'Envelope grado {degree}')
            bad=[q for q in rows if not q['numerical_ok']]
            if bad:ax.scatter([q['mass'] for q in bad],[q['limits']['observed'] if what=='limits' else q['p0'] for q in bad],marker='x',s=100,color='black',zorder=5,label=f'Flag numerici grado {degree}')
        ax.axvline(11,color='gray',ls=':');ax.axvspan(11,12,color='#fff0d0',alpha=.35);ax.set_yscale('log');ax.set_xlim(10.48,12.02);ax.set_xlabel('Massa ipotizzata [GeV]');ax.set_ylabel('Limite 95% CLs su r' if what=='limits' else 'p-value locale (eccesso)');ax.set_title(f'SR alta + CR | {len(a)}/32 configurazioni\n'+('Bande native 1σ/2σ PROVVISORIE: quantili superiori non validati' if what=='limits' else 'Eccessi nominali: nessuno dei due gradi supera uniformemente il bias Asimov'));ax.grid(alpha=.2);ax.legend(fontsize=8);save(fig,what)
    if b:
        fig,axes=plt.subplots(1,3,figsize=(15,4.5))
        for ax,label,title in zip(axes,['r0','rExpected','r10Expected'],['r vero = 0','r vero = limite atteso','r vero = 10 × limite atteso']):
            for degree in [1,2]:
                rows=sorted([q for q in b if q['degree']==degree and q['label']==label],key=lambda q:q['mass']);valid=[q for q in rows if q['valid']];bad=[q for q in rows if not q['valid'] and q['bias_sigma'] is not None]
                ax.plot([q['mass'] for q in rows],[q['bias_sigma'] if q['valid'] else np.nan for q in rows],'o-',color=colors[degree],label=f'Grado {degree}')
                ax.scatter([q['mass'] for q in bad],[q['bias_sigma'] for q in bad],color=colors[degree],marker='x',s=60)
            ax.axhspan(-.2,.2,color='gray',alpha=.15);ax.axhline(0,color='gray',lw=.6);ax.set_title(title);ax.set_xlabel('Massa [GeV]');ax.set_ylabel('(r̂-r vero) / σHesse(r)');ax.grid(alpha=.2);ax.legend();ax.set_xlim(10.48,12.02)
        fig.suptitle(f'Bias Asimov: generatore Chebyshev3, ristima μ nella CR | {len(b)}/96 fit\nCroci: fit/CR con flag; fascia ±0,2σ indicativa; errori Hesse condizionati alle μ fissate',fontsize=12);fig.tight_layout(rect=(0,0,1,.9));fig.savefig(R/'bias.png',dpi=160);fig.savefig(R/'bias.svg');plt.close(fig)
    if g:
        fig,ax=plt.subplots(figsize=(11,5.5))
        for degree in [1,2]:
            rows=sorted([q for q in g if q['degree']==degree],key=lambda q:q['mass']);x=np.array([q['mass'] for q in rows]);v=np.array([q['p'] for q in rows]);lo=np.array([q['low'] for q in rows]);hi=np.array([q['high'] for q in rows])
            if rows:
                nmin=min(q['n'] for q in rows);nmax=max(q['n'] for q in rows);label=f'{nmin}' if nmin==nmax else f'{nmin}–{nmax}'
                ax.errorbar(x+(-.008 if degree==1 else .008),v,yerr=[v-lo,hi-v],fmt='o-',color=colors[degree],capsize=2,label=f'Grado {degree}: {label} toy/punto')
        ax.axhline(.05,color='red',ls='--',label='0,05');ax.set_ylim(-.02,1.02);ax.set_xlim(10.48,12.02);ax.set_xlabel('Massa [GeV]');ax.set_ylabel('p GoF saturated, fit S+B con r libero');failed=max(q['CR_invalid'] for q in g);ax.set_title('GoF: B-only toys, ristima CR-only di μ prima di ciascun fit\nIntervalli MC 95%; toy correlati fra masse'+(f'; fino a {failed} toy/punto con flag CR' if failed else ''));ax.grid(alpha=.2);ax.legend(fontsize=9);save(fig,'gof')
    summary=dict(observed=a,bias=b,gof=g,counts=dict(observed=len(a),observed_ok=sum(q['numerical_ok'] for q in a),bias=len(b),bias_valid=sum(q['valid'] for q in b),gof_scores=sum(q['n'] for q in g),gof_target=3200))
    (R/'summary.json').write_text(json.dumps(summary,indent=2));(R/'progress.json').write_text(json.dumps(summary['counts'],indent=2))
    rows=''.join(f'<tr><td>{q["degree"]}</td><td>{q["mass"]:.1f}</td><td>{q["r"]:.4g}</td><td>{q["Z"]:.3f}</td><td>{q["p0"]:.3g}</td><td>{q["limits"]["observed"]:.4g}</td><td>{q["limits"]["expected_50"]:.4g}</td><td>{"OK" if q["numerical_ok"] else "FLAG"}</td><td><a href="degree{q["degree"]}/m{q["mass"]:.1f}/fit_SR.png">SR</a> / <a href="degree{q["degree"]}/m{q["mass"]:.1f}/fit_CR.png">CR</a></td></tr>' for q in a)
    body=f'''<!doctype html><html lang="it"><meta charset="utf-8"><meta http-equiv="refresh" content="60"><title>Scan SR alta + CR</title><style>body{{font:16px system-ui;max-width:1200px;margin:30px auto;color:#12314b}}img{{width:100%}}td,th{{padding:7px;border:1px solid #ccd}}table{{border-collapse:collapse}}.note{{padding:15px;background:#fff0d0}}</style><h1>Scan SR alta + CR: envelope grado 1 e 2</h1><p>Stato: {len(a)}/32 osservati; {len(b)}/96 bias Asimov; {summary['counts']['gof_scores']}/3200 score GoF.</p><p>SR 10,506667–12; CR 9,093333–10,506667 con Cheb4. Medie Υ stimate nella sola CR e poi fissate; normalizzazioni Υ condivise e libere, sigma/code MC fisse. Segnale DCB in entrambe le regioni, unico r. Due envelope separati, Cheb1/Bern1 e Cheb2/Bern2; indici indipendenti per era, 16 combinazioni. Stesso grado = stesso spazio polinomiale con domini dei coefficienti diversi, non famiglie indipendenti.</p><p class="note">Sopra 11: shape dalle formule parametriche originali; rese estrapolate con fit log-quadratico 10–11 ancorato a 11. Incertezze relative lnN mantenute a m11; variazione del modello di estrapolazione non inclusa come nuisance. Bias: generatore Cheb3, Asimov r0/rExpected/r10Expected, ripetendo la stima delle μ sul campione iniettato. GoF: 100 toy per grado, riutilizzati tra le 16 masse, con r libero nel test e μ ristimate per ciascun toy. Controlli e figure ancora in aggiornamento; i flag non sono risultati validati.</p><img src="limits.png"><img src="pvalues.png">'''
    for name in ['bias.png','gof.png','gof_distributions.png','signal_rate_extrapolation.png','extended_signal_shapes.png']:
        if (R/name).exists():body+=f'<img src="{name}">'
    completed=len(a)==32 and len(b)==96 and sum(q['n'] for q in g)==3200
    body=body.replace('Controlli e figure ancora in aggiornamento; i flag non sono risultati validati.', 'Calcoli richiesti completati, con limitazioni documentate sotto.' if completed else 'GoF ancora in aggiornamento; i flag non sono risultati validati.')
    if completed:body=body.replace('<meta http-equiv="refresh" content="60">','')
    body+=f'<h2>Lettura dei risultati e limiti della verifica</h2><p><b>Il bias Asimov non supporta la scelta definitiva di nessuno dei due gradi.</b> La soglia indicativa è |r̂−r vero|/σHesse &lt; 0,2. Il GoF qui usa un fit S+B con r libero: può dare un p-value accettabile anche quando la componente di segnale assorbe un errore del fondo. Non è quindi una verifica sostitutiva del bias. Le significanze sono locali e nominali, non corrette per lo scan.</p><p class="note"><b>Bande dei limiti provvisorie.</b> I limiti osservati attraversano CLs=0,05 entro la tolleranza verificata. Per molti quantili attesi superiori il valore coincide con il bordo numerico dell’intervallo di ricerca di Combine. Le bande native sono conservate e mostrate, ma non considerate quantitative finché il crossing dei quantili superiori non è validato. Un controllo con algoritmo alternativo ha prodotto NaN ed è stato escluso.</p><p>FitDiagnostics consente r negativo. A m=11,0 (grado 1) e m=12,0 (grado 2) esistono minimi distinti a r negativo e positivo: il minimo negativo del FitDiagnostics e il ramo fisico positivo usato dalla significanza danno risultati diversi. Un controllo con r≥0 e partenza positiva riproduce la Z nativa; la partenza dal bordo r=0 può invece restarvi bloccata. Le croci nere mantengono visibile questa dipendenza dal minimo. Anche gli errori Hesse possono variare fra inizializzazioni quasi equivalenti in likelihood.</p>'
    if (R/'validation.json').exists():
        audit=json.loads((R/'validation.json').read_text())
        br=''.join(f'<tr><td>{q["degree"]}</td><td>{q["injection"]}</td><td>{q["valid"]}/{q["n"]}</td><td>{q["within_0p2"]}/{q["valid"]}</td><td>{q["minimum"]:.3f} … {q["maximum"]:.3f}</td></tr>' for q in audit['bias_summary'])
        body+=f'<table><tr><th>Grado</th><th>Iniezione</th><th>Fit con qualità richiesta</th><th>Entro ±0,2σ fra i validi</th><th>Intervallo bias/σ</th></tr>{br}</table><p>Gli Asimov sono 96 campioni deterministici: non sono ensemble di pseudoesperimenti e non misurano copertura. I 26 fit con flag sono conservati, non trattati come bias validati. Il limite atteso usato per l’iniezione è quello disponibile alla generazione; il valore esatto di r vero è nel CSV.</p>'
    if g:
        gr=''.join(f'<tr><td>{q["degree"]}</td><td>{q["mass"]:.1f}</td><td>{q["k"]}/{q["n"]}</td><td>{q["p"]:.3f}</td><td>[{q["low"]:.3f}, {q["high"]:.3f}]</td><td>{q["CR_invalid"]}</td></tr>' for q in g)
        body+=f'<h2>Goodness of fit</h2><p>Statistica nativa Combine saturated; confronto stretto score(toy)&gt;score(dati). Intervalli binomiali Clopper–Pearson al 95%, sola precisione Monte Carlo. Nessun toy viene eliminato perché il suo score è sfavorevole. Ogni campione ripete il fit B-only della CR e il congelamento delle μ. Le normalizzazioni Υ e i continui sono poi profilati nel simultaneo. Il generatore è il B-only del rispettivo envelope, perciò questo GoF non testa tutte le possibili forme alternative.</p><p>I toy con CR non convergente restano nel denominatore e sono segnalati: i campi p_failure_low/high del JSON danno l’intervallo ottenuto trattando ciascuno di questi score come ignoto, da nessun superamento a tutti superamenti. È distinto dall’errore Monte Carlo. I punteggi GoF finiti non costituiscono una verifica multistart del minimo per tutti i 3200 fit.</p><table><tr><th>Grado</th><th>Massa</th><th>Superamenti/toy</th><th>p GoF</th><th>Intervallo MC 95%</th><th>CR con flag</th></tr>{gr}</table>'
    rates=json.loads((R/'signal_rate_extrapolation.json').read_text())
    rr=''.join(f'<tr><td>{q["era"]}</td><td>{q["N11"]:.5g}</td><td>{q["b"]:.6f}</td><td>{q["c"]:.6f}</td><td>{q["predictions"][-1]["rate"]:.5g}</td><td>{100*(q["predictions"][-1]["loglinear_alternative"]/q["predictions"][-1]["rate"]-1):+.1f}%</td></tr>' for q in rates)
    body+=f'<h2>Estensione del segnale sopra 11 GeV</h2><p>Per ciascuna era: N(m)=N(11) exp[b(m−11)+c(m−11)²], fit non pesato in log sulle rese nominali 10–11 GeV, con ancoraggio esatto a 11. Questa è la resa nominale a r=1 nel dominio originale 6–12. Le rese regionali sono N(m)×integraleDCB(regione)/integraleDCB(6–12). Non viene inventata una nuova sezione d’urto in unità fisiche.</p><table><tr><th>Era</th><th>N(11)</th><th>b [GeV⁻¹]</th><th>c [GeV⁻²]</th><th>N(12)</th><th>Alternativa log-lineare a 12</th></tr>{rr}</table><p>La scelta usa controlli lasciando fuori 10,6–11 o 10,9–11: RMS relativo sotto 1% per la log-quadratica. Questo NON è un errore sull’estrapolazione a 12. Le shape riusano le formule parametriche DCB originali per tutte le ere, con verifica sui punti già disponibili e nuisance ±1. Le incertezze relative lnN sono mantenute a 11; la differenza fra modelli di estrapolazione non è aggiunta come nuisance. I dati oltre 12 non sono disponibili, quindi la sensibilità al bordo è condizionata al modello e al dominio 6–12.</p><h2>Fit e valori per ciascuna ipotesi</h2><p>Plot SR/CR: dati, totale S+B, B-only, continuo, Υ(1S), Υ(2S), DCB e residui. Sono i fit con r libero anche negativo; non sempre coincidono con il ramo fisico della significanza. Parametri completi in degreeN/mM/parameters.csv e risultati ROOT nelle sottodirectory joint.</p>'
    body+=f'<table><tr><th>Grado</th><th>Massa</th><th>r̂</th><th>Z locale</th><th>p0</th><th>Limite obs</th><th>Atteso</th><th>Numerica</th><th>Fit</th></tr>{rows}</table><h2>Riproducibilità</h2><p>Gli input originali sono intatti. Script in scripts/upper_envelope_scan; workspace, comandi esatti, log e ricevute nelle sottodirectory. SR e CR disgiunte: 112+106 bin per era, 872 bin totali; 4067 eventi SR e 48038 CR. Gli errori delle μ della CR non sono propagati come nuisance: l’inferenza sui dati resta condizionata alle stime fissate, mentre bias e GoF ripetono il primo stadio. Il bias verifica tutta la procedura: con r iniettato diverso da zero il fit preliminare B-only della CR può assorbire contaminazione di segnale anche se il fit simultaneo successivo include correttamente la DCB in CR. Non attribuire quindi ogni bias al solo grado del continuo.</p><p><a href="summary.json">Risultati JSON</a> · <a href="validation.json">Controlli e limitazioni</a> · <a href="observed_results.csv">Fit/limiti/p0 CSV</a> · <a href="bias_results.csv">Bias CSV</a> · <a href="gof_toy_scores.csv">Score toy CSV</a> · <a href="signal_acceptances.csv">Rese e accettanze CSV</a> · <a href="signal_rate_extrapolation.json">Coefficienti e validazione estrapolazione</a></p><p>Aggiornato {datetime.datetime.now().astimezone().isoformat(timespec="seconds")}.</p></html>';(R/'index.html').write_text(body)
    print(summary['counts'],flush=True)
if __name__=='__main__':
    if '--watch' in sys.argv:
        for _ in range(240):
            try:report()
            except Exception as e:print(type(e).__name__,str(e),flush=True)
            time.sleep(30)
    else:report()
