from pathlib import Path
import json, html
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import chi2
P=Path.cwd();R=P/'runs/20261030_upper_sr_cr';A=json.loads((R/'results.json').read_text());eras=['2022','2022EE','2023','2023BPix']
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
def errors(n):
    alpha=1-.682689492137;low=np.zeros_like(n);s=n>0;low[s]=chi2.ppf(alpha/2,2*n[s])/2;high=chi2.ppf(1-alpha/2,2*(n+1))/2
    return np.array([n-low,high-n])
def axes_grid():
    fig=plt.figure(figsize=(13,9));outer=fig.add_gridspec(2,2,hspace=.35,wspace=.2);out=[]
    for cell in outer:
        gs=cell.subgridspec(2,1,height_ratios=[3,1],hspace=.06);ax=fig.add_subplot(gs[0]);rx=fig.add_subplot(gs[1],sharex=ax);out.append((ax,rx))
    return fig,out
def style(ax,rx,era,c):
    ax.set_title(era,fontweight='bold');ax.set_ylabel('Eventi / 13,33 MeV');ax.tick_params(labelbottom=False);ax.grid(alpha=.15)
    rx.axhline(0,color='gray',lw=.7);rx.set_ylabel('(D-M)/√M');rx.set_xlabel(r'$m_{ee}$ [GeV]');rx.set_xlim(*c['range']);rx.grid(alpha=.15)
def save(fig,name):
    fig.subplots_adjust(top=.88,bottom=.075)
    for ext in ['png','svg']:fig.savefig(R/(name+'.'+ext),dpi=160)
    plt.close(fig)
for a in A:
    for region in ['SR','CR']:
        fig,axes=axes_grid()
        for era,(ax,rx) in zip(eras,axes):
            c=next(c for c in a['channels'] if c['era']==era and c['region']==region)
            x=np.array([p['mass'] for p in c['points']]);n=np.array([p['observed'] for p in c['points']]);mu=np.array([p['expected'] for p in c['points']]);edges=np.linspace(*c['range'],c['bins']+1)
            ax.errorbar(x,n,yerr=errors(n),fmt='.',ms=3,color='black',elinewidth=.6,label='Dati')
            ax.stairs(mu,edges,color='#D55E00',lw=1.6,label='Fondo totale')
            for proc,col,label in [('dy','#0072B2',f"Chebyshev {a['degree'] if region=='SR' else 4}"),('upsilon1s','#E69F00','Υ(1S)'),('upsilon2s','#CC79A7','Υ(2S)')]:
                ax.stairs([p['components'][proc] for p in c['points']],edges,color=col,lw=1.1,label=label)
            ax.legend(fontsize=8,ncol=2);ax.set_ylim(bottom=0);rx.plot(x,(n-mu)/np.sqrt(mu),'.',color='#D55E00',ms=3);style(ax,rx,era,c)
        title='SR alta: 10,5067–12 GeV' if region=='SR' else 'CR: 9,0933–10,5067 GeV'
        fig.suptitle(f"{title} | Fit simultaneo solo fondo\nSR Chebyshev {a['degree']} + CR Chebyshev 4; Υ(1S), Υ(2S) condivise",fontsize=15,y=.98)
        note=f"FitDiagnostics: status {a['status']}, covQual {a['covQual']} | Binning originale | μ e rese Υ liberi, σ e code MC fisse"
        if a['boundaries']:note+=' | ATTENZIONE: 2 medie Υ(2S) al bordo'
        fig.text(.5,.018,note,ha='center',fontsize=9);save(fig,f"fit_degree{a['degree']}_{region}")
fig,axes=axes_grid();cols=['#999999','#0072B2','#D55E00']
for era,(ax,rx) in zip(eras,axes):
    for a,col in zip(A,cols):
        c=next(c for c in a['channels'] if c['era']==era and c['region']=='SR');x=np.array([p['mass'] for p in c['points']]);n=np.array([p['observed'] for p in c['points']]);mu=np.array([p['expected'] for p in c['points']]);edges=np.linspace(*c['range'],c['bins']+1)
        if a['degree']==1:ax.errorbar(x,n,yerr=errors(n),fmt='.',ms=3,color='black',elinewidth=.5,label='Dati')
        ax.stairs(mu,edges,color=col,lw=1.4,label=f"Totale: grado SR {a['degree']}")
        rx.plot(x,(n-mu)/np.sqrt(mu),color=col,lw=.6,alpha=.8)
    ax.legend(fontsize=8,ncol=2);ax.set_ylim(bottom=0);style(ax,rx,era,c)
fig.suptitle('SR alta: confronto dei tre fit simultanei con la stessa CR\nCurve del fondo totale; gradi diversi solo nel continuo della SR',fontsize=15,y=.98);save(fig,'comparison_SR')
summary=[]
for a in A:
    summary.append(dict(degree=a['degree'],status=a['status'],covQual=a['covQual'],D_SR=sum(c['deviance'] for c in a['channels'] if c['region']=='SR'),D_CR=sum(c['deviance'] for c in a['channels'] if c['region']=='CR'),D_total=a['D_GoF'],boundaries=a['boundaries']))
ms=json.loads((R/'multistart.json').read_text());assert all(a['exit_code']==0 and abs(a['delta_deviance'])<.01 for a in ms)
config=dict(status='exploratory_background_only_completed',SR_range=[10.506666666666668,12],CR_range=[9.093333333333334,10.506666666666668],SR_degrees=[1,2,3],CR_degree=4,eras=eras,signal_strength=0,lower_SR_included=False,events=52105,bins=872,gof_toys=0,bias_test=False,signal_scan=False,source='original m10.8 workspace, signal disabled in B-only fit',baseline_changed=False)
(R/'configuration.json').write_text(json.dumps(config,indent=2));(R/'summary.json').write_text(json.dumps(summary,indent=2))
trs=''.join(f"<tr><td>{a['degree']}</td><td>{a['D_SR']:.3f}</td><td>{a['D_CR']:.3f}</td><td>{a['D_total']:.3f}</td><td>{a['status']}/{a['covQual']}</td><td>{len(a['boundaries'])}</td></tr>" for a in summary)
gallery=''.join(f'<h2>Grado {a["degree"]}</h2><p><a href="degree{a["degree"]}/parameters.csv">Parametri pre/post fit (CSV)</a></p><img src="fit_degree{a["degree"]}_SR.png"><img src="fit_degree{a["degree"]}_CR.png">' for a in A)
body=f'''<!doctype html><html lang="it"><meta charset="utf-8"><title>SR alta + CR: primi fit</title><style>body{{font:17px system-ui;max-width:1200px;margin:35px auto;padding:0 20px;color:#12314b}}img{{width:100%}}table{{border-collapse:collapse}}td,th{{border:1px solid #ccd;padding:10px}}.note{{background:#fff3df;padding:15px}}</style>
<h1>SR alta + CR: confronto del fondo</h1><p>SR alta 10,506667–12 GeV, CR fissa 9,093333–10,506667 GeV. Quattro ere, 872 bin originali e 52.105 eventi (4.067 SR, 48.038 CR). La SR bassa non entra in questa prova. Nessuna modifica alla configurazione precedente.</p>
<p>Fit simultaneo Combine FitDiagnostics solo fondo (r=0). Continuo e sua normalizzazione indipendenti per regione/era. Υ(1S), Υ(2S) esplicite in entrambe; medie e normalizzazioni libere, condivise SR/CR per era. Sigma e code MC fisse. CR Chebyshev4 invariata. I fit precedenti forniscono soltanto inizializzazione, senza prior dai medesimi dati. Binning originale e nessun rebin.</p>
<table><tr><th>Grado SR</th><th>D SR</th><th>D CR</th><th>D totale</th><th>status/cov</th><th>Parametri al bordo</th></tr>{trs}</table>
<p>1→2: ΔD={summary[0]['D_total']-summary[1]['D_total']:.2f}; 2→3: ΔD={summary[1]['D_total']-summary[2]['D_total']:.2f}. Ogni incremento aggiunge quattro coefficienti, uno per era. Il grado 2 migliora nettamente il lineare; il grado 3 dà un miglioramento minore. Nel lineare le medie Υ(2S) 2022EE e 2023 sono al bordo; non accade per gradi 2/3.</p>
<p>Controlli: 3 fit status0/cov3, continui positivi agli estremi e punti stazionari, accettanze delle Υ post-fit confrontate con integrali ROOT, identità delle somme delle componenti, input originali immutati. Altre due inizializzazioni per ciascun grado: 5/6 status0/cov3; una del grado 2 ha status3/cov2, pur arrivando alla stessa NLL del fit nominale (ΔD=0,000234). Fit nominale e altra inizializzazione del grado 2 entrambi status0/cov3. Massimo |ΔD|={max(abs(a['delta_deviance']) for a in ms):.3g}. La devianza dal fit è coerente con Combine GoodnessOfFit osservato entro {max(abs(a['D_FD']-a['D_GoF']) for a in A):.3g}.</p>
<p class="note">D è la devianza saturata osservata, non un p-value GoF. Non sono stati generati toy né eseguiti test di bias, limiti o significanze. Queste prove motivano il grado 2 come primo modello da verificare; non fissano definitivamente il grado e non validano una ricerca di segnale. Il segnale presente nel workspace è disattivato (r=0): 10,8 nel nome dei comandi è solo compatibilità dell'input, non un'ipotesi testata.</p>
<img src="comparison_SR.png">{gallery}<p>Risultati completi: <a href="results.json">JSON</a>; <a href="multistart.json">controllo inizializzazioni</a>. Script in scripts/upper_sr_cr; comandi e log nelle directory dei singoli gradi.</p></html>'''
(R/'index.html').write_text(body)
print(json.dumps(summary,indent=2))
