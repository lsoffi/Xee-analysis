from pathlib import Path
import json, html
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import chi2
P=Path.cwd();R=P/'runs/20261030_upper_sr_cr_mu_fixed';A=json.loads((R/'results.json').read_text());eras=['2022','2022EE','2023','2023BPix']
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
        fig.suptitle(f"{title} | Fit simultaneo solo fondo\nSR Chebyshev {a['degree']} + CR Chebyshev 4; Υ(1S), Υ(2S): μ fissate al fit CR-only",fontsize=15,y=.98)
        note=f"FitDiagnostics: status {a['status']}, covQual {a['covQual']} | Binning originale | μ Υ dal fit CR-only, fissate; rese libere; σ e code MC fisse"
        if a['boundaries']:note+=' | ATTENZIONE: parametri al bordo'
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

OLD=json.loads((P/'runs/20261030_upper_sr_cr/results.json').read_text())
CR=json.loads((R/'cr_only/results.json').read_text())
summary=[]
for a,old in zip(A,OLD):
    summary.append(dict(degree=a['degree'],status=a['status'],covQual=a['covQual'],D_fixed=a['D_GoF'],D_floating=old['D_GoF'],difference=a['D_GoF']-old['D_GoF'],boundaries=a['boundaries']))
(R/'summary.json').write_text(json.dumps(summary,indent=2))
rows=''.join(f"<tr><td>{a['degree']}</td><td>{a['D_fixed']:.3f}</td><td>{a['D_floating']:.3f}</td><td>{a['difference']:.3f}</td><td>{a['status']}/{a['covQual']}</td></tr>" for a in summary)
murows=''
for era in eras:
    m=[next(a for a in CR['means'] if a['name']==f'Ymean_{proc}_{era}') for proc in ['upsilon1s','upsilon2s']]
    murows+=f"<tr><td>{era}</td><td>{m[0]['value']:.5f} ± {m[0]['error']:.5f}</td><td>{m[1]['value']:.5f} ± {m[1]['error']:.5f}</td></tr>"
gallery=''.join(f'<h2>SR grado {a["degree"]}</h2><a href="degree{a["degree"]}/parameters.csv">Parametri liberi pre/post</a><img src="fit_degree{a["degree"]}_SR.png"><img src="fit_degree{a["degree"]}_CR.png">' for a in A)
body=f'''<!doctype html><html lang="it"><meta charset="utf-8"><title>SR alta + CR: medie fissate dalla CR</title><style>body{{font:17px system-ui;max-width:1200px;margin:35px auto;padding:0 20px;color:#12314b}}img{{width:100%}}table{{border-collapse:collapse}}td,th{{border:1px solid #ccd;padding:10px}}.note{{background:#fff3df;padding:15px}}</style>
<h1>SR alta + CR: μ determinate solo nella CR</h1>
<p>Correzione richiesta: (1) fit B-only CR 9,093333–10,506667 GeV, Chebyshev4, μ e rese di Υ1/2 libere per era; (2) μ fissate a quelle stime nel fit simultaneo con SR alta10,506667–12. Le rese Υ restano libere e condivise. Sigma/code MC fisse; continuo indipendente per regione/era e libero. Non si aggiunge alcuna prior dal fit preliminare agli stessi dati.</p>
<table><tr><th>Era</th><th>μ Υ(1S) [GeV]</th><th>μ Υ(2S) [GeV]</th></tr>{murows}</table>
<p>Errori Hesse del solo fit CR, diagnostici: non propagati nel simultaneo. Tutte le medie sono interne ai limiti. Fit CR scelto tra tre inizializzazioni confrontando la devianza senza offset, D={CR['selected']['deviance']:.6f}, status0/cov3. Due inizializzazioni valide; la terza ha status3/cov2. Parametri e tentativi in <a href="cr_only/results.json">cr_only/results.json</a>.</p>
<table><tr><th>Grado SR</th><th>D con μ fissate</th><th>D con μ libere (precedente)</th><th>ΔD</th><th>status/cov</th></tr>{rows}</table>
<p>I tre fit nominali corretti hanno status0/cov3. Verificate positività del continuo, integrali Υ post-fit rispetto a ROOT, 8 μ esattamente costanti con valori della CR, 8 normalizzazioni Υ libere, e input invariati. Bin e dati invariati rispetto al confronto precedente. Grado2: differenza di devianza soltanto {summary[1]['difference']:.3f} rispetto a μ libere. Nessuna selezione definitiva del grado.</p>
<p class="note">Risultati condizionati alle μ stimate. Congelarle non ne propaga l'incertezza. Per una futura calibrazione a toy occorre ripetere la fase CR-only e il congelamento in ciascun pseudo-esperimento; non usare queste stime come misure indipendenti della stessa CR. Nessun toy, bias, limite o significanza eseguito in questa prova. SR bassa non inclusa; i risultati storici con μ libere sono conservati.</p>{gallery}</html>'''
(R/'index.html').write_text(body)
print(json.dumps(summary,indent=2))
