from pathlib import Path
import json,re,csv,time
import numpy as np,uproot
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
P=Path.cwd();R=P/'runs/20261030_cheb7_bias_scan';toyrows=[];summary=[];paired=[];eras=['2022','2022EE','2023','2023BPix']
def writecsv(name,rows):
 if not rows:return
 with (R/name).open('w') as f:w=csv.DictWriter(f,fieldnames=list(dict.fromkeys(k for r in rows for k in r)));w.writeheader();w.writerows(rows)
for study in sorted(R.glob('m*/study.json')):
 a=json.loads(study.read_text());d=study.parent
 for label in ['zero','injected']:
  truth=0 if label=='zero' else a['injection'];cache={}
  for model in ['fixed','sliding']:
   dd=d/label/model;receipt=dd/'combine.log.receipt.json';path=dd/'fitDiagnostics_toys.root'
   if not receipt.exists() or not path.exists():continue
   with uproot.open(path) as f:
    b=f['tree_fit_b'].arrays(library='np');s=f['tree_fit_sb'].arrays(library='np')
   matches=re.findall(r'Fit (B-only|S\+B), status =\s*(-?\d+), numBadNLL =\s*(\d+), covariance quality =\s*(-?\d+)',(dd/'combine.log').read_text());bc=[int(v) for st,_,_,v in matches if st=='B-only'];sc=[int(v) for st,_,_,v in matches if st=='S+B'];assert len(bc)==len(sc)==len(s['r'])==30,(dd,len(bc),len(sc),len(s['r']));rows=[]
   for i in range(30):
    rh=float(s['r'][i]);err=float(s['rErr'][i]);dn=float(s['nll_nll0'][i]);valid=int(b['fit_status'][i])==0 and int(s['fit_status'][i])==0 and np.isfinite(rh+dn) and abs(rh)<99.9 and dn<=.01;pullvalid=valid and bc[i]==3 and sc[i]==3 and np.isfinite(err) and err>0
    r=dict(mass=a['mass'],scenario=label,injection=truth,model=model,toy=i+1,r=rh,error=err,b_status=int(b['fit_status'][i]),s_status=int(s['fit_status'][i]),b_cov=bc[i],s_cov=sc[i],delta_nll=dn,Z=float(np.sqrt(max(0,-2*dn))) if rh>0 else 0,valid_point=bool(valid),valid_error=bool(pullvalid),pull=(rh-truth)/err if pullvalid else None);rows.append(r);toyrows.append(r)
   cache[model]=rows;v=[r for r in rows if r['valid_point']];ve=[r for r in rows if r['valid_error']];vals=np.array([r['r'] for r in v]);errors=np.array([r['error'] for r in ve]);n=len(v);ss=dict(mass=a['mass'],scenario=label,injection=truth,model=model,n_total=30,n_valid=n,n_valid_errors=len(ve),mean_r=float(vals.mean()) if n else None,bias=float(vals.mean()-truth) if n else None,bias_mc_error=float(vals.std(ddof=1)/np.sqrt(n)) if n>1 else None,sd_r=float(vals.std(ddof=1)) if n>1 else None,median_error=float(np.median(errors)) if len(errors) else None,mean_pull=float(np.mean([r['pull'] for r in ve])) if ve else None,sd_pull=float(np.std([r['pull'] for r in ve],ddof=1)) if len(ve)>1 else None,n_Zge3=sum(r['Z']>=3 for r in v),fraction_Zge3=sum(r['Z']>=3 for r in v)/n if n else None);summary.append(ss)
  if len(cache)==2:
   v=[(f,s) for f,s in zip(cache['fixed'],cache['sliding']) if f['valid_point'] and s['valid_point']];delta=np.array([f['r']-s['r'] for f,s in v]);paired.append(dict(mass=a['mass'],scenario=label,n_paired=len(v),mean_delta=float(delta.mean()) if len(v) else None,se_delta=float(delta.std(ddof=1)/np.sqrt(len(v))) if len(v)>1 else None))
writecsv('toys.csv',toyrows);writecsv('summary.csv',summary);writecsv('paired.csv',paired);(R/'summary.json').write_text(json.dumps(summary,indent=2));(R/'paired.json').write_text(json.dumps(paired,indent=2));status=dict(completed_fit_batches=len(summary),expected_fit_batches=84,toy_fits=len(toyrows),expected_toy_fits=2520,valid_points=sum(r['valid_point'] for r in toyrows),valid_errors=sum(r['valid_error'] for r in toyrows));(R/'progress.json').write_text(json.dumps(status,indent=2))
if not summary:print(status);raise SystemExit
fig,axes=plt.subplots(2,2,figsize=(13,8))
for i,label in enumerate(['zero','injected']):
 for model,color,offset in [('fixed','tab:red',-.012),('sliding','tab:blue',.012)]:
  a=sorted([s for s in summary if s['scenario']==label and s['model']==model and s['n_valid']>1],key=lambda s:s['mass']);x=np.array([s['mass'] for s in a]);axes[i,0].errorbar(x+offset,[s['bias'] for s in a],yerr=[s['bias_mc_error'] for s in a],fmt='o-',ms=3,lw=1,capsize=2,color=color,label=model);a=[s for s in a if s['median_error']];axes[i,1].errorbar(np.array([s['mass'] for s in a])+offset,[s['bias']/s['median_error'] for s in a],yerr=[s['bias_mc_error']/s['median_error'] for s in a],fmt='o-',ms=3,lw=1,capsize=2,color=color,label=model)
 for j,ax in enumerate(axes[i]):ax.axhline(0,color='gray');ax.set_xlabel('Massa [GeV]');ax.set_ylabel('Bias in r' if j==0 else 'Bias / mediana errore Hesse');ax.set_title('Nessun segnale iniettato' if label=='zero' else 'Segnale iniettato: r = limite atteso di riferimento');ax.legend();ax.grid(alpha=.2)
fig.suptitle('Generatore Chebyshev7 + Y; fit Chebyshev3 + Y + DCB\n30 toy/configurazione; barre = errore Monte Carlo sulla media; punti validi soltanto');fig.tight_layout();fig.savefig(R/'bias.png',dpi=150);fig.savefig(R/'bias.pdf');plt.close(fig)
fig,axes=plt.subplots(1,2,figsize=(12,4))
for model,color in [('fixed','tab:red'),('sliding','tab:blue')]:
 a=sorted([s for s in summary if s['scenario']=='zero' and s['model']==model and s['n_valid']],key=lambda s:s['mass']);axes[0].plot([s['mass'] for s in a],[s['fraction_Zge3'] for s in a],'o-',color=color,label=model)
 for label,ls in [('zero','-'),('injected','--')]:
  a=sorted([s for s in summary if s['scenario']==label and s['model']==model],key=lambda s:s['mass']);axes[1].plot([s['mass'] for s in a],[30-s['n_valid'] for s in a],'o'+ls,color=color,label=model+' '+label)
axes[0].set_ylabel('Frazione toy con Z locale >= 3');axes[0].set_title('Generazione senza segnale (solo fit validi)');axes[0].set_ylim(-.02,1.02);axes[1].set_ylabel('Fit esclusi / 30');axes[1].set_title('Problemi numerici o minimi incoerenti')
for ax in axes:ax.set_xlabel('Massa [GeV]');ax.grid(alpha=.2);ax.legend(fontsize=8)
fig.tight_layout();fig.savefig(R/'diagnostics.png',dpi=150);fig.savefig(R/'diagnostics.pdf');plt.close(fig)
html='<meta charset="utf-8"><meta http-equiv="refresh" content="60"><style>body{font:17px system-ui;max-width:1400px;margin:25px auto}img{width:100%}td,th{padding:7px;border-bottom:1px solid #ddd}</style><h1>Bias scan: generatore Chebyshev7</h1>'+f'<p>{status["completed_fit_batches"]}/84 batch; {status["toy_fits"]}/2520 fit. Quattro ere, CR9.0933–10.5067, SRfissa6–9.0933 oppure finestra95DCB. Generatore grado7 in SR/CR, fit grado3 in SR/CR. Stessi toy e stesse misure ausiliarie per i due fit. 30toy/config, errori Hesse, studio preliminare. Esclusioni numeriche esplicitate; non validazione della copertura.</p><a href="summary.csv">CSV</a> · <a href="toys.csv">Tutti i toy</a> · <a href="REPORT.md">Report</a><img src="bias.png"><img src="diagnostics.png"><table><tr><th>Massa</th><th>Scenario</th><th>Modello</th><th>Validi</th><th>Bias ± errore MC</th><th>Z>=3</th></tr>'
for s in summary:html+=f'<tr><td>{s["mass"]}</td><td>{s["scenario"]}</td><td>{s["model"]}</td><td>{s["n_valid"]}/30</td><td>{s["bias"]} ± {s["bias_mc_error"]}</td><td>{s["n_Zge3"]}</td></tr>'
(R/'index.html').write_text(html+'</table>');print(status)
