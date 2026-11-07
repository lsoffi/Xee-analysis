from pathlib import Path
import json,csv,re
import numpy as np,uproot
from scipy.stats import norm
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
P=Path('/Users/imb/Documents/Codex/Analisi-statistica-Xee');R=P/'runs/20261029_upsilon_control_bern2';q=json.loads((R/'parameters.json').read_text());audit=json.loads((R/'joint/audit.json').read_text());eras=['2022','2022EE','2023','2023BPix']
for mode,region,stage,filename,title in [('control','','b','fitDiagnostics_control.root','Solo regione di controllo: Bernstein 2 + Υ(1S) + Υ(2S)'),('joint','SR_','s','fitDiagnostics_local.root','Fit simultaneo SR + CR: regione di ricerca a 7,5 GeV'),('joint','CR_','s','fitDiagnostics_local.root','Fit simultaneo SR + CR: regione di controllo')]:
 fig=plt.figure(figsize=(14,11));outer=fig.add_gridspec(2,2,hspace=.3,wspace=.2)
 with uproot.open(R/mode/'fitdiagnostics'/filename) as root:
  for i,(era,cell) in enumerate(zip(eras,outer),1):
   dr=root[f'shapes_fit_{stage}/{region}ch{i}'];gs=cell.subgridspec(2,1,height_ratios=[3,1],hspace=.07);ax=fig.add_subplot(gs[0]);rx=fig.add_subplot(gs[1],sharex=ax);v,edges=dr['total'].to_numpy();bw=np.diff(edges);x=(edges[:-1]+edges[1:])/2;data=dr['data'].values('y')*bw;pred=v*bw
   expected=audit['counts'][('CR_' if mode=='control' else region)+f'ch{i}'];assert abs(data.sum()-expected)<.1,(mode,region,era,data.sum(),expected)
   parts={k:dr[k].values()*bw for k in ['dy','upsilon1s','upsilon2s']};parts['Zd']=dr['Zd'].values()*bw if 'Zd' in dr else np.zeros_like(pred);assert np.max(abs(sum(parts.values())-pred)/np.maximum(pred,1e-30))<3e-7
   ax.errorbar(x,data,yerr=np.sqrt(data),fmt='.',ms=2.5,lw=.6,color='black',label='Dati');ax.stairs(pred,edges,color='#c53030',label='Totale')
   if mode=='joint':ax.stairs(root[f'shapes_fit_b/{region}ch{i}/total'].values()*bw,edges,color='gray',ls='--',label='B-only')
   for k,c,l in [('dy','#0072B2','Bernstein 2'),('upsilon1s','#E69F00','Υ(1S)'),('upsilon2s','#CC79A7','Υ(2S)'),('Zd','#009E73','Gaussiana segnale')]:
    if k!='Zd' or region=='SR_':ax.stairs(parts[k],edges,color=c,label=l)
   ax.set_title(era);ax.set_ylabel('Eventi / 13,33 MeV');ax.legend(fontsize=7,ncol=2);ax.tick_params(labelbottom=False);ax.grid(alpha=.15);rx.axhline(0,color='gray',lw=.7);rx.plot(x,(data-pred)/np.sqrt(pred),'.',color='#c53030',ms=3);rx.set_ylabel('(D−M)/√M');rx.set_xlabel('Massa [GeV]');rx.grid(alpha=.15);rx.set_xlim(edges[0],edges[-1])
 fig.suptitle(title+'\nNormalizzazioni Υ floating; shape originali; binning originale',fontsize=13);fig.subplots_adjust(top=.90,bottom=.07);name='control_fit' if mode=='control' else 'joint_'+region.rstrip('_')+'_fit'
 for ext in ['png','pdf']:fig.savefig(R/f'{name}.{ext}',dpi=145)
 plt.close(fig)
base=json.loads((P/'runs/20261029_upsilon_control/results.json').read_text())['joint'];out={'baseline':base,'joint':{}}
with uproot.open(R/'joint/fitdiagnostics/fitDiagnostics_local.root') as root:
 t=root['tree_fit_sb'];out['joint']['r']=float(t['r'].array(library='np')[0]);delta=float(t['nll_nll0'].array(library='np')[0]);out['joint']['Z_from_FD']=float(np.sqrt(max(0,-2*delta))) if out['joint']['r']>0 else 0
with uproot.open(next((R/'joint/significance').glob('higgsCombine*.root'))) as root:out['joint']['Z']=float(root['limit']['limit'].array(library='np')[0]);out['joint']['p0']=float(norm.sf(out['joint']['Z']))
with uproot.open(next((R/'joint/limits_stability').glob('higgsCombine*.root'))) as root:
 t=root['limit'];v=t['limit'].array(library='np');quant=t['quantileExpected'].array(library='np');assert len(v)==6;out['joint']['limits']={k:float(v[np.argmin(abs(quant-p))]) for k,p in [('observed',-1),('expected_2p5',.025),('expected_16',.16),('expected_50',.5),('expected_84',.84),('expected_97p5',.975)]}
obs=out['joint']['limits']['observed'];log=(R/'joint/limits_stability/combine.log').read_text();pairs=re.findall(r'At r =\s*([\d.eE+-]+):[^\n]*CLs\s*=\s*([\d.eE+-]+)',log);rr,cls=min(pairs,key=lambda x:abs(float(x[0])-obs));out['joint']['CLs_crossing']=dict(r=float(rr),CLs=float(cls),verified=abs(float(cls)-.05)<.005 and abs(float(rr)/obs-1)<.02);assert out['joint']['CLs_crossing']['verified'];assert abs(out['joint']['Z']-out['joint']['Z_from_FD'])<.1
params=json.loads((P/'runs/20261029_signal_widths/parameters.json').read_text());out['CR_signal_fraction']={r['era']:float((norm.sf((8.8-r['mean'])/r['sigmaL'])-norm.sf((10.8-r['mean'])/r['sigmaL']))/(norm.cdf((12-r['mean'])/r['sigmaL'])-norm.cdf((6-r['mean'])/r['sigmaL']))) for r in params if r['mass']==7.5}
out['quality']={k:{st:{a:b for a,b in v.items() if a in ['status','covQual','edm']} for st,v in ss.items()} for k,ss in q.items()};(R/'results.json').write_text(json.dumps(out,indent=2));print(json.dumps({k:v for k,v in out.items() if k!='baseline'},indent=2))
with (R/'normalizations.csv').open('w') as f:
 w=csv.writer(f);w.writerow(['era','scale_1s','error_1s','scale_2s','error_2s','correlation_1s_2s'])
 for a in q['control']['b']['normalizations']:w.writerow([a['era'],a['upsilon1s']['scale'],a['upsilon1s']['error'],a['upsilon2s']['scale'],a['upsilon2s']['error'],a['correlation_1s_2s']])
fig,ax=plt.subplots(figsize=(8,5))
for i,record in enumerate([base,out['joint']]):
 lim=record['limits'];ax.plot([lim['expected_2p5'],lim['expected_97p5']],[i,i],color='#F6D743',lw=18,solid_capstyle='butt');ax.plot([lim['expected_16'],lim['expected_84']],[i,i],color='#69BE46',lw=18,solid_capstyle='butt');ax.plot(lim['expected_50'],i,'k|',ms=22,mew=2);ax.plot(lim['observed'],i,'o',color='#0072B2',label='Osservato' if i==0 else None)
ax.set_yticks([0,1],['Grado 2 + controllo Υ','Bernstein 2 + controllo Υ']);ax.set_xlabel('Limite 95% CLs su r');ax.set_title('Segnale gaussiano a 7,5 GeV — Chebyshev 2 / Bernstein 2');ax.grid(axis='x',alpha=.2);ax.legend();fig.tight_layout();fig.savefig(R/'limits_comparison.png',dpi=170);fig.savefig(R/'limits_comparison.pdf');plt.close(fig)
body='<!doctype html><meta charset="utf-8"><title>Controllo Υ</title><style>body{font:17px system-ui;max-width:1400px;margin:30px auto}img{width:100%}</style><h1>Misura delle Υ e fit simultaneo a 7,5 GeV</h1><p>CR8,8–10,8; SR6,373–7,920 GeV. Continui indipendenti, normalizzazioni Υ condivise, nessun doppio conteggio. Nessuna prior percentuale arbitraria. Studio preliminare: il modello del continuo in CR resta da validare.</p>'
body+='<p><a href="normalizations.csv">Normalizzazioni CSV</a> · <a href="parameters.json">Parametri e covarianze</a> · <a href="results.json">Limiti e p-value</a></p>'
for name in ['control_fit','joint_SR_fit','joint_CR_fit','limits_comparison']:body+=f'<h2>{name}</h2><img src="{name}.png">'
(R/'index.html').write_text(body)
