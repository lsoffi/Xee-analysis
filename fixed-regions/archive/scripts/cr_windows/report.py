from pathlib import Path
import json,re,math,csv
import numpy as np,uproot
from scipy.stats import beta,norm
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
P=Path.cwd();R=P/'runs/20261029_cr_windows';cc=json.loads((R/'control_comparison.json').read_text());tags=list(cc);labels=['8,8–10,8','9,0–10,6','9,093–10,507','9,2–10,4'];eras=['2022','2022EE','2023','2023BPix'];out={}
windows=json.loads((P/'runs/20261029_signal_widths/common95.json').read_text())['windows'];signal=json.loads((P/'runs/20261029_signal_widths/parameters.json').read_text())
def read(d):
 with uproot.open(next(d.glob('higgsCombine*.root'))) as f:a=f['limit'].arrays(library='np')
 if d==R/'baseline/degree2/gof/toys/fit2':
  values=[a['limit']]
  for start in [-.05,.05]:
   with uproot.open(next((d.parent/f'start{start}').glob('higgsCombine*.root'))) as f:b=f['limit'].arrays(library='np')
   assert np.array_equal(a['iToy'],b['iToy']);values.append(b['limit'])
  best=np.min(values,axis=0);(R/'multistart_audit.json').write_text(json.dumps(dict(improved=int(sum(a['limit']-best>.01)),maximum_improvement=float(max(a['limit']-best))),indent=2));a['limit']=best
 assert np.all(np.isfinite(a['limit']));return a

def prob(a,obs):
 k=int(sum(a>=obs));n=len(a);return dict(p=k/n,k=k,n=n,interval95=[float(beta.ppf(.025,k,n-k+1)) if k else 0,float(beta.ppf(.975,k+1,n-k)) if k<n else 1],upper95=1-.05**(1/n) if k==0 else None)
for tag,label in zip(tags,labels):
 q=cc[tag]['3'];lo,hi=q['range'];d=R/tag/'degree3';v=read(d/'tail_profile')['SR_upsilon_events'];tail=dict(best=float(v[0]),low=float(v[1]),high=float(v[2]),error_minus=float(v[0]-v[1]),error_plus=float(v[2]-v[0]));y2={};log=(d/'Y2_profile/combine.log').read_text()
 for era in eras:
  name='CMS_EXO25020_bkgScale_upsilon2s_'+era;m=re.search(r'^\s*'+name+r'\s*:\s*([+-]?[\d.]+)\s+(-[\d.]+)/\+([\d.]+)',log,re.M);assert m,(tag,era);best,low,high=map(float,m.groups());y2[era]=dict(best=best,low=max(0,best+low),high=best+high)
 override=d/'Y2_2022_scan/interval.json'
 if override.exists():y2['2022']=json.loads(override.read_text())
 obs={};gof={}
 for n in [2,3]:
  root=R/tag/f'degree{n}/gof';obs[n]=float(read(root/'observed')['limit'][0]);a=read(root/f'toys/fit{n}')['limit'];assert len(a)==500;gof[str(n)]=dict(deviance=obs[n],**prob(a,obs[n]))
 a=read(R/tag/'degree2/gof/toys/fit2');b=read(R/tag/'degree2/gof/toys/fit3');assert np.array_equal(a['iToy'],b['iToy']);delta=a['limit']-b['limit'];assert min(delta)>-.001;test=dict(delta_deviance=obs[2]-obs[3],extra_parameters=4,**prob(delta,obs[2]-obs[3]))
 disjoint=[]
 for w in windows:
  low=6+math.floor((w['low95_common_GeV']-6)*75+1e-8)/75;high=6+math.ceil((w['high95_common_GeV']-6)*75-1e-8)/75
  if high<=lo+1e-8:disjoint.append(w['mass_GeV'])
 frac={}
 for s in signal:
  if 6.5<=s['mass']<=10:
   value=norm.cdf((hi-s['mean'])/s['sigmaL'])-norm.cdf((lo-s['mean'])/s['sigmaL']);frac[str(s['mass'])]=max(frac.get(str(s['mass']),0),float(value))
 pilot=P/'runs/20261029_upsilon_control_cheb3/results.json' if tag=='baseline' else d/'results.json';pilot=json.loads(pilot.read_text())['joint'];out[tag]=dict(label=label,range=[lo,hi],recovered_GeV=2-(hi-lo),tail_prediction=tail,Y2_profile=y2,gof=gof,degree_test=test,largest_existing_disjoint_mass=max(disjoint),signal_CR_fraction_max_era=frac,pilot_7p5=pilot,normalizations=q)
(R/'summary.json').write_text(json.dumps(out,indent=2));print(json.dumps({k:{a:b for a,b in v.items() if a not in ['normalizations','signal_CR_fraction_max_era']} for k,v in out.items()},indent=2))
fig,ax=plt.subplots(figsize=(10,5));xx=np.arange(4);ys=np.array([out[t]['tail_prediction']['best'] for t in tags]);err=np.array([[out[t]['tail_prediction'][key] for t in tags] for key in ['error_minus','error_plus']]);ax.errorbar(xx,ys,yerr=err,fmt='o',capsize=5,color='#0072B2');ax.axhline(ys[0],color='gray',ls='--');ax.set_xticks(xx,labels);ax.set_xlabel('Control region [GeV]');ax.set_ylabel('Eventi Υ(1S)+Υ(2S) previsti nella SR');ax.set_title('Trasferimento CR → SR 6,373–7,920 GeV, quattro ere\nChebyshev grado 3; intervalli di profilo ΔNLL=0,5');ax.grid(alpha=.2);fig.tight_layout();fig.savefig(R/'tail_prediction.png',dpi=160);fig.savefig(R/'tail_prediction.pdf');plt.close(fig)
fig,axs=plt.subplots(2,2,figsize=(12,8),sharex=True)
for i,(ax,era) in enumerate(zip(axs.flat,eras)):
 for proc,color,offset in [('upsilon1s','#0072B2',-.07),('upsilon2s','#D55E00',.07)]:
  if proc=='upsilon1s':
   yy=[cc[t]['3']['eras'][i][proc]['scale'] for t in tags];ee=[cc[t]['3']['eras'][i][proc]['error'] for t in tags]
  else:
   yy=[out[t]['Y2_profile'][era]['best'] for t in tags];ee=[[out[t]['Y2_profile'][era]['best']-out[t]['Y2_profile'][era]['low'] for t in tags],[out[t]['Y2_profile'][era]['high']-out[t]['Y2_profile'][era]['best'] for t in tags]]
  ax.errorbar(xx+offset,yy,yerr=ee,fmt='o',capsize=3,label=proc,color=color)
 ax.set_title(era);ax.set_ylabel('Scala rispetto al nominale');ax.set_xticks(xx,labels,rotation=15);ax.grid(alpha=.15);ax.legend()
fig.suptitle('Normalizzazioni dal solo fit nella CR\nΥ(1S): Hesse; Υ(2S): profilo ΔNLL=0,5 con bordo fisico');fig.tight_layout(rect=[0,0,1,.93]);fig.savefig(R/'normalizations.png',dpi=160);fig.savefig(R/'normalizations.pdf');plt.close(fig)
fig,ax=plt.subplots(figsize=(9,4.5))
for n,col,offset in [('2','#0072B2',-.06),('3','#D55E00',.06)]:
 p=np.array([out[t]['gof'][n]['p'] for t in tags]);ax.errorbar(xx+offset,p,yerr=np.sqrt(p*(1-p)/500),fmt='o',capsize=4,color=col,label=f'Grado {n}')
ax.axhline(.05,color='gray',ls='--',label='0,05');ax.set_ylim(0,1);ax.set_xticks(xx,labels);ax.set_ylabel('p-value GoF nella CR');ax.set_xlabel('Control region [GeV]');ax.set_title('500 toy per modello e finestra — barre: errore Monte Carlo');ax.legend();ax.grid(alpha=.15);fig.tight_layout();fig.savefig(R/'gof.png',dpi=160);fig.savefig(R/'gof.pdf');plt.close(fig)
body='<meta charset="utf-8"><style>body{font:17px system-ui;max-width:1300px;margin:30px auto}img{width:100%}td,th{padding:10px;border-bottom:1px solid #ddd}</style><h1>Riduzione della control region Υ</h1><p>Quattro ere; binning originale e bordi estesi ai bin. FormaΥoriginale; normalizzazioni libere. CR grado3, confronto diagnostico grado2. SRpilota invariata6,373–7,920, segnale gaussiano a7,5. GoF/likelihood ratio con500toy frequentist per generatore; 4000dataset complessivi. Profilo della somma delle code costruito come riparametrizzazione esatta dei yield, non somma degli errori indipendenti. Intervalli ΔNLL=0,5, non validazione della copertura ai bordi.</p><table><tr><th>CR</th><th>Recupero GeV</th><th>Code Υ previste in SR</th><th>pGoF3</th><th>p confronto2→3</th><th>Ultima massa con finestra esistente disgiunta</th></tr>'
for t in tags:
 v=out[t];a=v['tail_prediction'];p=v['degree_test'];pt=f'{p["p"]:.3f}' if p['k'] else f'<{p["upper95"]:.3f} (limite95%)';body+=f'<tr><td>{v["label"]}</td><td>{v["recovered_GeV"]:.3f}</td><td>{a["best"]:.1f} -{a["error_minus"]:.1f}/+{a["error_plus"]:.1f}</td><td>{v["gof"]["3"]["p"]:.3f}</td><td>{pt}</td><td>{v["largest_existing_disjoint_mass"]:.1f}</td></tr>'
body+='</table><p>Le finestre di ricerca esistenti sono quelle comuni al95% della DCB condizionata a6–12: la loro disgiunzione è un requisito distinto dalla contaminazione della CR da segnale. La frazione gaussiana inCR per ciascuna massa/era è in summary.json; non estendere automaticamente i risultati a un segnaleDCB. Nessuna correlazione fra stime ottenute con CR sovrapposte è stata calcolata: differenze dei valori centrali non sono pull indipendenti.</p>'
for name in ['tail_prediction','normalizations','gof']:body+=f'<img src="{name}.png">'
for tag in tags:
 link='../20261029_upsilon_control_cheb3/index.html' if tag=='baseline' else f'{tag}/degree3/index.html';body+=f'<p><a href="{link}">{out[tag]["label"]}: fit e componenti, pilot7,5, limiti/p-value</a></p>'
body+='<p><a href="summary.json">Risultati completi</a> · <a href="control_comparison.json">Normalizzazioni, correlazioni e diagnostica</a></p>';(R/'index.html').write_text(body)
