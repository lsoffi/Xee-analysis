from pathlib import Path
import json,numpy as np,matplotlib,csv
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import beta
P=Path.cwd();R=P/'runs/20261030_sr_bernstein7_bias';A=[json.loads(l) for l in (R/'stable_asimov.jsonl').read_text().splitlines()];T=[json.loads(l) for l in (R/'final_toys.jsonl').read_text().splitlines()] if (R/'final_toys.jsonl').exists() else [];G=json.loads((R/'generators.json').read_text());GT=[json.loads(l) for l in (R/'gof.jsonl').read_text().splitlines()];eras=['2022','2022EE','2023','2023BPix'];edges=np.array(json.loads((P/'runs/20261030_cr_optimization/input.json').read_text())['edges'])[:339];COL={3:'#9467bd',4:'#0072B2',5:'#E69F00',6:'#009E73',7:'#D55E00'};summary=[];toy_summary=[]
for d in range(3,8):
 aa=[x for x in A if x['degree']==d];good=[x for x in aa if x['valid']];maxbias=max(abs(x['bias_sigma']) for x in good);summary.append(dict(degree=d,asimov_valid=len(good),asimov_total=len(aa),max_asimov_abs_bias_sigma=maxbias,passes_asimov=bool(len(good)==42 and maxbias<.2)))
 for mass in np.round(np.arange(6.5,8.51,.1),1):
  tt=[x for x in T if x['degree']==d and abs(x['mass']-mass)<1e-6];gg=[x for x in tt if x['valid']];ref=next(x for x in aa if abs(x['mass']-mass)<1e-6 and x['injection']=='zero');vals=np.array([x['rhat'] for x in gg]);rec=dict(degree=d,mass=float(mass),total=len(tt),valid=len(gg),error_valid=sum(x.get('error_valid',True) for x in gg),sigma_asimov=ref['error'])
  if len(gg)>1:rec.update(mean_r=float(vals.mean()),bias_sigma=float(vals.mean()/ref['error']),MCerror_sigma=float(vals.std(ddof=1)/np.sqrt(len(vals))/ref['error']),pull_mean=float(np.mean([x['bias_sigma'] for x in gg if x.get('error_valid',True)])) if any(x.get('error_valid',True) for x in gg) else None,std_r=float(vals.std(ddof=1)),wald68_conditional=float(np.mean([abs(x['rhat'])<=x['error'] for x in gg if x.get('error_valid',True)])) if any(x.get('error_valid',True) for x in gg) else None)
  toy_summary.append(rec)
(R/'summary.json').write_text(json.dumps(summary,indent=2));(R/'toy_summary.json').write_text(json.dumps(toy_summary,indent=2))
with (R/'toy_summary.csv').open('w') as f:w=csv.DictWriter(f,fieldnames=sorted(set(k for r in toy_summary for k in r)));w.writeheader();w.writerows(toy_summary)
good=[]
for k in range(100):
 fs=[x for x in GT if x['toy']==k]
 if len(fs)==4 and all(x['valid'] for x in fs):good.append(sum(x['D'] for x in fs))
D=sum(x['D'] for x in G);kk=sum(x>=D for x in good);gof=dict(observed=D,total=100,valid=len(good),exceedances_valid=kk,p_conditional=kk/len(good),p_failure_bounds=[kk/100,(kk+100-len(good))/100]);(R/'gof_summary.json').write_text(json.dumps(gof,indent=2))
for region,sl in [('SR',slice(0,232)),('CR',slice(232,338))]:
 fig,axes=plt.subplots(2,2,figsize=(11,7));xx=((edges[:-1]+edges[1:])/2)[sl]
 for ei,ax in enumerate(axes.flat):
  f=next(g for g in G if g['era']==eras[ei]);n=np.array(f['counts'])[sl];ax.errorbar(xx,n,yerr=np.sqrt(n),fmt='k.',ms=2.5,label='Dati');ax.plot(xx,np.array(f['mu'])[sl],label='Totale');ax.plot(xx,np.array(f['continuum'])[sl],'--',label='Continuo');ax.plot(xx,np.array(f['Y'][0])[sl],label='Υ(1S)');ax.plot(xx,np.array(f['Y'][1])[sl],label='Υ(2S)');ax.set_title(eras[ei]);ax.set_xlabel('Massa [GeV]');ax.set_ylabel('Eventi / bin originale');ax.grid(alpha=.2)
 axes[0,0].legend(fontsize=8);fig.suptitle(f'Generatore fitted: Bernstein 7 nella SR + CR definitiva Chebyshev 4 — {region}');fig.tight_layout();fig.savefig(R/f'generator_{region}.png',dpi=140);plt.close(fig)
fig,axes=plt.subplots(2,2,figsize=(13,8),sharex=True)
for row,inj in enumerate(['zero','expected']):
 for col in range(2):
  ax=axes[row,col]
  for d in range(3,8):
   aa=sorted([x for x in A if x['degree']==d and x['injection']==inj],key=lambda x:x['mass']);ax.plot([x['mass'] for x in aa],[x['bias_sigma'] if x['valid'] else np.nan for x in aa],'o-',ms=3,color=COL[d],label=f'Chebyshev {d}');bad=[x for x in aa if not x['valid']];ax.plot([x['mass'] for x in bad],[x['bias_sigma'] for x in bad],'x',color=COL[d],ms=6)
  ax.axhspan(-.2,.2,color='grey',alpha=.2);ax.axhline(0,color='k',lw=.6);ax.set_title(('r vero = 0' if inj=='zero' else 'r vero = limite atteso di riferimento')+(' · zoom' if col else ''));ax.grid(alpha=.2);ax.set_xlabel('Massa ipotizzata [GeV]');ax.set_ylabel('(r fit − r vero) / σ del fit')
  if col:ax.set_ylim(-.8,.8)
axes[0,0].legend(fontsize=8,ncol=2);fig.suptitle('Generatore Bernstein 7 — fit simultanei SR+CR, Chebyshev 3–7\nAsimov: bias di modello senza fluttuazioni; croci = fit non valido; fascia grigia ±0,2σ');fig.tight_layout();fig.savefig(R/'asimov_bias.png',dpi=160);fig.savefig(R/'asimov_bias.pdf');plt.close(fig)
if T:
 fig,axes=plt.subplots(1,2,figsize=(13,4.8),sharex=True)
 for col,ax in enumerate(axes):
  for d in range(3,8):
   aa=[x for x in toy_summary if x['degree']==d and 'bias_sigma' in x];ax.errorbar([x['mass']+(d-5)*.002 for x in aa],[x['bias_sigma'] for x in aa],yerr=[x['MCerror_sigma'] for x in aa],fmt='o-',ms=3,lw=1,color=COL[d],label=f'Chebyshev {d}')
  ax.axhspan(-.2,.2,color='grey',alpha=.2);ax.axhline(0,color='k',lw=.6);ax.set_xlabel('Massa ipotizzata [GeV]');ax.set_ylabel('Media r fit / σ Asimov');ax.grid(alpha=.2);ax.set_title('Tutto il range del bias' if not col else 'Zoom')
  if col:ax.set_ylim(-1,1)
 axes[0].legend(fontsize=8,ncol=2);completed=min(x['total'] for x in toy_summary);fig.suptitle(f'Bernstein 7, r vero = 0 — {completed}/30 toy completati per configurazione\nBarre: errore Monte Carlo sulla media; minimi convergenti, covarianze separate nel report; punti correlati fra masse');fig.tight_layout();fig.savefig(R/'toy_bias.png',dpi=160);plt.close(fig)
html='<meta charset="utf-8"><meta http-equiv="refresh" content="60"><style>body{font:17px system-ui;max-width:1300px;margin:auto}img{max-width:100%}</style><h1>Bernstein 7 → Chebyshev 3–7, SR+CR</h1><p><a href="REPORT.md">Report</a> · <a href="toy_summary.csv">Risultati toy</a></p><img src="asimov_bias.png"><img src="toy_bias.png"><img src="paired_bias.png"><h2>Generatore fitted ai dati</h2><img src="generator_SR.png"><img src="generator_CR.png">';(R/'index.html').write_text(html);print('Asimov',summary);print('GOF',gof);print('Toy fit completed',len(T))
