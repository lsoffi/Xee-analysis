from pathlib import Path
import json,numpy as np,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
P=Path.cwd();R=P/'runs/20261030_sr_bernstein7_bias';T=[json.loads(l) for l in (R/'final_toys.jsonl').read_text().splitlines()];A=[json.loads(l) for l in (R/'stable_asimov.jsonl').read_text().splitlines()];lookup={(x['degree'],x['mass'],x['toy']):x for x in T};out=[]
for d in range(3,8):
 for mass in np.round(np.arange(6.5,8.51,.1),1):
  tt=[x for x in T if x['degree']==d and abs(x['mass']-mass)<1e-6];pairs=[(x,lookup[(7,x['mass'],x['toy'])]) for x in tt if (7,x['mass'],x['toy']) in lookup];good=[(a,b) for a,b in pairs if a['valid'] and b['valid']];sig=next(a['error'] for a in A if a['degree']==d and abs(a['mass']-mass)<1e-6 and a['injection']=='zero');vals=np.array([(a['rhat']-b['rhat'])/sig for a,b in good]);rec=dict(degree=d,mass=float(mass),pairs_total=len(pairs),pairs_valid=len(good))
  if len(vals)>1:rec.update(mean_difference_sigma=float(vals.mean()),MCerror_sigma=float(vals.std(ddof=1)/np.sqrt(len(vals))))
  out.append(rec)
(R/'paired_summary.json').write_text(json.dumps(out,indent=2));col={3:'#9467bd',4:'#0072B2',5:'#E69F00',6:'#009E73',7:'#D55E00'};fig,axes=plt.subplots(1,2,figsize=(13,4.8),sharex=True)
for j,ax in enumerate(axes):
 for d in range(3,8):
  a=[x for x in out if x['degree']==d and 'mean_difference_sigma' in x];ax.errorbar([x['mass'] for x in a],[x['mean_difference_sigma'] for x in a],yerr=[x['MCerror_sigma'] for x in a],fmt='o-',ms=3,color=col[d],label=f'Chebyshev {d} − Chebyshev 7')
 ax.axhspan(-.2,.2,color='grey',alpha=.2);ax.axhline(0,color='k',lw=.6);ax.set_xlabel('Massa [GeV]');ax.set_ylabel('Media differenza r / σ Asimov del modello');ax.grid(alpha=.2)
 if j:ax.set_ylim(-.8,.8);ax.set_title('Zoom')
 else:ax.set_title('Tutti i gradi');ax.legend(fontsize=8)
fig.suptitle('Stessi toy Bernstein 7: confronto appaiato con il fit di chiusura\nMisura la differenza fra modelli; non sostituisce il bias assoluto o la copertura');fig.tight_layout();fig.savefig(R/'paired_bias.png',dpi=160);plt.close(fig)
print('paired counts',min(x['pairs_total'] for x in out),max(x['pairs_total'] for x in out))
