from pathlib import Path
import json,csv,hashlib
import numpy as np,uproot
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
P=Path.cwd();R=P/'runs/20261030_quick_no_upsilon';cases=json.loads((R/'cases.json').read_text());rows=[];toyrows=[]
fig,axes=plt.subplots(2,3,figsize=(13,7))
for j,a in enumerate(cases):
 for i,label in enumerate(['zero','injected']):
  d=R/('m'+a['mass'])/label;truth=0 if label=='zero' else a['injection'];res={}
  for model,folder in [('noY','fit'),('reference','reference')]:
   ar=uproot.open(d/folder/'fitDiagnostics_toys.root')['tree_fit_sb'].arrays(library='np');good=(ar['fit_status']==0)&np.isfinite(ar['r'])&np.isfinite(ar['rErr'])&(ar['rErr']>0);res[model]=(ar,good)
   for n in range(len(good)):toyrows.append(dict(mass=a['mass'],injection=truth,model=model,toy=n+1,r=float(ar['r'][n]),error=float(ar['rErr'][n]),status=int(ar['fit_status'][n]),badnll=int(ar['numbadnll'][n])))
  ar,g=res['noY'];ref,gr=res['reference'];paired=g&gr;delta=ar['r'][paired]-ref['r'][paired];asim=uproot.open(d/'asimov/fit/fitDiagnostics_toys.root')['tree_fit_sb'].arrays(library='np');bias=ar['r'][g]-truth
  row=dict(mass=float(a['mass']),label=label,injected=truth,attempted=len(g),good_noY=int(g.sum()),good_reference=int(gr.sum()),mean_r=float(ar['r'][g].mean()),mean_bias=float(bias.mean()),mean_bias_se=float(bias.std(ddof=1)/np.sqrt(len(bias))),median_error=float(np.median(ar['rErr'][g])),mean_bias_reference=float((ref['r'][gr]-truth).mean()),paired_shift=float(delta.mean()),paired_shift_se=float(delta.std(ddof=1)/np.sqrt(len(delta))),asimov_r=float(asim['r'][0]),asimov_error=float(asim['rErr'][0]),asimov_status=int(asim['fit_status'][0]),asimov_bias=float(asim['r'][0]-truth));rows.append(row)
  ax=axes[i,j];lo=min(ar['r'].min(),ref['r'].min(),truth);hi=max(ar['r'].max(),ref['r'].max(),truth);bins=np.linspace(lo-.05*(hi-lo),hi+.05*(hi-lo),13);ax.hist(ref['r'][gr],bins=bins,histtype='step',color='tab:blue',label='Con Y fissate vere');ax.hist(ar['r'][g],bins=bins,histtype='step',color='tab:red',label='Senza Y');ax.axvline(truth,color='black',ls='--',label='r iniettato');ax.set_title(f"m={a['mass']} GeV; r iniettato={truth:.4f}");ax.set_xlabel('r fitted');ax.set_ylabel('Toy');ax.grid(alpha=.2)
axes[0,0].legend(fontsize=8);fig.suptitle('Test preliminare: generazione Chebyshev3 + Y1S + Y2S + DCB\n30 toy per configurazione; stesse fluttuazioni per i due modelli di fit');fig.tight_layout();fig.savefig(R/'toy_distributions.png',dpi=150);fig.savefig(R/'toy_distributions.pdf');plt.close(fig)
fig,axes=plt.subplots(1,2,figsize=(11,4))
for label,color in [('zero','tab:blue'),('injected','tab:orange')]:
 rr=[r for r in rows if r['label']==label];x=np.array([r['mass'] for r in rr])+(-.025 if label=='zero' else .025)
 axes[0].errorbar(x,[r['mean_bias'] for r in rr],yerr=[r['mean_bias_se'] for r in rr],fmt='o',capsize=4,color=color,label=label)
 axes[1].errorbar(x,[r['paired_shift'] for r in rr],yerr=[r['paired_shift_se'] for r in rr],fmt='o',capsize=4,color=color,label=label)
for ax,title in zip(axes,['Bias medio r_fit - r_iniettato','Spostamento medio senza Y - con Y (toy appaiati)']):ax.axhline(0,color='gray');ax.set_title(title,fontsize=10);ax.set_xlabel('Massa [GeV]');ax.set_ylabel('r');ax.legend();ax.grid(alpha=.2)
fig.suptitle('Barre: errore Monte Carlo sulla media, non incertezza del singolo fit');fig.tight_layout();fig.savefig(R/'bias.png',dpi=150);fig.savefig(R/'bias.pdf');plt.close(fig)
(R/'summary.json').write_text(json.dumps(rows,indent=2))
for name,data in [('summary.csv',rows),('toys.csv',toyrows)]:
 with (R/name).open('w') as f:w=csv.DictWriter(f,fieldnames=list(data[0]));w.writeheader();w.writerows(data)
print(json.dumps(rows,indent=2));(R/'validation.json').write_text(json.dumps(dict(source_unchanged=all(hashlib.sha256(Path(a['source']).read_bytes()).hexdigest()==a['source_sha256'] for a in cases),toys_generated=180,toy_fits=len(toyrows),fit_status_nonzero=sum(r['status']!=0 for r in toyrows),boundary_r=sum(abs(r['r'])>4.99 for r in toyrows)),indent=2))
