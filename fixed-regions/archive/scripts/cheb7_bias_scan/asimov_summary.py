from pathlib import Path
import json,re,csv
import numpy as np,uproot
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
R=Path('runs/20261030_cheb7_bias_scan');rows=[]
for study in sorted(R.glob('m*/study.json')):
 a=json.loads(study.read_text())
 for scenario in ['zero','injected']:
  truth=0 if scenario=='zero' else a['injection']
  for model in ['fixed','sliding']:
   d=study.parent/scenario/'asimov'/model;p=d/'fitDiagnostics_asimov.root'
   if not p.exists() or not (d/'combine.log').exists():continue
   try:
    f=uproot.open(p);s=f['tree_fit_sb'].arrays(library='np');b=f['tree_fit_b'].arrays(library='np');rh=float(s['r'][0]);err=float(s['rErr'][0]);dn=float(s['nll_nll0'][0]);matches=re.findall(r'Fit (B-only|S\+B), status =\s*(-?\d+), numBadNLL =\s*(\d+), covariance quality =\s*(-?\d+)',(d/'combine.log').read_text());cov={st:int(v) for st,_,_,v in matches};rows.append(dict(mass=a['mass'],scenario=scenario,model=model,injected=truth,r=rh,bias=rh-truth,error=err,Z=float(np.sqrt(max(0,-2*dn))) if rh>0 else 0,delta_nll=dn,b_status=int(b['fit_status'][0]),s_status=int(s['fit_status'][0]),b_cov=cov.get('B-only',-1),s_cov=cov.get('S+B',-1)))
   except (OSError,KeyError,IndexError):continue
(R/'asimov_summary.json').write_text(json.dumps(rows,indent=2))
if rows:
 with (R/'asimov_summary.csv').open('w') as f:w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
fig,axes=plt.subplots(1,2,figsize=(12,4))
for model,color in [('fixed','tab:red'),('sliding','tab:blue')]:
 for scenario,ls in [('zero','-'),('injected','--')]:
  a=[r for r in rows if r['model']==model and r['scenario']==scenario];axes[0].plot([r['mass'] for r in a],[r['bias'] for r in a],ls+'o',ms=3,color=color,label=model+' '+scenario)
 a=[r for r in rows if r['model']==model and r['scenario']=='zero'];axes[1].plot([r['mass'] for r in a],[r['Z'] for r in a],'o-',ms=3,color=color,label=model)
axes[0].set_ylabel('r fitted - r iniettato');axes[1].set_ylabel('Z locale nominale senza segnale');axes[0].axhline(0,color='gray')
for ax in axes:ax.set_xlabel('Massa [GeV]');ax.grid(alpha=.2);ax.legend(fontsize=8)
fig.suptitle('Dataset Asimov da Chebyshev7 + Y: nessuna fluttuazione statistica');fig.tight_layout();fig.savefig(R/'asimov_bias.png',dpi=150);fig.savefig(R/'asimov_bias.pdf');plt.close(fig)
# Illustrate one mass, generator expectation and fitted components, both methods.
m='7.5';fig=plt.figure(figsize=(13,8));outer=fig.add_gridspec(2,2,hspace=.35,wspace=.2)
for j,method in enumerate(['fixed','sliding']):
 path=R/('m'+m)/'zero/asimov'/method/'fitDiagnostics_asimov.root'
 if not path.exists():continue
 f=uproot.open(path);data=pred=bg=sig=None
 for i in range(1,5):
  d=f[f'shapes_fit_s/SR_ch{i}'];v,e=d['total'].to_numpy();bw=np.diff(e);yy=d['data'].values('y')*bw;pp=v*bw;ss=d['Zd'].values()*bw;bb=pp-ss;data=yy if data is None else data+yy;pred=pp if pred is None else pred+pp;sig=ss if sig is None else sig+ss;bg=bb if bg is None else bg+bb
 x=(e[:-1]+e[1:])/2
 for zoom in [0,1]:
  gs=outer[zoom,j].subgridspec(2,1,height_ratios=[3,1],hspace=.06);ax=fig.add_subplot(gs[0]);rx=fig.add_subplot(gs[1],sharex=ax);ax.plot(x,data,'k.',ms=2,label='Generatore7 + Y, nessun segnale');ax.plot(x,pred,color='tab:red',label='Fit totale3 + Y + DCB');ax.plot(x,bg,color='tab:blue',ls='--',label='Fondo fitted');ax.plot(x,sig,color='tab:green',label='Segnale spurio fitted');ax.set_title(method+(' — zoom' if zoom else ' — intera SR'));ax.legend(fontsize=7);ax.set_ylabel('Eventi / bin');ax.tick_params(labelbottom=False);rx.axhline(0,color='gray');rx.plot(x,(data-pred)/np.sqrt(np.maximum(pred,1e-12)),'.',ms=2);rx.set_ylabel('(A-F)/sqrt(F)',fontsize=7);rx.set_xlabel('Massa [GeV]');rx.grid(alpha=.2)
  if zoom:
   low=max(7,e[0]);high=min(8.2,e[-1]);rx.set_xlim(low,high);mask=(x>=low)&(x<=high);ax.set_ylim(min(0,float(np.min(sig[mask]))),float(max(data[mask]))*1.2)
fig.suptitle('m=7.5 GeV: fit allo stesso fondo atteso, senza fluttuazioni e senza segnale\nQuattro ere sommate; la CR è presente in entrambi i fit');fig.subplots_adjust(top=.88,bottom=.07);fig.savefig(R/'asimov_fit_m7.5.png',dpi=150);fig.savefig(R/'asimov_fit_m7.5.pdf');plt.close(fig)
print('Asimov fits',len(rows))
