from pathlib import Path
import json,re,csv
import numpy as np,uproot
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
P=Path.cwd();R=P/'runs/20261030_cheb7_fit7_comparison';OLD=P/'runs/20261030_cheb7_bias_scan';oldsum=json.loads((OLD/'summary.json').read_text());summary=[];toyrows=[];asimov=[];eras=['2022','2022EE','2023','2023BPix']
for d in sorted(R.glob('m*')):
 if not d.is_dir():continue
 mass=float(d.name[1:]);study=json.loads((OLD/d.name/'study.json').read_text())
 for label in ['zero','injected']:
  truth=0 if label=='zero' else study['injection']
  for kind in ['toys','asimov']:
   dd=d/label/kind;receipt=dd/'receipt.json';path=dd/('fitDiagnostics_toys.root' if kind=='toys' else 'fitDiagnostics_asimov.root')
   if not receipt.exists() or not path.exists():continue
   with uproot.open(path) as f:b=f['tree_fit_b'].arrays(library='np');s=f['tree_fit_sb'].arrays(library='np')
   matches=re.findall(r'Fit (B-only|S\+B), status =\s*(-?\d+), numBadNLL =\s*(\d+), covariance quality =\s*(-?\d+)',(dd/'combine.log').read_text());bc=[int(v) for st,_,_,v in matches if st=='B-only'];sc=[int(v) for st,_,_,v in matches if st=='S+B'];N=len(s['r']);assert len(bc)==len(sc)==N,(dd,len(bc),len(sc),N)
   track=uproot.open(next(dd.glob('higgsCombine*.root')))['limit'].arrays(library='np');last=np.flatnonzero(track['quantileExpected']==-1);assert len(last)==N,(dd,len(last),N);rows=[]
   for i in range(N):
    minima=[]
    for era in eras:
     for suffix in ['', '_CR']:
      coeff=[1]+[float(track[f'trackedParam_local_c{j}_{era}{suffix}'][last[i]]) for j in range(1,8)];roots=np.polynomial.chebyshev.chebroots(np.polynomial.chebyshev.chebder(coeff));pts=[-1,1]+[v.real for v in roots if abs(v.imag)<1e-8 and abs(v.real)<1];minima.append(float(min(np.polynomial.chebyshev.chebval(pts,coeff))))
    pos=min(minima);rh=float(s['r'][i]);err=float(s['rErr'][i]);dn=float(s['nll_nll0'][i]);valid=int(b['fit_status'][i])==0 and int(s['fit_status'][i])==0 and np.isfinite(rh+dn) and abs(rh)<99.9 and dn<=.01 and pos>=-1e-6;ve=valid and bc[i]==3 and sc[i]==3 and np.isfinite(err) and err>0
    row=dict(mass=mass,scenario=label,injection=truth,model='fixed7',toy=i+1,r=rh,error=err,b_status=int(b['fit_status'][i]),s_status=int(s['fit_status'][i]),b_cov=bc[i],s_cov=sc[i],delta_nll=dn,Z=float(np.sqrt(max(0,-2*dn))) if rh>0 else 0,continuum_min=pos,valid_point=bool(valid),valid_error=bool(ve),pull=(rh-truth)/err if ve else None);rows.append(row)
   if kind=='asimov':asimov.extend(rows);continue
   toyrows.extend(rows);v=[r for r in rows if r['valid_point']];ve=[r for r in rows if r['valid_error']];vals=np.array([r['r'] for r in v]);n=len(v);summary.append(dict(mass=mass,scenario=label,injection=truth,model='fixed7',n_total=N,n_valid=n,n_valid_errors=len(ve),mean_r=float(vals.mean()) if n else None,bias=float(vals.mean()-truth) if n else None,bias_mc_error=float(vals.std(ddof=1)/np.sqrt(n)) if n>1 else None,sd_r=float(vals.std(ddof=1)) if n>1 else None,median_error=float(np.median([r['error'] for r in ve])) if ve else None,mean_pull=float(np.mean([r['pull'] for r in ve])) if ve else None,sd_pull=float(np.std([r['pull'] for r in ve],ddof=1)) if len(ve)>1 else None,n_Zge3=sum(r['Z']>=3 for r in v),fraction_Zge3=sum(r['Z']>=3 for r in v)/n if n else None))
for name,rows in [('summary',summary),('toys',toyrows),('asimov_summary',asimov)]:
 (R/(name+'.json')).write_text(json.dumps(rows,indent=2))
 if rows:
  with (R/(name+'.csv')).open('w') as f:w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
status=dict(completed_batches=len(summary),expected_batches=42,toy_fits=len(toyrows),expected_toy_fits=1260,valid_points=sum(r['valid_point'] for r in toyrows),valid_errors=sum(r['valid_error'] for r in toyrows),asimov_fits=len(asimov),expected_asimov_fits=42,negative_continua=sum(r['continuum_min']< -1e-6 for r in toyrows));(R/'progress.json').write_text(json.dumps(status,indent=2))
combined=[dict(a,model='fixed3' if a['model']=='fixed' else 'sliding3') for a in oldsum]+summary;(R/'comparison.json').write_text(json.dumps(combined,indent=2))
styles=[('fixed3','tab:red','Grado3, range fisso'),('sliding3','tab:blue','Grado3, sliding'),('fixed7','tab:green','Grado7, range fisso')]
fig,axes=plt.subplots(2,2,figsize=(14,9))
for i,label in enumerate(['zero','injected']):
 for model,color,title in styles:
  a=sorted([s for s in combined if s['scenario']==label and s['model']==model and s['n_valid']>1],key=lambda s:s['mass']);x=[s['mass'] for s in a];axes[i,0].errorbar(x,[s['bias'] for s in a],yerr=[s['bias_mc_error'] for s in a],fmt='o-',ms=3,lw=1,capsize=2,color=color,label=title);a=[s for s in a if s['median_error']];axes[i,1].errorbar([s['mass'] for s in a],[s['bias']/s['median_error'] for s in a],yerr=[s['bias_mc_error']/s['median_error'] for s in a],fmt='o-',ms=3,lw=1,capsize=2,color=color,label=title)
 for j,ax in enumerate(axes[i]):ax.axhline(0,color='gray');ax.set_xlabel('Massa [GeV]');ax.set_ylabel('Bias in r' if j==0 else 'Bias / mediana errore Hesse');ax.set_title('Nessun segnale iniettato' if label=='zero' else 'Segnale iniettato al limite atteso di riferimento');ax.legend(fontsize=8);ax.grid(alpha=.2)
fig.suptitle('Stessi toy da Chebyshev7 + Y: confronto dei tre fit\nMedie solo sui fit validi; barre = errore MC sulla media; campioni selezionati diversi se falliscono i fit');fig.tight_layout();fig.savefig(R/'bias_overlay.png',dpi=150);fig.savefig(R/'bias_overlay.pdf');plt.close(fig)
fig,axes=plt.subplots(1,2,figsize=(13,4.5))
for model,color,title in styles:
 a=sorted([s for s in combined if s['model']==model and s['scenario']=='zero' and s['n_valid']],key=lambda s:s['mass']);axes[0].plot([s['mass'] for s in a],[s['fraction_Zge3'] for s in a],'o-',ms=3,color=color,label=title)
 for label,ls in [('zero','-'),('injected','--')]:
  a=sorted([s for s in combined if s['model']==model and s['scenario']==label],key=lambda s:s['mass']);axes[1].plot([s['mass'] for s in a],[s['n_valid']/s['n_total'] for s in a],ls,color=color,label=title+' '+label)
axes[0].set_ylabel('Frazione toy validi con Z locale >=3');axes[0].set_title('Generazione senza segnale');axes[1].set_ylabel('Frazione fit utilizzabili');axes[1].set_title('Qualita numerica (linea continua r=0, tratteggio iniezione)')
for ax in axes:ax.set_xlabel('Massa [GeV]');ax.set_ylim(-.03,1.03);ax.grid(alpha=.2);ax.legend(fontsize=7)
fig.tight_layout();fig.savefig(R/'quality_overlay.png',dpi=150);fig.savefig(R/'quality_overlay.pdf');plt.close(fig)
html='<meta charset="utf-8"><meta http-equiv="refresh" content="60"><style>body{font:17px system-ui;max-width:1400px;margin:25px auto}img{width:100%}</style><h1>Fit grado7 a range fisso: confronto appaiato</h1>'+f'<p>{status["completed_batches"]}/42 batch, {status["toy_fits"]}/1260 fit toy; {status["valid_points"]} utilizzabili. Stessi toy del test precedente, nessuna nuova generazione. Grado7 anche nella CR. Curve medie sui soli fit validi: controllare le esclusioni prima di interpretare differenze.</p><a href="REPORT.md">Report</a> · <a href="summary.csv">CSV</a><img src="bias_overlay.png"><img src="quality_overlay.png">';(R/'index.html').write_text(html);print(status)
