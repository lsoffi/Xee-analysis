from pathlib import Path
import json,csv
import numpy as np,uproot
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
P=Path.cwd();R=P/'runs/20261029_polynomial_comparison';eras=['2022','2022EE','2023','2023BPix'];models=[('Chebyshev 2','20261029_upsilon_control','#0072B2','-'),('Bernstein 2','20261029_upsilon_control_bern2','#009E73',':'),('Chebyshev 3','20261029_upsilon_control_cheb3','#D55E00','-'),('Bernstein 3','20261029_upsilon_control_bern3','#CC79A7','--')];rows=[];arrays={}
for label,folder,color,style in models:
 d=P/'runs'/folder;a=json.loads((d/'results.json').read_text());v=a['joint']
 stable=d/'joint/limits_stability'
 if stable.exists():
  import re
  with uproot.open(next(stable.glob('higgsCombine*.root'))) as f:
   vv=f['limit']['limit'].array(library='np');qq=f['limit']['quantileExpected'].array(library='np');v['limits']={k:float(vv[np.argmin(abs(qq-q))]) for k,q in [('observed',-1),('expected_2p5',.025),('expected_16',.16),('expected_50',.5),('expected_84',.84),('expected_97p5',.975)]}
  pairs=re.findall(r'At r =\s*([\d.eE+-]+):[^\n]*CLs\s*=\s*([\d.eE+-]+)',(stable/'combine.log').read_text());rr,cls=min(pairs,key=lambda t:abs(float(t[0])-v['limits']['observed']));assert abs(float(cls)-.05)<.005
 rows.append(dict(model=label,r=v['r'],Z=v['Z'],p0=v['p0'],**v['limits'],quality=a['quality']));arrays[label]={}
 for region in ['SR','CR']:
  for part in ['total','dy']:
   arrays[label][region+part]=[]
   with uproot.open(d/'joint/fitdiagnostics/fitDiagnostics_local.root') as f:
    for i in range(1,5):
     y,e=f[f'shapes_fit_s/{region}_ch{i}/{part}'].to_numpy();arrays[label][region+part].append((y*np.diff(e),e))
for region in ['SR','CR']:
 for part in ['total','dy']:
  fig,axs=plt.subplots(2,2,figsize=(12,8),sharex=True)
  for ax,i,era in zip(axs.flat,range(4),eras):
   for label,folder,color,style in models:
    y,e=arrays[label][region+part][i];ax.plot((e[:-1]+e[1:])/2,y,lw=2,ls=style,color=color,label=label)
   ax.set_title(era);ax.set_xlabel('Massa [GeV]');ax.set_ylabel('Eventi / 13,33 MeV');ax.legend(fontsize=9);ax.grid(alpha=.15)
  fig.suptitle(('Totale post-fit S+B' if part=='total' else 'Solo continuum post-fit S+B')+f' — {region} — segnale a 7,5 GeV');fig.tight_layout(rect=[0,0,1,.95]);fig.savefig(R/f'overlay_{region}_{part}.png',dpi=160);fig.savefig(R/f'overlay_{region}_{part}.pdf');plt.close(fig)
# Quantify agreement rather than claim identical models from rounded numbers.
agreement={}
for n in [2,3]:
 agreement[str(n)]={}
 for key in ['SRtotal','CRtotal','SRdy','CRdy']:
  agreement[str(n)][key]=[float(max(abs(a[0]-b[0])/np.maximum(abs(a[0]),1e-10))) for a,b in zip(arrays[f'Chebyshev {n}'][key],arrays[f'Bernstein {n}'][key])]
(R/'shape_agreement.json').write_text(json.dumps(agreement,indent=2));(R/'summary.json').write_text(json.dumps(rows,indent=2))
with (R/'summary.csv').open('w') as f:
 keys=[k for k in rows[0] if k!='quality'];w=csv.DictWriter(f,keys,extrasaction='ignore');w.writeheader();w.writerows(rows)
fig,ax=plt.subplots(figsize=(9,5))
for i,row in enumerate(rows):
 ax.plot([row['expected_2p5'],row['expected_97p5']],[i,i],color='#F6D743',lw=17);ax.plot([row['expected_16'],row['expected_84']],[i,i],color='#69BE46',lw=17);ax.plot(row['expected_50'],i,'k|',ms=23,mew=2);ax.plot(row['observed'],i,'o',color='#0072B2',label='Osservato' if i==0 else None)
ax.set_yticks(range(4),[x['model'] for x in rows]);ax.set_xlabel('Limite 95% CLs su r');ax.set_title('7,5 GeV — SR + controllo Υ\nMediana attesa: nero; bande attese: 68% e 95%');ax.grid(axis='x',alpha=.2);ax.legend();fig.tight_layout();fig.savefig(R/'limits.png',dpi=160);fig.savefig(R/'limits.pdf');plt.close(fig)
body='<meta charset="utf-8"><style>body{font:17px system-ui;max-width:1300px;margin:auto}img{width:100%}td,th{padding:10px}</style><h1>Chebyshev e Bernstein: gradi 2 e 3</h1><p>Massa7,5; SR6,373–7,920, CR8,8–10,8 GeV. Gaussiana sigma-core DCB, quattro ere, binning originale, Υ originali con normalizzazioni condivise floating. Bernstein positivo, primo coefficiente fissato a1, altri rapporti in [1e-6,20]. Stesso numero di parametri di forma del Chebyshev; famiglie non necessariamente identiche ai bordi dei domini. Limiti grado2 dal controllo numerico strategy2/tol1e-5/rMax1, primari conservati. Nessun nuovo bias test/GOF.</p><table><tr><th>Modello</th><th>Z locale</th><th>p0</th><th>Limite osservato</th><th>Atteso</th></tr>'
for row,model in zip(rows,models):
 body+=f'<tr><td>{row["model"]}</td><td>{row["Z"]:.6f}</td><td>{row["p0"]:.6g}</td><td>{row["observed"]:.6g}</td><td>{row["expected_50"]:.6g}</td></tr>'
body+='</table><h2>Fit e componenti</h2>'
for label,folder,_,_ in models:body+=f'<p><a href="../{folder}/index.html">{label}: fit SR/CR, componenti e parametri</a></p>'
for name in ['overlay_SR_total','overlay_SR_dy','overlay_CR_total','overlay_CR_dy','limits']:body+=f'<h2>{name}</h2><img src="{name}.png">'
(R/'index.html').write_text(body);print(json.dumps(rows,indent=2));print('shape agreement',json.dumps(agreement))
