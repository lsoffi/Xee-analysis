from pathlib import Path
import json,csv,numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
R=Path('runs/20261030_cheb7_fit7_comparison');OLD=Path('runs/20261030_cheb7_bias_scan');new=json.loads((R/'toys.json').read_text());old={}
with (OLD/'toys.csv').open() as f:
 for r in csv.DictReader(f):old[(float(r['mass']),r['scenario'],int(r['toy']),r['model'])]=r
rows=[];styles=[('fixed3','tab:red','Grado3 fisso'),('sliding3','tab:blue','Grado3 sliding'),('fixed7','tab:green','Grado7 fisso')]
for mass,label in sorted(set((r['mass'],r['scenario']) for r in new)):
 selected=[]
 for r in new:
  if (r['mass'],r['scenario'])!=(mass,label):continue
  f=old[(mass,label,r['toy'],'fixed')];s=old[(mass,label,r['toy'],'sliding')]
  if r['valid_point'] and f['valid_point']=='True' and s['valid_point']=='True':selected.append((r,f,s))
 if not selected:continue
 truth=selected[0][0]['injection'];n=len(selected)
 for model,index in [('fixed7',0),('fixed3',1),('sliding3',2)]:
  vals=np.array([float(t[index]['r']) for t in selected]);rows.append(dict(mass=mass,scenario=label,model=model,n_common=n,injection=truth,bias=float(vals.mean()-truth),bias_mc_error=float(vals.std(ddof=1)/np.sqrt(n)) if n>1 else None))
(R/'common_sample.json').write_text(json.dumps(rows,indent=2))
if rows:
 with (R/'common_sample.csv').open('w') as f:w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
fig,axes=plt.subplots(1,2,figsize=(13,4.5))
for ax,label in zip(axes,['zero','injected']):
 for model,color,title in styles:
  a=[r for r in rows if r['scenario']==label and r['model']==model and r['n_common']>1];ax.errorbar([r['mass'] for r in a],[r['bias'] for r in a],yerr=[r['bias_mc_error'] for r in a],fmt='o-',ms=3,capsize=2,color=color,label=title)
 ax.axhline(0,color='gray');ax.set_title('Nessun segnale' if label=='zero' else 'Segnale iniettato');ax.set_xlabel('Massa [GeV]');ax.set_ylabel('Bias in r');ax.grid(alpha=.2);ax.legend(fontsize=8)
fig.suptitle('Confronto sullo stesso sottoinsieme di toy validi in tutti e tre i fit\nSelezione condizionata alla convergenza: non misura il comportamento dei toy falliti');fig.tight_layout();fig.savefig(R/'paired_bias_overlay.png',dpi=150);fig.savefig(R/'paired_bias_overlay.pdf');plt.close(fig)
