from pathlib import Path
import json,csv,hashlib
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
P=Path.cwd();R=P/'runs/20261030_sr6_cr4_scan';O=P/'tmp/pdfs/sr_cr_report';eras=['2022','2022EE','2023','2023BPix'];rows=[]
for k in range(60,91):
 m=k/10;a=json.loads((R/f'm{m:.1f}/signal/audit.json').read_text());b=json.loads((R/f'm{m:.1f}/control/audit.json').read_text())
 for e,x,y in zip(eras,a['eras'],b['eras']):
  fs=x['processes']['Zd']['fraction'];fc=y['processes']['Zd']['fraction'];rows.append(dict(mass=m,era=e,SR_fraction=fs,CR_fraction=fc,unused_fraction=1-fs-fc,CR_fraction_of_used=fc/(fs+fc),full_rate_r1=x['processes']['Zd']['full_rate_r1']))
(O/'signal_acceptance.json').write_text(json.dumps(rows,indent=2))
fig,axs=plt.subplots(1,2,figsize=(12,4.1));cols=['#0072B2','#D55E00','#009E73','#CC79A7']
for e,col in zip(eras,cols):
 a=[v for v in rows if v['era']==e];m=[v['mass'] for v in a];axs[0].plot(m,[100*v['CR_fraction'] for v in a],label=e,color=col);axs[1].plot(m,[100*v['SR_fraction'] for v in a],color=col)
axs[0].set_ylabel('Segnale nella CR / segnale in 6–12 [%]');axs[1].set_ylabel('Segnale nella SR / segnale in 6–12 [%]')
for ax in axs:ax.set_xlabel('Ipotesi di massa [GeV]');ax.grid(alpha=.2)
axs[0].legend();axs[0].set_title('Contaminazione della CR: modellata nella likelihood');axs[1].set_title('Accettanza nella SR: nessuna rinormalizzazione a 100%');fig.tight_layout();fig.savefig(O/'acceptance.png',dpi=190);plt.close(fig)
fig,ax=plt.subplots(figsize=(12,2.6));ax.broken_barh([(6,3.0933333333)],(1,0.55),facecolors='#0072B2');ax.broken_barh([(9.0933333333,1.4133333333)],(1,0.55),facecolors='#E69F00');ax.broken_barh([(10.5066666667,1.4933333333)],(1,0.55),facecolors='#dddddd');ax.text(7.55,1.27,'SR: 232 bin',ha='center',va='center',color='white',weight='bold');ax.text(9.8,1.27,'CR: 106 bin',ha='center',va='center',weight='bold');ax.text(11.25,1.27,'Non usata',ha='center',va='center');ax.plot([6,9],[.6,.6],color='#009E73',lw=6);ax.text(7.5,.3,'Ipotesi testate: 6,0–9,0 GeV, passo 0,1',ha='center');ax.set_xticks([6,9.0933333333,10.5066666667,12],['6','9,0933','10,5067','12']);ax.tick_params(axis='x',labelsize=10);ax.set_ylim(0,1.9);ax.set_xlim(5.9,12.1);ax.set_yticks([]);ax.set_xlabel('Massa ricostruita m(ee) [GeV]');ax.spines[['top','left','right']].set_visible(False);fig.tight_layout();fig.savefig(O/'regions.png',dpi=190);plt.close(fig)
# Audit of nuisance boundaries from the final parameter tables, without new fits.
flags=[]
for d in R.glob('m*'):
 p=d/'joint/parameters.csv'
 if not p.exists():continue
 for a in csv.DictReader(p.open()):
  if a['stage']!='s' or not ('Ymean_' in a['name'] or 'bkgScale_upsilon' in a['name']):continue
  v,lo,hi=map(float,[a['value'],a['minimum'],a['maximum']])
  if min(v-lo,hi-v)<1e-4*(hi-lo):flags.append(dict(mass=d.name,name=a['name'],value=v))
(O/'upsilon_boundaries.json').write_text(json.dumps(flags,indent=2));print('Signal fractions',len(rows),'Y boundary flags',len(flags))
f=json.loads((P/'runs/20261030_sr_bernstein7_bias/summary.json').read_text());b=json.loads((P/'runs/20261030_sr_cheb7_bern_bias/summary.json').read_text());print('reverse keys',b[0].keys())
fig,ax=plt.subplots(figsize=(10,4.3));ax.plot([x['degree'] for x in f],[x['max_asimov_abs_bias_sigma'] for x in f],'o-',label='Bernstein7 → Chebyshev',color='#0072B2')
key=next(k for k in b[0] if 'max' in k and 'bias' in k)
ax.plot([x['degree'] for x in b],[x[key] for x in b],'s-',label='Chebyshev7 → Bernstein',color='#D55E00');ax.axhline(.2,color='gray',ls='--',label='Criterio esplorativo 0,2σ');ax.set_yscale('symlog',linthresh=.02);ax.set_ylim(0,10);ax.set_xticks(range(3,8));ax.set_xlabel('Grado del polinomio nella SR');ax.set_ylabel('Massimo |bias| / σ(r)');ax.set_title('Selezione storica condizionale: 21 masse 6,5–8,5 GeV, due iniezioni\nNuisance del segnale nominali; per Bernstein3–5 alcuni fit non validi');ax.grid(alpha=.2);ax.legend();fig.tight_layout();fig.savefig(O/'selection_degrees.png',dpi=180);plt.close(fig)
