from pathlib import Path
import json,csv,math,datetime
from zoneinfo import ZoneInfo
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
base=Path(__file__).resolve().parent.parent
out=Path('/eos/user/s/soffi/www/Xee/261106/scan_SR_CR/partial_reports');out.mkdir(exist_ok=True)
rows=[]
for folder in ['limits_recommended_envelope/results','high_c4l4_limits_bias/limits']:
 for p in (base/folder).glob('m*.json'):
  d=json.load(open(p));r={'mass':d['mass'],'source':str(p)}
  for key,q in [('observed',-1),('minus2',.025),('minus1',.16),('median',.5),('plus1',.84),('plus2',.975)]:
   vv=[v['limit'] for v in d['quantiles'] if abs(v['quantile']-q)<1e-5];assert len(vv)==1 and math.isfinite(vv[0]) and vv[0]>0;r[key]=vv[0]
  assert r['minus2']<=r['minus1']<=r['median']<=r['plus1']<=r['plus2'];rows.append(r)
rows.sort(key=lambda r:r['mass']);lookup={round(r['mass']*10):r for r in rows};assert len(lookup)==len(rows)
stamp=datetime.datetime.now(ZoneInfo('Europe/Madrid')).strftime('%d %b %Y, %H:%M CET')
fig,ax=plt.subplots(figsize=(12,6.4));first=True
for ids,label,x in [(range(4,7),'Low SR: C3+L3',.5),(range(12,28),'Middle SR: C4+L4',1.95),(range(41,91),'High SR: C4+L4',6.5)]:
 ids=list(ids);xx=np.array(ids)/10
 vals={key:np.array([lookup.get(i,{}).get(key,np.nan) for i in ids]) for key in ['observed','minus2','minus1','median','plus1','plus2']}
 ax.fill_between(xx,vals['minus2'],vals['plus2'],color='#f5cc42',label='Expected ±2σ' if first else None)
 ax.fill_between(xx,vals['minus1'],vals['plus1'],color='#75bd63',label='Expected ±1σ' if first else None)
 ax.plot(xx,vals['median'],'k--',lw=1.5,label='Median expected' if first else None)
 ax.plot(xx,vals['observed'],'ko-',ms=3,lw=1.2,label='Observed' if first else None)
 ax.text(x,.97,label,transform=ax.get_xaxis_transform(),ha='center',va='top',fontsize=10)
 first=False
for lo,hi,label in [(.65,1.155,'ω/φ CR'),(2.8,4,'Charmonium CR')]:
 ax.axvspan(lo,hi,color='.92',zorder=-2);ax.text((lo+hi)/2,.5,label,rotation=90,ha='center',va='center',transform=ax.get_xaxis_transform(),color='.45',fontsize=10)
missing=[i/10 for i in list(range(4,7))+list(range(12,28))+list(range(41,91)) if i not in lookup]
for m in missing:ax.plot(m,.015,marker='x',color='#a50f15',transform=ax.get_xaxis_transform(),clip_on=False)
ax.set(yscale='log',xlim=(.28,9.12),xlabel='Signal mass hypothesis [GeV]',ylabel='95% CL upper limit on signal strength r',title=f'Sliding-window limits: low + middle + high SR ({len(rows)}/69 points)')
ax.grid(which='both',axis='y',alpha=.18);ax.legend(loc='lower left',fontsize=9)
fig.text(.08,.047,stamp+' | Missing masses: '+(', '.join(map(str,missing)) if missing else 'none')+' GeV. No interpolation across missing points or CRs.',fontsize=9)
fig.text(.08,.018,'High SR provisional; observed limits near Υ under numerical investigation. Low-SR GoF tension at 0.6 GeV.',fontsize=9,color='#8b2525')
fig.tight_layout(rect=(0,.08,1,1))
for name in ['limits','limits_all_sr']:
 for ext in ['png','pdf']:fig.savefig(out/f'{name}.{ext}',dpi=170)
with (out/'limits_all_sr.csv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
(out/'limits_all_sr.json').write_text(json.dumps({'updated':stamp,'points':len(rows),'missing':missing,'rows':rows},indent=2))
print(stamp,len(rows),missing)
