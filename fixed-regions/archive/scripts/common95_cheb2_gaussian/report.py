import sys,pathlib,json,re,csv
import numpy as np,uproot
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import norm
P=pathlib.Path('/Users/imb/Documents/Codex/Analisi-statistica-Xee');M=sys.argv[1];R=P/'runs/20261029_common95_cheb2_gaussian'/('m'+M);a=json.loads((R/'audit.json').read_text());q=json.loads((R/'quality.json').read_text());result={'mass':float(M),'range':a['range'],'bins':a['bins'],'quality':{k:{j:v for j,v in b.items() if j!='parameters'} for k,b in q.items()}}
for stage in ['b','s']:
 vals={v['name']:v['value'] for v in q[stage]['parameters']};mins=[]
 for era in ['2022','2022EE','2023','2023BPix']:
  c1=vals['local_c1_'+era];c2=vals['local_c2_'+era];ys=[1-c1+c2,1+c1+c2]
  if c2>0 and abs(c1)<4*c2:ys.append(1-c2-c1*c1/(8*c2))
  mins.append(min(ys))
 result['quality'][stage]['continuum_min']=min(mins)
result['fit_ok']=all(v['status']==0 and v['covQual']==3 and v['edm']<.01 and v['continuum_min']>0 for v in result['quality'].values())
fig=plt.figure(figsize=(14,11));outer=fig.add_gridspec(2,2,hspace=.30,wspace=.20)
with uproot.open(R/'fitdiagnostics/fitDiagnostics_local.root') as f:
 t=f['tree_fit_sb'];result['r']={k:float(t[k].array(library='np')[0]) for k in ['r','rErr']}
 cmd=json.loads((R/'fitdiagnostics/command.json').read_text());rlo=float(cmd[cmd.index('--rMin')+1]);rhi=float(cmd[cmd.index('--rMax')+1]);result['r_at_boundary']=min(result['r']['r']-rlo,rhi-result['r']['r'])<1e-4*(rhi-rlo);result['fit_ok']=result['fit_ok'] and not result['r_at_boundary']
 for i,(cell,ea) in enumerate(zip(outer,a['eras']),1):
  gs=cell.subgridspec(2,1,height_ratios=[3,1],hspace=.07);ax=fig.add_subplot(gs[0]);rx=fig.add_subplot(gs[1],sharex=ax)
  dr=f[f'shapes_fit_s/ch{i}'];v,edges=dr['total'].to_numpy();bw=np.diff(edges);x=(edges[:-1]+edges[1:])/2;data=dr['data'].values('y')*bw
  assert len(data)==a['bins'] and abs(sum(data)-ea['events'])<.05,(M,i,len(data),sum(data),ea['events'])
  shapes={k:dr[k].values()*bw for k in ['total','dy','Zd','upsilon1s','upsilon2s']};b=f[f'shapes_fit_b/ch{i}/total'].values()*bw
  ax.errorbar(x,data,yerr=np.sqrt(data),fmt='.',ms=2.5,color='black',lw=.6,label='Dati');ax.stairs(b,edges,color='gray',ls='--',label='B-only');ax.stairs(shapes['total'],edges,color='#c53030',label='S+B')
  for key,col,label in [('dy','#0072B2','Chebyshev 2'),('upsilon1s','#E69F00','Υ(1S)'),('upsilon2s','#CC79A7','Υ(2S)'),('Zd','#009E73','Gaussiana segnale')]:ax.stairs(shapes[key],edges,color=col,label=label,lw=1)
  ax.set_title(ea['era']);ax.set_ylabel('Eventi / 13,33 MeV');ax.legend(fontsize=7,ncol=2);ax.tick_params(labelbottom=False);ax.grid(alpha=.15)
  rx.axhline(0,color='gray',lw=.7);rx.plot(x,(data-shapes['total'])/np.sqrt(np.maximum(shapes['total'],1e-9)),'.',color='#c53030',ms=3);rx.set_ylabel('(D−M)/√M',fontsize=9);rx.set_xlabel('Massa [GeV]');rx.set_xlim(*a['range']);rx.grid(alpha=.15)
  with (R/f'shapes_{ea["era"]}.csv').open('w') as h:
   w=csv.writer(h);w.writerow(['mass','data','B','SB','cheb2','upsilon1s','upsilon2s','signal']);w.writerows(zip(x,data,b,shapes['total'],shapes['dy'],shapes['upsilon1s'],shapes['upsilon2s'],shapes['Zd']))
qtext='; '.join(f"{st}: status {v['status']}, cov {v['covQual']}" for st,v in result['quality'].items())
fig.suptitle(f"Chebyshev 2 + Gaussiana + Υ(1S) + Υ(2S), m = {M} GeV\nFinestra {a['range'][0]:.3f}–{a['range'][1]:.3f} GeV — {qtext}",fontsize=13);fig.subplots_adjust(top=.90,bottom=.07)
for ext in ['png','pdf']:fig.savefig(R/f'fit.{ext}',dpi=145)
plt.close(fig)
for role in ['limits','significance']:
 receipt=R/role/'receipt.json'
 if not receipt.exists() or json.loads(receipt.read_text())['exit_code']!=0:result[role+'_error']='execution failed or timeout';continue
 try:
  with uproot.open(next((R/role).glob('higgsCombine*.root'))) as f:
   v=f['limit']['limit'].array(library='np')
   if role=='significance':result.update(Z=float(v[0]),p0=float(norm.sf(v[0])))
   else:
    quant=f['limit']['quantileExpected'].array(library='np');assert len(v)==6
    result['limits']={k:float(v[np.argmin(abs(quant-t))]) for k,t in [('observed',-1),('expected_2p5',.025),('expected_16',.16),('expected_50',.5),('expected_84',.84),('expected_97p5',.975)]}
  if role=='limits':
   log=(R/role/'combine.log').read_text();matches=re.findall(r'At r =\s*([\d.eE+-]+):[^\n]*CLs\s*=\s*([\d.eE+-]+)',log);obs=result['limits']['observed'];r,cls=min(matches,key=lambda z:abs(float(z[0])-obs));result['limit_crossing']=dict(r=float(r),CLs=float(cls),verified=abs(float(cls)-.05)<.005 and abs(float(r)/obs-1)<.02)
 except Exception as e:result[role+'_error']=str(e)
(R/'results.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))
