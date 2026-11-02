from pathlib import Path
import sys,json,re
import numpy as np,uproot
from scipy.stats import norm
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
P=Path.cwd();M=sys.argv[1];R=P/'runs/20261030_sr7_cr3_scan'/('m'+M);J=R/'joint';q=json.loads((J/'quality.json').read_text());audit=json.loads((J/'audit.json').read_text());sa=json.loads((R/'signal/audit.json').read_text());ca=json.loads((R/'control/audit.json').read_text());eras=['2022','2022EE','2023','2023BPix'];result=dict(mass=float(M),range=sa['range'],CR_range=ca['range'],bins=sa['bins'],quality={k:{j:v for j,v in vv.items() if j!='parameters'} for k,vv in q.items()},CR_signal_fraction={e:a['processes']['Zd']['fraction'] for e,a in zip(eras,ca['eras'])})
with uproot.open(J/'fitdiagnostics/fitDiagnostics_local.root') as f:
 t=f['tree_fit_sb'];r=float(t['r'].array(library='np')[0]);dn=float(t['nll_nll0'].array(library='np')[0]);result.update(r=r,Z_from_FD=float(np.sqrt(max(0,-2*dn))) if r>0 else 0)
 for region in ([] if len(sys.argv)>2 and sys.argv[2]=='numeric' else ['SR','CR']):
  fig=plt.figure(figsize=(13,10));outer=fig.add_gridspec(2,2,hspace=.32,wspace=.2)
  for i,(era,cell) in enumerate(zip(eras,outer),1):
   dr=f[f'shapes_fit_s/{region}_ch{i}'];v,e=dr['total'].to_numpy();bw=np.diff(e);x=(e[:-1]+e[1:])/2;data=dr['data'].values('y')*bw;pred=v*bw;assert abs(data.sum()-audit['counts'][region+f'_ch{i}'])<.1
   parts={p:dr[p].values()*bw for p in ['dy','upsilon1s','upsilon2s','Zd']};assert np.max(abs(sum(parts.values())-pred)/np.maximum(abs(pred),1e-20))<5e-7
   gs=cell.subgridspec(2,1,height_ratios=[3,1],hspace=.06);ax=fig.add_subplot(gs[0]);rx=fig.add_subplot(gs[1],sharex=ax);ax.errorbar(x,data,yerr=np.sqrt(data),fmt='.',ms=2,color='black',lw=.5,label='Dati');ax.stairs(pred,e,color='#D55E00',label='Totale S+B');ax.stairs(f[f'shapes_fit_b/{region}_ch{i}/total'].values()*bw,e,ls='--',color='gray',label='B-only')
   for proc,col,lab in [('dy','#0072B2','Chebyshev '+('7' if region=='SR' else '3')),('upsilon1s','#E69F00','Υ(1S)'),('upsilon2s','#CC79A7','Υ(2S)'),('Zd','#009E73','Segnale DCB')]:ax.stairs(parts[proc],e,color=col,label=lab)
   ax.set_title(era);ax.set_ylabel('Eventi / 13,33 MeV');ax.tick_params(labelbottom=False);ax.legend(fontsize=7,ncol=2);ax.grid(alpha=.15);rx.axhline(0,color='gray',lw=.7);rx.plot(x,(data-pred)/np.sqrt(np.maximum(pred,1e-10)),'.',color='#D55E00',ms=2.5);rx.set_ylabel('(D−M)/√M');rx.set_xlabel('Massa [GeV]');rx.set_xlim(e[0],e[-1]);rx.grid(alpha=.15)
  fig.suptitle(f'Massa {M} GeV — {region} — fit simultaneo SR Chebyshev7 / CR Chebyshev3 + DCB + Υ\nDCB presente anche nella CR; binning originale');fig.subplots_adjust(top=.90,bottom=.07)
  for ext in ['png','pdf']:fig.savefig(R/f'fit_{region}.{ext}',dpi=140)
  plt.close(fig)
with uproot.open(next((J/'significance').glob('higgsCombine*.root'))) as f:result['Z']=float(f['limit']['limit'].array(library='np')[0]);result['p0']=float(norm.sf(result['Z']))
with uproot.open(next((J/'limits').glob('higgsCombine*.root'))) as f:
 t=f['limit'];vv=t['limit'].array(library='np');qq=t['quantileExpected'].array(library='np');assert len(vv)==6;result['limits']={k:float(vv[np.argmin(abs(qq-v))]) for k,v in [('observed',-1),('expected_2p5',.025),('expected_16',.16),('expected_50',.5),('expected_84',.84),('expected_97p5',.975)]}
obs=result['limits']['observed'];pairs=re.findall(r'At r =\s*([\d.eE+-]+):[^\n]*CLs\s*=\s*([\d.eE+-]+)',(J/'limits/combine.log').read_text());rr,cls=min(pairs,key=lambda x:abs(float(x[0])-obs));result['limit_crossing']=dict(r=float(rr),CLs=float(cls),verified=abs(float(cls)-.05)<.005 and abs(float(rr)/obs-1)<.02)
result['fit_ok']=all(v['status']==0 and v['covQual']==3 and v['continuum_min']>0 for v in q.values());result['significance_consistent']=abs(result['Z']-result['Z_from_FD'])<.1;result['numerical_ok']=result['fit_ok'] and result['significance_consistent'] and abs(r)<99.9 and result['limit_crossing']['verified'];(R/'results.json').write_text(json.dumps(result,indent=2));print(M,result['Z'],result['limits']['observed'],result['numerical_ok'],flush=True)
