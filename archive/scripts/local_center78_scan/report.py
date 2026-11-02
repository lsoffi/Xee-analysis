import sys
MASS=sys.argv[1]
MODEL=sys.argv[2]
import pathlib,json,math,re,csv
import numpy as np,uproot
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import chi2,norm
P=pathlib.Path('/Users/imb/Documents/Codex/Analisi-statistica-Xee');R=P/'runs/20261029_local_center78_scan'/MODEL/('m'+MASS);R.mkdir(parents=True,exist_ok=True);audit=json.loads((R/'audit.json').read_text());quality=json.loads((R/'quality.json').read_text());result={};fig=plt.figure(figsize=(14,12));outer=fig.add_gridspec(2,2,hspace=.3,wspace=.18);totals={};stats=[]
def dev(d,m):return float(2*np.sum(m-d+np.where(d>0,d*np.log(np.maximum(d,1)/m),0)))
with uproot.open(R/'fitdiagnostics/fitDiagnostics_local.root') as f:
 t=f['tree_fit_sb'];result['r']={k:float(t[k].array(library='np')[0]) for k in ['r','rLoErr','rHiErr','rErr']}
 for i,(cell,a) in enumerate(zip(outer,audit['eras']),1):
  era=a['era'];gs=cell.subgridspec(3,1,height_ratios=[3,1,1],hspace=.08);ax=fig.add_subplot(gs[0]);res=fig.add_subplot(gs[1],sharex=ax);arr={}
  for st in ['b','s']:
   dr=f[f'shapes_fit_{st}/ch{i}'];v,edges=dr['total'].to_numpy();bw=np.diff(edges);x=(edges[:-1]+edges[1:])/2;d=dr['data'].values('y')*bw
   assert len(d)==76 and abs(d.sum()-a['events'])<.01
   arr[st]={k:dr[k].values()*bw for k in ['total','dy','Zd']}
  b=arr['b']['total'];s=arr['s']['total'];stats.append(dict(era=era,B_deviance=dev(d,b),SB_deviance=dev(d,s),B_df=(74 if MODEL=='exp' else 73),B_approx_p=float(chi2.sf(dev(d,b),74 if MODEL=='exp' else 73)),signal_events=float(arr['s']['Zd'].sum())))
  row=dict(mass=x,data=d,B=b,SB=s,continuum_SB=arr['s']['dy'],signal_SB=arr['s']['Zd']);
  for k,v in row.items():totals[k]=totals.get(k,0)+v if k!='mass' else v
  with (R/f'shapes_{era}.csv').open('w') as h:
   w=csv.writer(h);w.writerow(row);w.writerows(zip(*row.values()))
  ax.errorbar(x,d,yerr=np.sqrt(d),fmt='.',color='black',ms=3,lw=.7,label='Dati');ax.stairs(b,edges,color='#296ca1',ls='--',label='Solo fondo');ax.stairs(s,edges,color='#c53030',label='S+B');ax.stairs(arr['s']['dy'],edges,color='#dd9e22',ls=':',label='Fondo nel fit S+B');ax.set_title(era+f" — devianza B = {dev(d,b):.1f} / {74 if MODEL=='exp' else 73}");ax.set_ylabel('Eventi / 13,33 MeV');ax.legend(fontsize=8);ax.tick_params(labelbottom=False);ax.grid(alpha=.15)
  res.axhline(0,color='black',lw=.7);res.axhspan(-2,2,color='gray',alpha=.1);res.plot(x,(d-b)/np.sqrt(b),'.',color='#296ca1',label='B');res.plot(x,(d-s)/np.sqrt(s),'.',color='#c53030',label='S+B');res.set_ylim(-4.5,4.5);res.set_xlim(7.293333333333333,8.306666666666667);res.tick_params(labelbottom=False);res.set_ylabel('(D−M)/√M');res.grid(alpha=.15)
  sax=fig.add_subplot(gs[2],sharex=ax);sax.stairs(arr['s']['Zd'],edges,color='#188548',lw=1.5);sax.axhline(0,color='gray',lw=.7);sax.set_ylabel('Segnale / bin');sax.set_xlabel('Massa [GeV]');sax.grid(alpha=.15);sax.text(.02,.88,'DCB post-fit (con segno)',transform=sax.transAxes,fontsize=8)

fig.suptitle(f'Fit locale 7,2933–8,3067 GeV: {MODEL} + DCB workspace a {MASS} GeV\nQuattro ere combinate; binning originale; risoluzione MC',fontsize=14);fig.subplots_adjust(top=.90,bottom=.07)
for ext in ['png','pdf']:fig.savefig(R/f'fits.{ext}',dpi=145)
plt.close(fig)
fig,(ax,rx,sx)=plt.subplots(3,1,figsize=(10,9),sharex=True,gridspec_kw={'height_ratios':[3,1,1]});x=totals['mass'];d=totals['data'];ax.errorbar(x,d,yerr=np.sqrt(d),fmt='.',color='black',label='Dati, somma ere');ax.stairs(totals['B'],edges,color='#296ca1',ls='--',label='B-only');ax.stairs(totals['SB'],edges,color='#c53030',label='S+B');ax.legend();ax.set_ylabel('Eventi / 13,33 MeV');ax.set_title('Fit simultaneo nelle quattro ere — somma per visualizzazione')
for st,c in [('B','#296ca1'),('SB','#c53030')]:rx.plot(x,(d-totals[st])/np.sqrt(totals[st]),'.',color=c,label=st)
rx.axhline(0,color='gray');rx.set_ylabel('(D−M)/√M');sx.stairs(totals['signal_SB'],edges,color='#23954a',label='Componente segnale: negativa');sx.axhline(0,color='gray');sx.set_ylabel('Eventi segnale / bin');sx.set_xlabel('Massa [GeV]');sx.legend();sx.set_xlim(7.293333333333333,8.306666666666667)
for a in [ax,rx,sx]:a.grid(alpha=.15)
fig.tight_layout()
for ext in ['png','pdf']:fig.savefig(R/f'fit_sum.{ext}',dpi=145)
plt.close(fig)
for role in ['significance','limits']:
 with uproot.open(next((R/role).glob('higgsCombine*.root'))) as f:
  v=f['limit']['limit'].array(library='np')
  if role=='significance':result.update(Z=float(v[0]),p0=float(norm.sf(v[0])))
  else:
   q=f['limit']['quantileExpected'].array(library='np');result['limits']={k:float(v[np.argmin(abs(q-t))]) for k,t in [('observed',-1),('expected_2p5',.025),('expected_16',.16),('expected_50',.5),('expected_84',.84),('expected_97p5',.975)]};assert len(q)==6
log=(R/'limits/combine.log').read_text();matches=re.findall(r'At r =\s*([\d.eE+-]+):[^\n]*CLs\s*=\s*([\d.eE+-]+)',log);obs=result['limits']['observed'];r,cls=min(matches,key=lambda a:abs(float(a[0])-obs));result['limit_crossing']=dict(r=float(r),CLs=float(cls),verified=abs(float(cls)-.05)<.01 and abs(float(r)/obs-1)<.02);assert result['limit_crossing']['verified']
result['goodness']=stats;result['B_deviance']=sum(a['B_deviance'] for a in stats);result['B_df']=296 if MODEL=='exp' else 292;result['B_approx_p']=float(chi2.sf(result['B_deviance'],result['B_df']));result['SB_deviance']=sum(a['SB_deviance'] for a in stats);result['signal_events']=sum(a['signal_events'] for a in stats);result['quality']={k:{a:v for a,v in val.items() if a!='parameters'} for k,val in quality.items()};(R/'results.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))
