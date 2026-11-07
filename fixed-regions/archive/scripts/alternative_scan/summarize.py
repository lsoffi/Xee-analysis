import pathlib,json,csv,re,time,argparse,math,traceback
import numpy as np,uproot
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import norm
P=pathlib.Path('/Users/imb/Documents/Codex/Analisi-statistica-Xee');R=P/'runs/20261029_alternative_background_scan';ap=argparse.ArgumentParser();ap.add_argument('--watch',action='store_true');args=ap.parse_args();(R/'summary').mkdir(exist_ok=True);(R/'plots').mkdir(exist_ok=True)
LABEL={'logquad':'Potenza con curvatura','expquad':'Esponenziale quadratico','doubleexp':'Somma di due esponenziali','bernstein7':'Bernstein 7 identificabile'};COLORS={'logquad':'#176caa','expquad':'#df7319','doubleexp':'#9d409c','bernstein7':'#278044'};ERAS=['2022','2022EE','2023','2023BPix'];cache={};cp=R/'summary/processed.json'
if cp.exists():cache=json.loads(cp.read_text())
def csvwrite(path,rows):
 if not rows:return
 keys=list(dict.fromkeys(k for r in rows for k in r));tmp=path.with_suffix(path.suffix+'.tmp')
 with tmp.open('w') as f:w=csv.DictWriter(f,fieldnames=keys);w.writeheader();w.writerows(rows)
 tmp.replace(path)
def process(rp):
 t=json.loads(rp.read_text());wd=rp.parent;row={k:t[k] for k in ['family','mass','mass_label','role','exit_code','elapsed_seconds']};row['directory']=str(wd.relative_to(R));row['receipt_mtime']=rp.stat().st_mtime;row['valid']=False
 if t['exit_code']!=0:row['error']=t.get('error');return row
 log=(wd/'combine.log').read_text(errors='replace');row['warning_lines']=len(re.findall(r'\[WARNING\]|WARNING:',log))
 if t['role'] in ['limits','significance']:
  roots=list(wd.glob('higgsCombine*.root'));assert len(roots)==1
  with uproot.open(roots[0]) as f:
   values=f['limit']['limit'].array(library='np')
   if t['role']=='significance':
    assert len(values)==1;z=float(values[0]);row.update(Z=z,p0=float(norm.sf(z)),log10_p0=float(norm.logsf(z)/np.log(10)),valid=bool(np.isfinite(z)))
   else:
    qs=f['limit']['quantileExpected'].array(library='np');targets={'observed':-1,'expected_2p5':.025,'expected_16':.16,'expected_50':.5,'expected_84':.84,'expected_97p5':.975}
    for k,q in targets.items():
     v=values[np.isclose(qs,q,atol=.004)];row[k]=float(v[0]) if len(v)==1 else None
    row['valid']=all(row[k] is not None and np.isfinite(row[k]) and row[k]>0 for k in targets)
    ex=[row[k] for k in list(targets)[1:]];row['expected_ordered']=all(v is not None for v in ex) and all(a<=b for a,b in zip(ex,ex[1:]));row['fit_failure_messages']=len(re.findall(r'failed to converge|Minimization failed|fit failed|Failed to find',log,re.I))
  return row
 matches=re.findall(r'Fit (B-only|S\+B), status =\s*(-?\d+), numBadNLL =\s*(\d+), covariance quality =\s*(-?\d+)',log)
 for stage,status,bad,cov in matches:
  st='b' if stage=='B-only' else 's';row.update({st+'_status':int(status),st+'_covQual':int(cov),st+'_bad_nll':int(bad)})
 with uproot.open(wd/'fitDiagnostics_alternative.root') as f:
  tr=f['tree_fit_sb'];pars=[];bounds=[]
  for stage in ['tree_fit_b','tree_fit_sb']:
   tree=f[stage]
   for key in tree.keys():
    if key.startswith(('ALT_','CMS_','lumi')) or key=='r':
     v=float(tree[key].array(library='np')[0]);pars.append(dict(stage=stage,parameter=key,value=v))
     if key.startswith('ALT_') and stage=='tree_fit_sb':
      tail=key.split('_')[-1];fam=t['family'];lo,hi=(-100,100)
      if fam=='logquad':lo,hi=(-100,50) if tail=='c1' else (-200,200)
      elif fam=='expquad':lo,hi=(-30,30)
      elif fam=='doubleexp':lo,hi={'kslow':(.001,30),'gap':(.001,60),'fraction':(0,1)}[tail]
      if min(v-lo,hi-v)<=1e-4*(hi-lo):bounds.append(dict(parameter=key,value=v,low=lo,high=hi))
  for k in ['r','rLoErr','rHiErr','rErr','nll_nll0']:row[k]=float(tr[k].array(library='np')[0])
  row['Z_from_FD']=math.sqrt(max(0,-2*row['nll_nll0'])) if row['r']>0 else 0;row['r_interval_at_bound']=row['r']-abs(row['rLoErr'])<=-99.999 or row['r']+abs(row['rHiErr'])>=99.999;row['r_bestfit_at_bound']=abs(row['r'])>=99.999;row['shape_bounds']=bounds;row['valid']=row.get('b_status')==0 and row.get('s_status')==0 and row.get('b_covQual')==3 and row.get('s_covQual')==3
  csvwrite(wd/'parameters_from_trees.csv',pars);fig=plt.figure(figsize=(14,10));outer=fig.add_gridspec(2,2,hspace=.28,wspace=.2);tables=[];normalizations=[];deviance=0;ndata=0
  for i,(cell,era) in enumerate(zip(outer,ERAS),1):
   gs=cell.subgridspec(2,1,height_ratios=[3,1],hspace=.05);ax=fig.add_subplot(gs[0]);res=fig.add_subplot(gs[1],sharex=ax);arrays={}
   for st in ['fit_b','fit_s']:
    dr=f[f'shapes_{st}/ch{i}'];v,edges=dr['total'].to_numpy();w=np.diff(edges);x=(edges[:-1]+edges[1:])/2;d=dr['data'].values('y')*w;table=dict(mass_GeV=x,low=edges[:-1],high=edges[1:],data=d)
    assert len(d)==450 and abs(edges[0]-6)<1e-9 and abs(edges[-1]-12)<1e-9
    for key in ['total','dy','upsilon1s','upsilon2s','Zd']:
     table[key]=dr[key].values()*w;normalizations.append(dict(era=era,stage=st,component=key,events=float(table[key].sum())))
    arrays[st]=table;csvwrite(wd/f'{era}_{st}_shapes.csv',[{k:float(v[n]) for k,v in table.items()} for n in range(450)])
   b=arrays['fit_b'];s=arrays['fit_s'];D=s['data'];B=b['total'];S=s['total'];ndata+=float(D.sum());mask=(B>0)&(D>=0);term=B[mask]-D[mask];pos=D[mask]>0;term[pos]+=D[mask][pos]*np.log(D[mask][pos]/B[mask][pos]);deviance+=2*float(term.sum())
   ax.errorbar(x,D,yerr=np.sqrt(D),fmt='.',ms=1.9,lw=.4,color='black',label='Dati',zorder=6);ax.stairs(B,edges,baseline=None,color='gray',ls='--',lw=1.3,label='Totale B-only')
   for key,color,label,ls in [('total','#c62828','Totale S+B','-'),('dy','#246eae','Continuo S+B','--'),('upsilon1s','#e88a11','Υ(1S)','-.'),('upsilon2s','#9442aa','Υ(2S)',':'),('Zd','#159447','Segnale','-')]:
    y=s[key]
    if np.any(y<0):y=np.abs(y);label='|'+label+'|, negativo';ls=':'
    if np.any(y>0):ax.stairs(np.where(y>0,y,np.nan),edges,baseline=None,color=color,label=label,ls=ls,lw=1.1)
   ax.set_yscale('log');ax.set_ylim(.01,max(D)*3);ax.set_xlim(6,12);ax.set_title(era);ax.set_ylabel('Eventi / 13,33 MeV');ax.tick_params(labelbottom=False);ax.grid(alpha=.12,which='both');ax.legend(fontsize=6.5,ncol=2,loc='lower left')
   res.axhline(0,color='gray',lw=.7);res.axhspan(-2,2,color='gray',alpha=.1);res.plot(x,(D-B)/np.sqrt(np.maximum(B,1e-12)),'.',ms=1.7,color='gray',label='B-only');res.plot(x,(D-S)/np.sqrt(np.maximum(S,1e-12)),'.',ms=1.7,color='#c62828',label='S+B');res.set_ylim(-6,6);res.set_ylabel('(D−M)/√M',fontsize=9);res.set_xlabel('Massa [GeV]');res.grid(alpha=.12)
  row['data_events']=ndata;row['B_poisson_deviance']=deviance;row['diagnostic_warning']=not row['valid'] or bool(bounds) or row['r_interval_at_bound'];csvwrite(wd/'normalizations.csv',normalizations)
  fig.suptitle(f"{LABEL[t['family']]} — massa ipotizzata {t['mass_label']} GeV\nFit 6–12 GeV, quattro ere, Υ libere; B: {row.get('b_status')}/{row.get('b_covQual')}, S+B: {row.get('s_status')}/{row.get('s_covQual')} (status/covQual)",fontsize=12,y=.99);fig.subplots_adjust(top=.91,bottom=.065,left=.065,right=.98)
  for ext in ['png','pdf']:fig.savefig(wd/('fits.'+ext),dpi=135)
  plt.close(fig)
 return row

def check_limit(r):
 expected=[r.get(k) for k in ['expected_2p5','expected_16','expected_50','expected_84','expected_97p5']]
 r['valid_expected']=all(v is not None and np.isfinite(v) and v>0 for v in expected) and all(a<=b for a,b in zip(expected,expected[1:]))
 log=(R/r['directory']/'combine.log').read_text(errors='replace');matches=re.findall(r'At r =\s*([\d.eE+-]+):[^\n]*CLs\s*=\s*([\d.eE+-]+)',log);obs=r.get('observed');near=min(matches,key=lambda m:abs(float(m[0])-obs)) if matches and obs is not None else None
 r['cls_at_reported_limit']=float(near[1]) if near else None;r['cls_r_relative_distance']=abs(float(near[0])-obs)/max(abs(obs),1e-9) if near else None
 r['valid_observed']=obs is not None and np.isfinite(obs) and obs>0 and near is not None and abs(float(near[1])-.05)<=.01 and r['cls_r_relative_distance']<.02
 r['valid']=r['valid_expected'] and r['valid_observed'];return r

def aggregate():
 primary=sorted([dict(r) for r in cache.values()],key=lambda a:(a['family'],a['mass'],a['role']));
 for r in primary:
  if r['role']=='limits' and r.get('exit_code')==0:check_limit(r)
 csvwrite(R/'summary/primary_all_jobs.csv',primary)
 selected={(r['family'],r['mass'],r['role']):r for r in primary}
 for rp in sorted((R/'limit_repairs').glob('*/*/receipt.json')):
  if json.loads(rp.read_text())['exit_code']!=0:continue
  r=check_limit(process(rp))
  if r['valid']:r['replaces_primary_directory']=selected[r['family'],r['mass'],r['role']]['directory'];r['numerically_recalculated']=True;selected[r['family'],r['mass'],r['role']]=r
 rows=sorted(selected.values(),key=lambda a:(a['family'],a['mass'],a['role']));
 for r in rows:
  if r['role']=='significance' and r.get('valid'):r['log10_p0']=float(norm.logsf(r['Z'])/np.log(10))
 csvwrite(R/'summary/all_jobs.csv',rows)
 for role in ['limits','significance','fitdiagnostics']:csvwrite(R/f'summary/{role}.csv',[r for r in rows if r['role']==role])
 baseline=P/'runs/20261029_original_mass_scan/summary';origlim=list(csv.DictReader((baseline/'limits.csv').open()));origsig=list(csv.DictReader((baseline/'significance_pvalues.csv').open()))
 for fam in LABEL:
  a=[r for r in rows if r['family']==fam and r['role']=='limits' and r.get('valid_expected')]
  if a:
   fig,ax=plt.subplots(figsize=(9,5.5),layout='constrained');good=[r for r in a if r.get('expected_ordered')];x=[r['mass'] for r in good]
   ax.fill_between(x,[r['expected_2p5'] for r in good],[r['expected_97p5'] for r in good],color='#f8cf32',label='Atteso ±2σ');ax.fill_between(x,[r['expected_16'] for r in good],[r['expected_84'] for r in good],color='#66bd63',label='Atteso ±1σ');ax.plot([r['mass'] for r in a],[r['expected_50'] for r in a],'k--',label='Atteso mediano');ax.plot([r['mass'] for r in a],[r['observed'] if r.get('valid_observed') else np.nan for r in a],'ko-',ms=2.5,label='Osservato');ax.set_yscale('log');ax.set_xlabel('Massa ipotizzata [GeV]');ax.set_ylabel('Limite superiore 95% CL su r');ax.set_title(f'{LABEL[fam]} — {sum(r.get("valid_observed",False) for r in a)}/51 osservati verificati; {len(a)}/51 attesi\nStudio esplorativo, range 6–12 GeV, Υ libere');ax.grid(alpha=.2,which='both');ax.legend(fontsize=8)
   for ext in ['png','pdf']:fig.savefig(R/f'plots/limits_{fam}.{ext}',dpi=170)
   plt.close(fig)
 for quantity,title,ylabel in [('observed','Limiti osservati','Limite superiore 95% CL su r'),('expected_50','Limiti attesi mediani','Limite superiore 95% CL su r'),('p0','p-value locali osservati','p₀ locale')]:
  fig,ax=plt.subplots(figsize=(10,6),layout='constrained');role='significance' if quantity=='p0' else 'limits';anydata=False
  for fam in LABEL:
   a=[r for r in rows if r['family']==fam and r['role']==role and r.get('valid_expected' if quantity=='expected_50' else 'valid')]
   if a:anydata=True;ax.plot([round(6+i*.1,1) for i in range(51)],[max({round(r['mass'],1):r[quantity] for r in a}.get(round(6+i*.1,1),np.nan),1e-300) for i in range(51)],marker='o',ms=2.2,lw=1.2,ls=':' if fam=='doubleexp' else '-',color=COLORS[fam],label=f'{LABEL[fam]} ({len(a)}/51)'+(' — degenerata' if fam=='doubleexp' else ''))
  if quantity=='p0':
   ax.plot([float(r['mass']) for r in origsig],[float(r['p0_local']) for r in origsig],'k--',lw=1.3,label='Envelope originale (scan precedente)')
   pvals=[r['p0'] for r in rows if r['role']=='significance' and r.get('valid')];bottom=max(min([max(p,1e-300) for p in pvals]+[1e-8])/5,1e-300);ax.set_ylim(bottom=bottom,top=1)
  else:ax.plot([float(r['mass']) for r in origlim],[float(r[quantity]) for r in origlim],'k--',lw=1.3,label='Envelope originale (scan precedente)')
  ax.set_yscale('log');ax.set_xlim(6,11.15);ax.set_xlabel('Massa ipotizzata [GeV]');ax.set_ylabel(ylabel);ax.set_title(title+' — confronto delle famiglie\nFit combinato 6–12 GeV; risultati esplorativi, senza validazione del bias');ax.grid(alpha=.15,which='both');ax.legend(fontsize=8,loc='best')
  for ext in ['png','pdf']:fig.savefig(R/f'plots/{"pvalues" if quantity=="p0" else "limits_"+quantity}_comparison.{ext}',dpi=170)
  plt.close(fig)
 fig,ax=plt.subplots(figsize=(10,6),layout='constrained');floor=1e-12
 for fam in LABEL:
  a=[r for r in rows if r['family']==fam and r['role']=='significance' and r.get('valid')];x=np.array([r['mass'] for r in a]);y=np.array([r['p0'] for r in a]);ax.plot(x,np.where(y>=floor,y,np.nan),marker='o',ms=2.5,lw=1.2,ls=':' if fam=='doubleexp' else '-',color=COLORS[fam],label=LABEL[fam]+(' — degenerata' if fam=='doubleexp' else ''));mask=y<floor;ax.scatter(x[mask],np.full(mask.sum(),floor*1.25),marker='v',s=14,color=COLORS[fam])
 ax.plot([float(r['mass']) for r in origsig],[float(r['p0_local']) for r in origsig],'k--',label='Envelope originale (precedente)')
 for z in [1,2,3,4,5]:ax.axhline(norm.sf(z),color='gray',lw=.7,ls=':');ax.text(11.02,norm.sf(z),f'{z}σ',fontsize=8,va='bottom')
 ax.set_yscale('log');ax.set_ylim(floor,1);ax.set_xlim(6,11.15);ax.set_xlabel('Massa ipotizzata [GeV]');ax.set_ylabel('p₀ locale');ax.set_title('p-value locali — dettaglio della scala superiore\nTriangoli: valori sotto 10⁻¹², fuori scala; risultati esplorativi');ax.legend(fontsize=8,loc='lower left');ax.grid(alpha=.15,which='both')
 for ext in ['png','pdf']:fig.savefig(R/f'plots/pvalues_comparison_zoom.{ext}',dpi=170)
 plt.close(fig)
 stats=dict(processed=len(rows),expected=612,execution_failures=sum(r['exit_code']!=0 for r in rows),postprocess_errors=sum(bool(r.get('postprocess_error')) for r in rows),limits_valid=sum(r['role']=='limits' and r.get('valid_observed',False) for r in rows),limits_expected_valid=sum(r['role']=='limits' and r.get('valid_expected',False) for r in rows),significance_valid=sum(r['role']=='significance' and r.get('valid',False) for r in rows),fitdiagnostics_valid=sum(r['role']=='fitdiagnostics' and r.get('valid',False) for r in rows),diagnostic_warnings=sum(r.get('diagnostic_warning',False) for r in rows));(R/'summary/status.json').write_text(json.dumps(stats,indent=2));(R/'summary/selected.json').write_text(json.dumps(rows,indent=2));cp.write_text(json.dumps(cache,indent=2));print(json.dumps(stats),flush=True)

iteration=0
while True:
 fresh=0
 for rp in sorted((R/'fits').glob('*/*/*/receipt.json')):
  key=str(rp.relative_to(R))
  if key in cache and cache[key].get('receipt_mtime')==rp.stat().st_mtime:continue
  try:r=process(rp)
  except Exception as e:
   plt.close('all');t=json.loads(rp.read_text());r={k:t[k] for k in ['family','mass','mass_label','role','exit_code']};r.update(valid=False,postprocess_error=str(e),receipt_mtime=rp.stat().st_mtime,directory=str(rp.parent.relative_to(R)));(rp.parent/'postprocess_error.log').write_text(traceback.format_exc())
  cache[key]=r;fresh+=1
 if fresh or iteration==0:aggregate()
 if not args.watch or len(cache)>=612:break
 iteration+=1;time.sleep(5)
