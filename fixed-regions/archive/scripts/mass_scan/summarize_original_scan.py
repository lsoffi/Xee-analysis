#!/usr/bin/env python3
import argparse,csv,json,math,pathlib,re,time,traceback
import numpy as np,uproot
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import norm
P=pathlib.Path('/Users/imb/Documents/Codex/Analisi-statistica-Xee');run=P/'runs/20261029_original_mass_scan';ap=argparse.ArgumentParser();ap.add_argument('--watch',action='store_true');args=ap.parse_args()
cache={};indexpath=run/'summary/processed.json'
if indexpath.exists():cache=json.loads(indexpath.read_text())
def writecsv(path,rows):
 if not rows:return
 keys=list(dict.fromkeys(k for r in rows for k in r));tmp=path.with_suffix(path.suffix+'.tmp')
 with tmp.open('w') as f:
  w=csv.DictWriter(f,fieldnames=keys);w.writeheader();w.writerows(rows)
 tmp.replace(path)
def valtree(tree,name):
 if name not in tree:return None
 a=tree[name].array(library='np');return float(a[0]) if len(a) else None
def process(receipt):
 t=json.loads(receipt.read_text());wd=receipt.parent;row={k:t[k] for k in ['mass','mass_label','role','era','exit_code','elapsed_seconds']};row['directory']=str(wd.relative_to(run));log=(wd/'combine.log').read_text(errors='replace') if (wd/'combine.log').exists() else '';row['warnings']=len(re.findall(r'\[WARNING\]|WARNING:',log));row['error']=t.get('error')
 if t['exit_code']!=0:row['valid']=False;return row
 if t['role'] in ['limits','significance']:
  roots=list(wd.glob('higgsCombine*.root'))
  if len(roots)!=1:raise RuntimeError('Missing or ambiguous result ROOT '+str(wd))
  with uproot.open(roots[0]) as f:
   values=f['limit']['limit'].array(library='np')
   if t['role']=='significance':
    if len(values)!=1:raise RuntimeError('Unexpected significance entries')
    z=float(values[0]);row.update(Z_local=z,p0_local=float(norm.sf(z)),log10_p0_local=float(norm.logsf(z)/np.log(10)),valid=bool(np.isfinite(z)))
   else:
    qs=f['limit']['quantileExpected'].array(library='np');targets={'observed':-1,'expected_2p5':.025,'expected_16':.16,'expected_50':.5,'expected_84':.84,'expected_97p5':.975}
    for k,q in targets.items():
     found=values[np.isclose(qs,q,atol=.004)];row[k]=float(found[0]) if len(found)==1 else None
    row['valid']=all(row[k] is not None and np.isfinite(row[k]) and row[k]>0 for k in targets)
    expected=[row[k] for k in ['expected_2p5','expected_16','expected_50','expected_84','expected_97p5']];row['expected_ordered']=all(v is not None for v in expected) and all(a<=b for a,b in zip(expected,expected[1:]))
    row['log_fit_failures']=len(re.findall(r'failed to converge|Minimization failed|fit failed|Failed to find',log,re.I))
  return row
 path=next(wd.glob('fitDiagnostics*.root'));matches=re.findall(r'Fit (B-only|S\+B), status =\s*(-?\d+), numBadNLL =\s*(\d+), covariance quality =\s*(-?\d+)',log)
 for stage,status,bad,cov in matches:label='b' if stage=='B-only' else 's';row.update({label+'_status':int(status),label+'_bad_nll':int(bad),label+'_covQual':int(cov)})
 with uproot.open(path) as f:
  b=f.get('tree_fit_b');s=f.get('tree_fit_sb')
  if s is not None:
   for k in ['r','rErr','rLoErr','rHiErr','nll_min','nll_nll0']:row[k]=valtree(s,k)
  row['valid']=row.get('b_status')==0 and row.get('s_status')==0 and row.get('b_covQual')==3 and row.get('s_covQual')==3
  pars=[]
  for tag,tree in [('fit_b',b),('fit_s',s)]:
   if tree is not None:
    for k in tree.keys():
     if k.startswith(('CMS_','lumi','t','a')) or k=='r':pars.append(dict(stage=tag,parameter=k,value=valtree(tree,k)))
  writecsv(wd/'parameters.csv',pars)
  channels=[k.split(';')[0] for k,cl in f['shapes_fit_s'].classnames(recursive=False).items() if cl.startswith('TDirectory')]
  era_names={'ch1':'2022','ch2':'2022EE','ch3':'2023','ch4':'2023BPix'};fig,axs=plt.subplots(len(channels),2,figsize=(12,3.3*len(channels)),squeeze=False,layout='constrained');yields=[]
  for i,ch in enumerate(channels):
   for j,tag in enumerate(['fit_b','fit_s']):
    ax=axs[i,j];dr=f['shapes_'+tag+'/'+ch];h=dr['total'];values,edges=h.to_numpy();width=np.diff(edges);centers=(edges[:-1]+edges[1:])/2;g=dr['data'];data=g.values('y')*width
    if len(data)!=len(values):raise RuntimeError('Data/shape bin count mismatch')
    ax.errorbar(centers,data,yerr=np.sqrt(np.maximum(data,0)),fmt='.',ms=2,color='black',elinewidth=.4,label='Dati',zorder=5)
    table={'mass_GeV':centers,'low':edges[:-1],'high':edges[1:],'data':data}
    for key,col,label,ls in [('total','#c62828','Totale','-'),('dy','#1766ba','Continuo originale','--'),('upsilon1s','#ef8a00','Υ(1S), norma fissa','-.'),('upsilon2s','#9c27b0','Υ(2S)',':'),('Zd','#188341','Segnale','-')]:
     y=dr[key].values()*width;table[key]=y;yields.append(dict(channel=ch,stage=tag,component=key,events=float(y.sum())))
     if np.any(y<0):label='|'+label+'| (negativo)';draw=np.abs(y);ls='-.'
     else:draw=y
     if np.any(draw>0):ax.stairs(np.where(draw>0,draw,np.nan),edges,color=col,label=label,ls=ls,lw=1.3,baseline=None)
    writecsv(wd/(ch+'_'+tag+'_shapes.csv'),[{k:float(v[n]) for k,v in table.items()} for n in range(len(values))])
    ax.set_yscale('log');ax.set_xlim(6,12);ax.set_ylim(.01,max(float(np.max(data)),float(np.max(values*width)),1)*3);ax.set_title(era_names.get(ch,ch)+(' — B' if tag=='fit_b' else ' — S+B')+f"; status {row.get('b_status' if tag=='fit_b' else 's_status','?')}, covQual {row.get('b_covQual' if tag=='fit_b' else 's_covQual','?')}",fontsize=10);ax.set_ylabel('Eventi / bin');ax.set_xlabel('Massa [GeV]');ax.grid(alpha=.12,which='both');ax.legend(fontsize=6.5,ncol=2,loc='lower left')
  fig.suptitle(f"Workspace originale — massa {t['mass_label']} GeV, {t['era']}\nFit 6–12 GeV; Υ(1S) congelata solo in FitDiagnostics",fontsize=11)
  for ext in ['png','pdf']:fig.savefig(wd/('fits.'+ext),dpi=140)
  plt.close(fig);writecsv(wd/'normalizations_from_shapes.csv',yields)
 return row

def aggregate(final=False):
 rows=sorted(cache.values(),key=lambda r:(r['mass'],r['role'],r['era']));writecsv(run/'summary/all_jobs.csv',rows);writecsv(run/'summary/fit_diagnostics.csv',[r for r in rows if r['role']=='fitdiagnostics']);limits=[r for r in rows if r['role']=='limits'];sigs=[r for r in rows if r['role']=='significance'];writecsv(run/'summary/limits.csv',limits);writecsv(run/'summary/significance_pvalues.csv',sigs)
 good=[r for r in limits if r.get('valid')]
 if good:
  x=np.array([r['mass'] for r in good]);f,ax=plt.subplots(figsize=(9,5.5),layout='constrained');band=[r for r in good if r.get('expected_ordered')]
  if band:
   xb=[r['mass'] for r in band];ax.fill_between(xb,[r['expected_2p5'] for r in band],[r['expected_97p5'] for r in band],color='#f8cf32',label='Atteso ±2σ');ax.fill_between(xb,[r['expected_16'] for r in band],[r['expected_84'] for r in band],color='#66bd63',label='Atteso ±1σ')
  ax.plot(x,[r['expected_50'] for r in good],'k--',lw=1.5,label='Atteso mediano');ax.plot(x,[r['observed'] for r in good],'ko-',ms=3,lw=1.3,label='Osservato');ax.set_yscale('log');ax.set_xlabel('Massa ipotetica [GeV]');ax.set_ylabel('Limite superiore 95% CL su r');ax.grid(alpha=.2,which='both');ax.legend(fontsize=9);ax.set_title(f'AsymptoticLimits — workspace originali, tutte le ere\n{len(good)}/51 punti; range di fit 6–12 GeV; Υ(1S) libera')
  for ext in ['png','pdf']:f.savefig(run/'plots'/('limits_vs_mass.'+ext),dpi=180)
  plt.close(f)
 good=[r for r in sigs if r.get('valid')]
 if good:
  x=[r['mass'] for r in good];ys=[max(r['p0_local'],1e-300) for r in good];f,ax=plt.subplots(figsize=(9,5.5),layout='constrained');ax.plot(x,ys,'ko-',ms=3,lw=1.2,label='p₀ locale osservato');bottom=max(min(ys)/5,1e-300)
  for z in [1,2,3,4,5]:
   pv=norm.sf(z)
   if pv>=bottom:ax.axhline(pv,ls='--',color='#999999',lw=.8);ax.text(max(x),pv,f' {z}σ',fontsize=8,va='bottom',ha='left')
  ax.set_yscale('log');ax.set_ylim(bottom,1);ax.set_xlabel('Massa ipotetica [GeV]');ax.set_ylabel('p₀ locale');ax.grid(alpha=.15,which='both');ax.set_title(f'Significance — workspace originali, tutte le ere\n{len(good)}/51 punti; asintotica locale, senza correzione look-elsewhere');ax.legend(fontsize=9)
  for ext in ['png','pdf']:f.savefig(run/'plots'/('pvalue_vs_mass.'+ext),dpi=180)
  plt.close(f)
  f,ax=plt.subplots(figsize=(9,5),layout='constrained');ax.plot(x,[r['Z_local'] for r in good],'ko-',ms=3);ax.set_xlabel('Massa ipotetica [GeV]');ax.set_ylabel('Significanza locale [σ]');ax.grid(alpha=.2);ax.set_title('Significanza osservata — originali, Υ(1S) libera');f.savefig(run/'plots/significance_vs_mass.png',dpi=180);f.savefig(run/'plots/significance_vs_mass.pdf');plt.close(f)
 stats=dict(processed=len(rows),expected=357,execution_failures=sum(r['exit_code']!=0 for r in rows),diagnostic_quality_warnings=sum(r['role']=='fitdiagnostics' and not r.get('valid',False) for r in rows),limit_valid=sum(r.get('valid',False) for r in limits),significance_valid=sum(r.get('valid',False) for r in sigs),postprocess_errors=sum(bool(r.get('postprocess_error')) for r in rows))
 (run/'summary/status.json').write_text(json.dumps(stats,indent=2));(run/'summary/processed.json').write_text(json.dumps(cache,indent=2));return stats
iteration=0
while True:
 fresh=0
 for receipt in sorted((run/'fits').rglob('receipt.json')):
  key=str(receipt.relative_to(run))
  if key in cache:continue
  try:row=process(receipt)
  except Exception as e:
   plt.close('all');t=json.loads(receipt.read_text());row={k:t[k] for k in ['mass','mass_label','role','era','exit_code','elapsed_seconds']};row.update(valid=False,postprocess_error=str(e),directory=str(receipt.parent.relative_to(run)));(receipt.parent/'postprocess_error.log').write_text(traceback.format_exc())
  cache[key]=row;fresh+=1
 if fresh or iteration==0:
  stats=aggregate();print(json.dumps(stats),flush=True)
 if not args.watch or len(cache)>=357:break
 iteration+=1;time.sleep(5)
