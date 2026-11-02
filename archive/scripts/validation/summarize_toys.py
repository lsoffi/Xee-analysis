import pathlib,json,csv,re,time,argparse
import numpy as np,uproot
from scipy.stats import norm
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
P=pathlib.Path('/Users/imb/Documents/Codex/Analisi-statistica-Xee/runs/20261029_background_validation');ap=argparse.ArgumentParser();ap.add_argument('--watch',action='store_true');args=ap.parse_args();cache={}
def writecsv(path,rows):
 if not rows:return
 with path.open('w') as f:
  w=csv.DictWriter(f,fieldnames=list(dict.fromkeys(k for r in rows for k in r)));w.writeheader();w.writerows(rows)
def wilson(k,n):
 if not n:return (None,None)
 v=k/n;d=1+1/n;c=(v+.5/n)/d;h=np.sqrt(v*(1-v)/n+.25/n**2)/d;return c-h,c+h
while True:
 for receipt in P.glob('toys/m*/*/*/batch_*/*/receipt.json'):
  wd=receipt.parent
  if wd.name not in ['original','degree5'] or str(wd) in cache:continue
  meta=json.loads((wd.parent/'scenario.json').read_text());rr=json.loads(receipt.read_text());rows=[]
  if rr['exit_code']:
   cache[str(wd)]=[dict(**meta,model=wd.name,execution_failed=True)];continue
  with uproot.open(wd/'fitDiagnostics_toys.root') as f:
   b=f['tree_fit_b'].arrays(library='np');s=f['tree_fit_sb'].arrays(library='np');n=len(s['r'])
   log=(wd/'combine.log').read_text(errors='replace');matches=re.findall(r'Fit (B-only|S\+B), status =\s*(-?\d+), numBadNLL =\s*(\d+), covariance quality =\s*(-?\d+)',log);bc=[int(c) for st,_,_,c in matches if st=='B-only'];sc=[int(c) for st,_,_,c in matches if st=='S+B']
   for i in range(n):
    mu=float(s['r'][i]);lo=abs(float(s['rLoErr'][i]));hi=abs(float(s['rHiErr'][i]));d=float(s['nll_nll0'][i]);truth=meta['injection'];covb=bc[i] if len(bc)==n else -1;covs=sc[i] if len(sc)==n else -1;bounded=abs(mu)>99.99 or mu-lo<=-99.99 or mu+hi>=99.99
    valid=int(b['fit_status'][i])==0 and int(s['fit_status'][i])==0 and covb==3 and covs==3 and np.isfinite(mu+lo+hi+d) and lo>0 and hi>0 and not bounded and d<=.01
    denom=lo if mu>=truth else hi
    rows.append(dict(mass=meta['mass'],generator=meta['generator'],injection=truth,seed=meta['seed'],toy=i+1,model=wd.name,r=mu,rErr=float(s['rErr'][i]),low=mu-lo,high=mu+hi,pull=(mu-truth)/denom if denom>0 else None,covered68=bool(mu-lo<=truth<=mu+hi),Z=np.sqrt(max(0,-2*d)) if mu>0 else 0,b_status=int(b['fit_status'][i]),s_status=int(s['fit_status'][i]),b_covQual=covb,s_covQual=covs,interval_at_bound=bounded,delta_nll=d,valid=bool(valid),valid_point=bool(int(b['fit_status'][i])==0 and int(s['fit_status'][i])==0 and np.isfinite(mu+d) and abs(mu)<99.99 and d<=.01)))
  cache[str(wd)]=rows
 rows=[r for group in cache.values() for r in group];writecsv(P/'toys_all.csv',rows);summ=[]
 for m in [7.5,9.5]:
  for g in ['cheb7','bernstein']:
   for inj in ([0,.06] if m==7.5 else [0,3.,25.]):
    for model in ['original','degree5']:
     allr=[r for r in rows if r.get('mass')==m and r.get('generator')==g and r.get('injection')==inj and r.get('model')==model];v=[r for r in allr if r.get('valid')];n=len(v)
     if not allr:continue
     point=[r for r in allr if r.get('valid_point')];npnt=len(point)
     a=dict(mass=m,generator=g,injection=inj,model=model,n_total=len(allr),n_valid=n,n_excluded=len(allr)-n,n_point_valid=npnt,point_fraction_Z_ge_2=sum(r['Z']>=2 for r in point)/npnt if npnt else None,point_fraction_Z_ge_3=sum(r['Z']>=3 for r in point)/npnt if npnt else None)
     if n:
      rh=np.array([r['r'] for r in v]);pull=np.array([r['pull'] for r in v]);covered=sum(r['covered68'] for r in v);clo,chi=wilson(covered,n)
      a.update(mean_r=float(rh.mean()),sd_r=float(rh.std(ddof=1)) if n>1 else None,bias=float(rh.mean()-inj),bias_mc_error=float(rh.std(ddof=1)/np.sqrt(n)) if n>1 else None,mean_pull=float(pull.mean()),sd_pull=float(pull.std(ddof=1)) if n>1 else None,coverage68=covered/n,coverage68_lo=clo,coverage68_hi=chi,median_rErr=float(np.median([r['rErr'] for r in v])),fraction_Z_ge_1=sum(r['Z']>=1 for r in v)/n,fraction_Z_ge_2=sum(r['Z']>=2 for r in v)/n,fraction_Z_ge_3=sum(r['Z']>=3 for r in v)/n)
     summ.append(a)
 writecsv(P/'toy_summary.csv',summ);(P/'toy_summary.json').write_text(json.dumps(summ,indent=2));(P/'postprocess_status.json').write_text(json.dumps(dict(completed_fit_batches=len(cache),expected_fit_batches=200,entries=len(rows),valid=sum(r.get('valid',False) for r in rows)),indent=2))
 if summ:
  for injname in ['zero','signal']:
   fig,axs=plt.subplots(2,2,figsize=(12,8),layout='constrained')
   for ax,(m,g) in zip(axs.flat,[(m,g) for m in [7.5,9.5] for g in ['cheb7','bernstein']]):
    inj=0 if injname=='zero' else (.06 if m==7.5 else 3.)
    for model,col in [('original','#246eae'),('degree5','#d35226')]:
     v=[r['pull'] for r in rows if r.get('valid') and r['mass']==m and r['generator']==g and r['injection']==inj and r['model']==model]
     if v:ax.hist(v,bins=np.linspace(-6,6,25),density=True,histtype='step',lw=1.5,color=col,label=f'{model}, n={len(v)}')
    xx=np.linspace(-6,6,300);ax.plot(xx,norm.pdf(xx),'k--',lw=1,label='N(0,1)');ax.set_xlim(-6,6);ax.set_title(f'm={m}, generatore {g}, r iniettato={inj}');ax.set_xlabel('(r̂ − r vero) / errore verso il valore vero');ax.set_ylabel('Densità');ax.legend(fontsize=8);ax.grid(alpha=.15)
   fig.suptitle('Test di bias e intervalli — '+('nessun segnale' if injname=='zero' else 'signal injection'))
   for ext in ['png','pdf']:fig.savefig(P/f'pulls_{injname}.{ext}',dpi=170)
   plt.close(fig)
 print(json.dumps(dict(batches=len(cache),entries=len(rows))),flush=True)
 if not args.watch or len(cache)>=200:break
 time.sleep(20)
