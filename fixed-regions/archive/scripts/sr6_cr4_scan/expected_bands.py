from pathlib import Path
import subprocess,json,re,time,sys,concurrent.futures
import numpy as np,uproot
from scipy.stats import norm
from scipy.optimize import brentq
P=Path.cwd();R=P/'runs/20261030_sr6_cr4_scan';C=str(P/'bin/combine');TARGET={q:float((norm.ppf(q)+norm.isf(.05*q))**2) for q in [.84,.975]}
def job(mass):
 d=R/('m'+mass)/'joint';dest=d/'expected_profile';dest.mkdir(exist_ok=True)
 if (dest/'result.json').exists():return
 start=time.time();original=json.loads((d.parent/'results.json').read_text());basecmd=json.loads((d/'limits/command.json').read_text());cache={0.:0.};evaluations=[]
 # Use previous profiled q_A points as initial brackets only. Final crossings are checked with fresh singlePoint evaluations.
 for r,q in re.findall(r'At r =\s*([\d.eE+-]+):[^\n]*q_A\s*=\s*([\d.eE+-]+)',(d/'limits/combine.log').read_text()):cache[float(r)]=float(q)
 def evaluate(x,force=False):
  x=float(x)
  if not force and x in cache:return cache[x]
  wd=dest/('r_'+format(x,'.10g').replace('.','p'));wd.mkdir(exist_ok=True);valuefile=wd/'value.json'
  if valuefile.exists():v=json.loads(valuefile.read_text());cache[x]=v['qA'];return v['qA']
  cmd=basecmd+['--singlePoint',format(x,'.14g')];cmd[cmd.index('--rMax')+1]=str(max(float(cmd[cmd.index('--rMax')+1]),2*x));cmd[cmd.index('--cminDefaultMinimizerStrategy')+1]='1';cmd[cmd.index('--cminDefaultMinimizerTolerance')+1]='0.001';(wd/'command.json').write_text(json.dumps(cmd));tt=time.time()
  with (wd/'combine.log').open('w') as out:res=subprocess.run(cmd,cwd=wd,stdout=out,stderr=subprocess.STDOUT,timeout=600)
  files=list(wd.glob('higgsCombine*.root'));assert res.returncode==0 and files,(mass,x,'call failed');a=uproot.open(files[0])['limit'].arrays(library='np');cls=float(a['limit'][np.argmin(abs(a['quantileExpected']-.5))]);qa=float(norm.isf(cls/2)**2);assert np.isfinite(qa) and cls>0,(mass,x,cls)
  v=dict(r=x,qA=qa,cls_median=cls,seconds=time.time()-tt,exit_code=res.returncode);valuefile.write_text(json.dumps(v,indent=2));cache[x]=qa;evaluations.append(v);print('POINT',mass,x,qa,flush=True);return qa
 result={}
 try:
  for quantile in [.84,.975]:
   target=TARGET[quantile];ordered=sorted(cache);below=[x for x in ordered if cache[x]<=target];lo=max(below);above=[x for x in ordered if x>lo and cache[x]>=target]
   if above:hi=min(above)
   else:
    hi=max(max(ordered)*1.5,original['limits']['expected_50']*3,lo*1.5)
    for k in range(16):
     if evaluate(hi)>=target:break
     hi*=1.7
    else:raise RuntimeError('Upper crossing not bracketed')
   # Brent solver profiles ALL discrete combinations through native Combine at new points.
   root=brentq(lambda x:evaluate(x)-target,lo,hi,xtol=1e-5,rtol=.001,maxiter=25);qa=evaluate(root,True);cls=float(norm.sf(np.sqrt(qa)-norm.ppf(quantile))/quantile)
   if abs(cls-.05)>.0005:
    ordered=sorted(cache);lo=max(x for x in ordered if cache[x]<target);hi=min(x for x in ordered if x>lo and cache[x]>target);root=brentq(lambda x:evaluate(x)-target,lo,hi,xtol=1e-6,rtol=.0002,maxiter=20);qa=evaluate(root,True);cls=float(norm.sf(np.sqrt(qa)-norm.ppf(quantile))/quantile)
   assert abs(cls-.05)<.0005,(mass,quantile,cls);result[str(quantile)]=dict(limit=root,qA=qa,target_qA=target,CLs=cls,verified=True)
  out=dict(mass=float(mass),quantiles=result,seconds=time.time()-start,new_evaluations=len(evaluations),method='Native Combine AsymptoticLimits singlePoint, complete discrete profiling; adaptive bracket and verified CLs=.05 for each upper expected quantile')
  (dest/'result.json').write_text(json.dumps(out,indent=2));print('DONE BANDS',mass,out,flush=True)
 except Exception as e:(dest/'failure.json').write_text(json.dumps(dict(error=str(e))));print('BANDS FAIL',mass,str(e),flush=True)
if __name__=='__main__':
 masses=sys.argv[1:] or ['7.5']+[f'{k/10:.1f}' for k in range(60,91) if k!=75]
 with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:list(pool.map(job,masses))
 print('EXPECTED BANDS END',flush=True)
