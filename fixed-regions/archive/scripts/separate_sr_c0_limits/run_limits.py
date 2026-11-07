"""Only native AsymptoticLimits for separate SR+CR analyses; resumable and locked."""
from pathlib import Path
import argparse,concurrent.futures,datetime,fcntl,json,math,os,signal,subprocess,time,traceback
P=Path(__file__).resolve().parents[2];S=P/'scripts/separate_sr_c0_limits';R=P/'runs/20261031_separate_sr_c0_limits'
PY=str(P/'bin/python');C=str(P/'bin/combine')
def now():return datetime.datetime.now().astimezone().isoformat()
def atomic(path,data):
 p=Path(path);p.parent.mkdir(parents=True,exist_ok=True);q=p.with_suffix(p.suffix+'.tmp.'+str(os.getpid()));q.write_text(json.dumps(data,indent=2));q.replace(p)
def scale(region,m):
 previous=P/f'runs/20261031_separate_sr_cheb_limits/{region}/m{m}/results.json'
 if previous.exists():
  a=json.loads(previous.read_text()); e=a.get('limits',{}).get('expected_50')
  if e and math.isfinite(e) and e>0:return float(e),str(previous)
 candidates=[P/f'runs/20261030_sr6_cr4_scan/m{m}/results.json'] if region=='low' else [P/f'runs/20261030_upper_envelope_scan/degree3/m{m}/results.json']
 candidates.append(P/f'runs/20261031_shared_sr_cr_scan/observed/m{m}/results.json')
 for p in candidates:
  if p.exists():
   e=json.loads(p.read_text()).get('limits',{}).get('expected_50')
   if e and math.isfinite(e) and e>0:return float(e),str(p)
 return 10.,'fallback numerical bracket scale; no inferred limit'
def job(point):
 region=point['region'];m=f"{point['mass']:.1f}";D=R/region/('m'+m);J=D/'limits';D.mkdir(parents=True,exist_ok=True)
 with (D/'job.lock').open('a') as lock:
  fcntl.flock(lock,fcntl.LOCK_EX)
  if (D/'completion.json').exists():return region,m,'already_closed'
  try:
   w=D/'workspace.root';deadline=time.time()+900
   while not(w.exists() and w.with_suffix('.audit.json').exists()):
    if time.time()>deadline:raise RuntimeError('Workspace not ready')
    time.sleep(2)
   E,source=scale(region,m);bound=min(1e6,max(1.,5*E));J.mkdir(parents=True,exist_ok=True)
   cmd=[C,'-M','AsymptoticLimits','-d',str(w),'-m',m,'-n','_separate','--rMin','0','--rMax',str(bound),'--rRelAcc','.001','--rAbsAcc','.00001','--cminDefaultMinimizerStrategy','0','--cminDefaultMinimizerTolerance','.001','--cminDefaultMinimizerPrecision','1e-12','--X-rtd','MINIMIZER_freezeDisassociatedParams','--setParameters','r=0','-v','1']
   atomic(D/'protocol.json',dict(region=region,mass=float(m),degree=point['degree'],method='AsymptoticLimits',initial_rMax=bound,bracket_scale=E,bracket_scale_source=source,no_discrete_envelope=True))
   receipt=J/'combine.receipt.json'
   if not receipt.exists():
    atomic(J/'combine.command.json',cmd);start=time.time()
    with (J/'combine.log').open('w') as log:
     proc=subprocess.Popen(cmd,cwd=J,stdout=log,stderr=subprocess.STDOUT,env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1'),start_new_session=True)
     atomic(J/'combine.process.json',dict(pid=proc.pid,started_at=now(),command=cmd))
     try:rc=proc.wait(timeout=3600)
     except subprocess.TimeoutExpired:
      os.killpg(proc.pid,signal.SIGTERM)
      try:proc.wait(timeout=5)
      except subprocess.TimeoutExpired:os.killpg(proc.pid,signal.SIGKILL);proc.wait()
      rc=124
    atomic(receipt,dict(exit_code=rc,seconds=time.time()-start,finished_at=now(),timeout_seconds=3600))
   a=subprocess.run([PY,str(S/'extract_limits.py'),region,m],cwd=P,capture_output=True,text=True,timeout=120)
   (D/'extraction.log').write_text(a.stdout+a.stderr)
   if a.returncode:raise RuntimeError('Extraction failed')
   result=json.loads((D/'results.json').read_text());atomic(D/'completion.json',dict(status=result['status'],flags=result['flags'],finished_at=now()))
   print(region,m,result['status'],result.get('limits',{}),result['flags'],flush=True);return region,m,result['status']
  except Exception as e:
   atomic(D/'worker_failure.json',dict(error=str(e),traceback=traceback.format_exc(),finished_at=now()));print('FAILED',region,m,str(e),flush=True);return region,m,'failed'
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--region',choices=['low','high']);ap.add_argument('--masses',nargs='+',type=float);ap.add_argument('--workers',type=int,default=8);ap.add_argument('--exclude-pilots',action='store_true');ap.add_argument('--name',default='runner');args=ap.parse_args();assert 1<=args.workers<=10
 manifest=json.loads((R/'manifest.json').read_text());points=manifest['points']
 if args.region:points=[a for a in points if a['region']==args.region]
 if args.masses:points=[a for a in points if a['mass'] in args.masses]
 if args.exclude_pilots:points=[a for a in points if (a['region'],a['mass']) not in [('low',7.5),('high',11.)]]
 points=sorted(points,key=lambda a:(0 if (a['region'],a['mass']) in [('low',7.5),('high',11.)] else 1,a['mass']))
 atomic(R/(args.name+'.json'),dict(status='running',pid=os.getpid(),started_at=now(),workers=args.workers,points=points))
 with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as pool:results=list(pool.map(job,points))
 atomic(R/(args.name+'.json'),dict(status='completed',finished_at=now(),workers=args.workers,results=results))
