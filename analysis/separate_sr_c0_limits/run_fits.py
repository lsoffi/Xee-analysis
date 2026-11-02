"""Native FitDiagnostics for the requested C0 model; independent from limits."""
from run_limits import P,S,R,PY,C,atomic,now,scale
from pathlib import Path
import argparse,concurrent.futures,fcntl,json,math,os,signal,subprocess,time,traceback

def execute(cmd,D,timeout=1800):
 D.mkdir(parents=True,exist_ok=True);receipt=D/'combine.receipt.json'
 if receipt.exists():return json.loads(receipt.read_text())
 atomic(D/'combine.command.json',cmd);start=time.time()
 with (D/'combine.log').open('w') as log:
  proc=subprocess.Popen(cmd,cwd=D,stdout=log,stderr=subprocess.STDOUT,env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1'),start_new_session=True)
  atomic(D/'combine.process.json',dict(pid=proc.pid,started_at=now(),command=cmd))
  try:rc=proc.wait(timeout=timeout)
  except subprocess.TimeoutExpired:
   os.killpg(proc.pid,signal.SIGTERM)
   try:proc.wait(timeout=5)
   except subprocess.TimeoutExpired:os.killpg(proc.pid,signal.SIGKILL);proc.wait()
   rc=124
 a=dict(exit_code=rc,seconds=time.time()-start,finished_at=now(),timeout_seconds=timeout);atomic(receipt,a);return a

def job(point):
 side=point['region'];m=f"{point['mass']:.1f}";D=R/side/('m'+m);D.mkdir(parents=True,exist_ok=True)
 with (D/'fits.lock').open('a') as lock:
  fcntl.flock(lock,fcntl.LOCK_EX)
  if (D/'fits_completion.json').exists():return side,m,'already_closed'
  try:
   w=D/'workspace.root';deadline=time.time()+900
   while not(w.exists() and w.with_suffix('.audit.json').exists()):
    if time.time()>deadline:raise RuntimeError('Workspace not ready')
    time.sleep(2)
   numerical=json.loads((P/'configs/numerical_scales.json').read_text())[f'{side}/m{m}']
   E=numerical['fit_expected_anchor'];src=str(P/'configs/numerical_scales.json')
   bound=max(100.,5*E)
   atomic(D/'fit_protocol.json',dict(method='FitDiagnostics',side=side,mass=float(m),strategy=2,tolerance=.001,r_range=[-bound,bound],anchors=[0,E],expected_seed_source=src,robustFit=0,minos='none'))
   for label,start in [('zero',0.),('expected',E)]:
    cmd=[C,'-M','FitDiagnostics','-d',str(w),'-m',m,'-n','_c0','--cminDefaultMinimizerStrategy','2','--cminDefaultMinimizerTolerance','.001','--cminDefaultMinimizerPrecision','1e-12','--X-rtd','MINIMIZER_freezeDisassociatedParams','--saveShapes','--saveNormalizations','--keepFailures','--robustFit','0','--minos','none','--rMin',str(-bound),'--rMax',str(bound),'--preFitValue',str(start),'-v','1']
    if label=='expected':cmd+=['--skipBOnlyFit']
    execute(cmd,D/('fitdiagnostics_anchor_'+label))
   parser=S/'extract_fits.py';deadline=time.time()+300
   while not (parser.exists() and (S/'fit_extractor.ready').exists()):
    if time.time()>deadline:raise RuntimeError('Fit extractor not ready')
    time.sleep(2)
   a=subprocess.run([PY,str(parser),side,m],cwd=P,capture_output=True,text=True,timeout=240)
   (D/'fit_extraction.log').write_text(a.stdout+a.stderr)
   if a.returncode:raise RuntimeError('Fit extraction failed')
   result=json.loads((D/'fit_results.json').read_text());atomic(D/'fits_completion.json',dict(status=result.get('status','complete'),flags=result.get('flags',[]),finished_at=now()))
   print(side,m,'fits_complete',result.get('flags',[]),flush=True);return side,m,'fits_complete'
  except Exception as e:
   atomic(D/'fit_worker_failure.json',dict(error=str(e),traceback=traceback.format_exc(),finished_at=now()));print('FIT FAILED',side,m,str(e),flush=True);return side,m,'failed'

if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--masses',nargs='+',type=float);ap.add_argument('--workers',type=int,default=4);ap.add_argument('--exclude-pilots',action='store_true');ap.add_argument('--name',default='fit_runner');a=ap.parse_args();assert 1<=a.workers<=10
 points=json.loads((R/'manifest.json').read_text())['points']
 if a.masses:points=[v for v in points if v['mass'] in a.masses]
 if a.exclude_pilots:points=[v for v in points if (v['region'],v['mass']) not in [('low',7.5),('high',11.)]]
 points=sorted(points,key=lambda v:(0 if (v['region'],v['mass']) in [('low',7.5),('high',11.)] else 1,0 if v['mass']<9.1 or v['mass']>10.5 else 1,v['mass']))
 atomic(R/(a.name+'.json'),dict(status='running',pid=os.getpid(),started_at=now(),workers=a.workers,points=points))
 with concurrent.futures.ThreadPoolExecutor(max_workers=a.workers) as pool:out=list(pool.map(job,points))
 atomic(R/(a.name+'.json'),dict(status='completed',finished_at=now(),workers=a.workers,results=out))
