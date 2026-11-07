"""One native Significance attempt per C0 point, seeded from selected FD S+B."""
from pathlib import Path
import argparse,concurrent.futures,datetime,fcntl,hashlib,json,math,os,signal,subprocess,time,traceback

P=Path(__file__).resolve().parents[2];S=P/'scripts/separate_sr_c0_limits';R=P/'runs/20261031_separate_sr_c0_limits'
PY=str(P/'bin/python');C=str(P/'bin/combine')
def now():return datetime.datetime.now().astimezone().isoformat()
def atomic(path,obj):
    path=Path(path);temp=path.with_suffix(path.suffix+'.tmp.'+str(os.getpid()));temp.write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n');temp.replace(path)
def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def job(point):
    side=point['region'];mass=float(point['mass']);D=R/side/f'm{mass:.1f}';J=D/'significance';J.mkdir(exist_ok=True)
    with (D/'significance.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX)
        if (D/'significance_results.json').exists():return dict(side=side,mass=mass,status='already_closed')
        try:
            workspace=D/'workspace.root';audit=json.loads(workspace.with_suffix('.audit.json').read_text());assert sha(workspace)==audit['workspace_sha256']
            fd=json.loads((D/'fit_results.json').read_text());protocol=json.loads((D/'fit_protocol.json').read_text());bounds=protocol['r_range']
            initial={v['name']:v['value'] for v in fd['fit']['s']['parameters']};initial['r']=max(0.,fd['fit']['s']['r'])
            assert all(math.isfinite(v) for v in initial.values())
            params=','.join(f'{name}={value:.17g}' for name,value in sorted(initial.items()))
            cmd=[C,'-M','Significance','-d',str(workspace),'-m',f'{mass:.1f}','-n','_c0','--redefineSignalPOIs','r',
                 '--rMin',str(bounds[0]),'--rMax',str(bounds[1]),'--setParameters',params,
                 '--cminDefaultMinimizerStrategy','2','--cminDefaultMinimizerTolerance','0.00001','--cminDefaultMinimizerPrecision','1e-12',
                 '--X-rtd','MINIMIZER_freezeDisassociatedParams','-v','1']
            receipt=J/'combine.receipt.json'
            if not receipt.exists():
                if (J/'combine.process.json').exists():raise RuntimeError('Preserving an earlier unfinished native attempt; no automatic retry')
                atomic(J/'protocol.json',dict(method='Significance',side=side,mass=mass,workspace_sha256=audit['workspace_sha256'],
                    source_FD=fd['fit']['s']['source'],source_FD_sha256=sha(fd['fit']['s']['source']),initial_values=initial,
                    r_range=bounds,strategy=2,tolerance=1e-5,precision=1e-12,parameter_freezes=[],attempts=1))
                atomic(J/'combine.command.json',cmd);started=time.time()
                with (J/'combine.log').open('w') as log:
                    proc=subprocess.Popen(cmd,cwd=J,stdout=log,stderr=subprocess.STDOUT,
                        env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1'),start_new_session=True)
                    atomic(J/'combine.process.json',dict(pid=proc.pid,started_at=now(),command=cmd))
                    try:rc=proc.wait(timeout=1800)
                    except subprocess.TimeoutExpired:
                        os.killpg(proc.pid,signal.SIGTERM)
                        try:proc.wait(timeout=5)
                        except subprocess.TimeoutExpired:os.killpg(proc.pid,signal.SIGKILL);proc.wait()
                        rc=124
                atomic(receipt,dict(exit_code=rc,seconds=time.time()-started,finished_at=now(),timeout_seconds=1800))
            extraction=subprocess.run([PY,str(S/'extract_significance.py'),side,f'{mass:.1f}'],cwd=P,capture_output=True,text=True,timeout=120)
            (J/'extraction.log').write_text(extraction.stdout+extraction.stderr)
            if extraction.returncode:raise RuntimeError('Native significance extraction failed: '+extraction.stderr[-1500:])
            assert sha(workspace)==audit['workspace_sha256']
            result=json.loads((D/'significance_results.json').read_text())
            print(side,f'{mass:.1f}',result['status'],result['Z'],result['p0'],result['flags'],flush=True)
            return dict(side=side,mass=mass,status=result['status'],Z=result['Z'],p0=result['p0'],flags=result['flags'])
        except Exception as error:
            atomic(J/'worker_failure.json',dict(error=str(error),traceback=traceback.format_exc(),finished_at=now()))
            print('FAILED',side,mass,str(error),flush=True);return dict(side=side,mass=mass,status='failed',error=str(error))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--workers',type=int,default=4);p.add_argument('--pilots',action='store_true');p.add_argument('--name',default='significance_runner');a=p.parse_args()
    assert 1<=a.workers<=4
    points=json.loads((R/'manifest.json').read_text())['points'];pilots={('low',7.5),('high',11.)}
    if a.pilots:points=[q for q in points if (q['region'],q['mass']) in pilots]
    points=sorted(points,key=lambda q:(0 if (q['region'],q['mass']) in pilots else 1,q['mass'],q['region']))
    atomic(R/(a.name+'.json'),dict(status='running',pid=os.getpid(),started_at=now(),workers=a.workers,points=points))
    with concurrent.futures.ThreadPoolExecutor(max_workers=a.workers) as pool:out=list(pool.map(job,points))
    atomic(R/(a.name+'.json'),dict(status='completed',finished_at=now(),workers=a.workers,results=out))
