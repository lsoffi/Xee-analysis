from pathlib import Path
import subprocess,concurrent.futures,json,time,os
P=Path.cwd();R=P/'runs/20261030_upper_envelope_scan';S=P/'scripts/upper_envelope_scan';PY=str(P/'bin/python')
def job(args):
    degree,toy=args;D=R/f'gof/degree{degree}/toy{toy:03d}';D.mkdir(parents=True,exist_ok=True)
    with (D/'worker.log').open('w') as h:rc=subprocess.run([PY,str(S/'toy_worker.py'),'gof',str(degree),str(toy)],cwd=P,stdout=h,stderr=subprocess.STDOUT,timeout=7200,env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1')).returncode
    (D/'worker_receipt.json').write_text(json.dumps(dict(exit_code=rc)));print('GOF',degree,toy,rc,flush=True)
    return rc
with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:
    for rc in pool.map(job,[(d,k) for k in range(51,101) for d in [1,2]]):pass
print('GOF DONE',flush=True)
