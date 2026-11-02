from pathlib import Path
import subprocess,concurrent.futures,json,time,os
P=Path.cwd();R=P/'runs/20261030_upper_envelope_scan';S=P/'scripts/upper_envelope_scan';PY=str(P/'bin/python')
def job(args):
    degree,mass=args;D=R/f'bias/degree{degree}/m{mass}';D.mkdir(parents=True,exist_ok=True)
    source=R/f'degree{degree}/m{mass}/results.json'
    while not source.exists():time.sleep(5)
    a=json.loads(source.read_text());median=a['limits']['expected_50']
    for label,value in [('r0',0),('rExpected',median),('r10Expected',10*median)]:
        wd=D/label;wd.mkdir(exist_ok=True)
        if (wd/f'm{mass}/result.json').exists():continue
        for name,cmd in [('generate',[PY,str(S/'generator.py'),'bias',str(degree),mass,str(value),label]),('fit',[PY,str(S/'toy_worker.py'),'bias',str(degree),mass,label])]:
            with (wd/(name+'.log')).open('w') as h:rc=subprocess.run(cmd,cwd=P,stdout=h,stderr=subprocess.STDOUT,timeout=1800,env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1')).returncode
            (wd/(name+'.receipt.json')).write_text(json.dumps(dict(exit_code=rc)))
            if rc:print('BIAS FAIL',degree,mass,label,name,flush=True);return
        print('BIAS',degree,mass,label,flush=True)
with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:list(pool.map(job,[(d,f'{m/10:.1f}') for m in range(105,121) for d in [1,2]]))
print('BIAS DONE',flush=True)
