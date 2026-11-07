"""Resumable, CPU-bounded scheduler for shared-SR bias and GoF campaigns."""
from pathlib import Path
import argparse, concurrent.futures as cf, json, math, os, signal, subprocess, time, fcntl

P = Path(__file__).resolve().parents[2]
R = P / 'runs/20261031_shared_sr_cr_scan'
S = P / 'scripts/shared_sr_cr_scan'
PY = str(P / 'bin/python')
ENV = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', VECLIB_MAXIMUM_THREADS='1', PYTHONDONTWRITEBYTECODE='1')


def atomic(path, obj):
    path = Path(path);path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix+'.tmp');tmp.write_text(json.dumps(obj,indent=2));tmp.replace(path)


def subprocess_job(command, directory, timeout):
    directory.mkdir(parents=True, exist_ok=True)
    start=time.time();atomic(directory/'runner_command.json',command)
    with (directory/'runner.log').open('w') as log:
        proc=subprocess.Popen(command,cwd=P,env=ENV,stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
        try:code=proc.wait(timeout=timeout)
        except subprocess.TimeoutExpired:
            os.killpg(proc.pid,signal.SIGTERM)
            try:proc.wait(timeout=5)
            except subprocess.TimeoutExpired:os.killpg(proc.pid,signal.SIGKILL);proc.wait()
            code=124
    result=dict(exit_code=code,seconds=time.time()-start)
    atomic(directory/'runner_receipt.json',result)
    return result


def result_path(task):
    kind,mass,label,toy=task
    return (R/f'bias/fits/m{mass}/{label}/toy{toy:03d}/result.json' if kind=='bias'
            else R/f'gof/toy{toy:03d}/m{mass}/result.json')


def externally_locked(task):
    path=result_path(task).parent/'worker.lock'
    if not path.exists():return False
    with path.open('a') as lock:
        try:fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        except BlockingIOError:return True
        fcntl.flock(lock,fcntl.LOCK_UN)
    return False


def worker_budget(args):
    if args.worker_budget_file and args.worker_budget_file.exists():
        try:return max(1,min(8,int(json.loads(args.worker_budget_file.read_text())['workers'])))
        except (OSError,ValueError,KeyError,TypeError):pass
    return args.workers


def run_task(task):
    kind,mass,label,toy=task
    out=result_path(task).parent
    if result_path(task).exists():return dict(task=task,existing=True)
    if kind=='bias':
        cmd=[PY,str(S/'toys_worker.py'),'--mass',mass,'--injection',label,'--toy',str(toy),'--point-only']
        timeout=4800 if toy==1 else 3600
    else:
        cmd=[PY,str(S/'toys_gof_worker.py'),'--mass',mass,'--toy',str(toy)]
        timeout=2400
    value=subprocess_job(cmd,out,timeout)
    value.update(task=task,result_exists=result_path(task).exists())
    if not value['result_exists']:
        failure=dict(mass=float(mass),toy=toy,label=label,injection=label,
                fit_valid=False,valid=False,exception='runner_failed_before_result',runner=value,
                failure_flags=['worker_exception'],protocol='shared_sr_two_fixed_anchors_v1' if kind=='bias' else 'shared_sr_gof_rfree_v1')
        if kind=='gof':
            failure.update(score=None,nominal_valid=False,nominal_diagnostics=dict(available=False,attempted=False,outcome_complete=True,valid=False,flags=['GoF_worker_failed_before_native_result'],runner_exit_code=value['exit_code']))
        atomic(result_path(task),failure)
    return value


def median_ready(mass):
    path=R/f'observed/m{mass}/results.json'
    try:
        result=json.loads(path.read_text());x=float(result['limits']['expected_50'])
        return math.isfinite(x) and x>0
    except (OSError,ValueError,KeyError,TypeError):return False


def generation_ready(kind,mass,label):
    path=R/f'bias/generation/m{mass}/{label}/truth.json' if kind=='bias' else R/'gof/generation/truth.json'
    return path.exists()


def generate(kind,mass,label):
    if kind=='bias':
        cmd=[PY,str(S/'generator_toys.py'),'--mass',mass,'--label',label]
        out=R/f'bias/generation/m{mass}/{label}'
    else:
        cmd=[PY,str(S/'generator_toys.py'),'--gof'];out=R/'gof/generation'
    return subprocess_job(cmd,out,5400)


def main(args):
    masses=[f'{float(m):.1f}' for m in args.mass]
    # Fast point estimates first; control profiles are fixed toy1, never selected by outcome.
    toyorder=list(range(1,args.bias_toys+1))  # All toy IDs run point-only; profiles have a final queue.
    tasks=[]
    for toy in range(1,max(args.bias_toys,args.gof_toys)+1):
        for mass in masses:
            if args.mode in ['all','bias'] and toy<=args.bias_toys:
                for label in ['r0','rExpected']:tasks.append(('bias',mass,label,toyorder[toy-1]))
            if args.mode in ['all','gof'] and toy<=args.gof_toys:tasks.append(('gof',mass,None,toy))
    state=R/f'toys_progress_{args.name}.json';active={};genactive={};failedgen=set();finished=0;started=time.time()
    value=dict(status='starting',workers=args.workers,target_bias=sum(t[0]=='bias' for t in tasks),target_gof=sum(t[0]=='gof' for t in tasks),completed={})
    with cf.ThreadPoolExecutor(max_workers=8) as pool,cf.ThreadPoolExecutor(max_workers=1) as genpool:
        while True:
            for future in list(active):
                if future.done():
                    task=active.pop(future)
                    try:result=future.result();print('DONE',result,flush=True)
                    except Exception as error:
                        print('TASKERROR',task,repr(error),flush=True)
                        kind,mass,label,toy=task
                        failure=dict(mass=float(mass),toy=toy,label=label,injection=label,fit_valid=False,valid=False,exception=repr(error),failure_flags=['runner_exception'])
                        if kind=='gof':
                            failure.update(score=None,nominal_valid=False,nominal_diagnostics=dict(available=False,attempted=False,outcome_complete=True,valid=False,flags=['GoF_scheduler_exception'],exception=repr(error)))
                        atomic(result_path(task),failure)
                    finished+=1
            for future in list(genactive):
                if future.done():
                    key=genactive.pop(future)
                    try:
                        result=future.result();print('GEN',key,result,flush=True)
                        if result['exit_code']:failedgen.add(key)
                    except Exception as error:failedgen.add(key);print('GENERROR',key,repr(error),flush=True)
            pending=[task for task in tasks if not result_path(task).exists() and task not in active.values()]
            if not pending and not active:
                break
            ready_masses={m for m in masses if median_ready(m)}
            budget=worker_budget(args)
            for task in pending:
                if len(active)>=budget:break
                if externally_locked(task):continue
                kind,mass,label,toy=task
                if mass not in ready_masses:continue
                if generation_ready(kind,mass,label):active[pool.submit(run_task,task)]=task
                elif not genactive:
                    key=(kind,mass if kind=='bias' else '7.5',label)
                    if key in failedgen:continue
                    if kind=='gof' and not (R/'observed/m7.5/joint/workspace_postfit_b.root').exists():continue
                    genactive[genpool.submit(generate,*key)]=key
            counts={kind:sum(result_path(t).exists() for t in tasks if t[0]==kind) for kind in ['bias','gof']}
            value=dict(status='running',workers=budget,ready_masses=len(ready_masses),
                       target_bias=sum(t[0]=='bias' for t in tasks),target_gof=sum(t[0]=='gof' for t in tasks),
                       completed=counts,active=[list(t) for t in active.values()],
                       generation_active=[list(t) for t in genactive.values()],failed_generation=[list(t) for t in failedgen],
                       elapsed_seconds=time.time()-started,updated_at=time.strftime('%Y-%m-%dT%H:%M:%S%z'))
            atomic(state,value)
            if not active and not genactive and failedgen and all(m in ready_masses for m in masses):
                print('BLOCKED GENERATION',failedgen,flush=True);break
            time.sleep(3)
    value['status']='completed' if all(result_path(t).exists() for t in tasks) else 'incomplete'
    value['completed']={kind:sum(result_path(t).exists() for t in tasks if t[0]==kind) for kind in ['bias','gof']}
    value['active']=[];value['elapsed_seconds']=time.time()-started;atomic(state,value)
    print('POINT AND GOF QUEUE',value,flush=True)
    if value['status']=='completed' and args.mode in ['all','bias'] and not args.no_profiles:
        value['status']='control_profiles_running';atomic(state,value)
        cmd=[PY,str(S/'toys_profiles_run.py'),'--mass']+masses+['--workers',str(worker_budget(args))]
        if args.worker_budget_file:cmd+=['--worker-budget-file',str(args.worker_budget_file)]
        receipt=subprocess_job(cmd,R/'bias/profiles/final_queue',43200)
        value['profile_queue_receipt']=receipt
        value['status']='completed' if receipt['exit_code']==0 else 'control_profiles_incomplete'
        atomic(state,value)
    print('FINAL',value,flush=True)
    return 0 if value['status']=='completed' else 1


if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--workers',type=int,default=4)
    p.add_argument('--mass',nargs='+',default=[f'{k/10:.1f}' for k in range(60,121)])
    p.add_argument('--mode',choices=['all','bias','gof'],default='all')
    p.add_argument('--bias-toys',type=int,default=50)
    p.add_argument('--gof-toys',type=int,default=100)
    p.add_argument('--name',default='main')
    p.add_argument('--worker-budget-file',type=Path)
    p.add_argument('--no-profiles',action='store_true')
    args=p.parse_args()
    assert 1<=args.workers<=8 and 1<=args.bias_toys<=50 and 1<=args.gof_toys<=100
    raise SystemExit(main(args))
