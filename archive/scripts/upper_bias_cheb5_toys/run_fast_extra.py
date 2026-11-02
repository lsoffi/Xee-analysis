"""Two additional Hesse-only workers, traversing untouched toy IDs in reverse."""
from pathlib import Path
import concurrent.futures
import datetime
import fcntl
import json
import os
import signal
import subprocess
import time

P=Path(__file__).resolve().parents[2]
R=P/'runs/20261030_upper_bias_cheb5_toys'
S=P/'scripts/upper_bias_cheb5_toys'
PROTOCOL='two_fixed_anchors_v2'
ENV=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1')

def atomic(path,value):
    path.parent.mkdir(parents=True,exist_ok=True)
    tmp=path.with_suffix(path.suffix+'.EXTRA.tmp')
    tmp.write_text(json.dumps(value,indent=2,allow_nan=False));tmp.replace(path)

def ready(path):
    try:return json.loads(path.read_text()).get('protocol')==PROTOCOL
    except (FileNotFoundError,json.JSONDecodeError):return False

def count_protocol():
    return sum(ready(p) for p in (R/'fits').glob('m*/*/toy*/result.json'))

def owned_process_groups(root):
    """Find only descendants of this runner's worker, including Combine sessions."""
    try:
        lines=subprocess.check_output(['ps','-axo','pid=,ppid='],text=True).splitlines()
        pairs=[tuple(map(int,line.split())) for line in lines]
    except (OSError,subprocess.SubprocessError,ValueError):
        return {root}
    descendants={root}
    changed=True
    while changed:
        old=len(descendants)
        descendants.update(pid for pid,parent in pairs if parent in descendants)
        changed=len(descendants)!=old
    groups={root}
    for pid in descendants:
        try:
            if os.getpgid(pid)==pid:groups.add(pid)
        except ProcessLookupError:pass
    return groups

def terminate_own(process):
    groups=owned_process_groups(process.pid)
    for group in groups:
        try:os.killpg(group,signal.SIGTERM)
        except ProcessLookupError:pass
    try:process.wait(timeout=10)
    except subprocess.TimeoutExpired:
        for group in groups:
            try:os.killpg(group,signal.SIGKILL)
            except ProcessLookupError:pass
        process.wait()
    return sorted(groups)

def job(mass,injection,toy):
    assert 2<=toy<=50
    directory=R/f'fits/m{mass}/{injection}/toy{toy:03d}'
    directory.mkdir(parents=True,exist_ok=True)
    result=directory/'result.json'
    if ready(result):return dict(mass=mass,injection=injection,toy=toy,state='already_v2')
    command=[str(P/'bin/python'),str(S/'toy_worker.py'),'--mass',mass,'--injection',injection,'--toy',str(toy),'--interval','hesse68','--fit-timeout','600','--upgrade-multistart']
    atomic(directory/'EXTRA_runner_command.json',command)
    start=time.time();groups=[]
    with (directory/'EXTRA_runner.log').open('w') as log:
        process=subprocess.Popen(command,cwd=P,env=ENV,stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
        try:code=process.wait(timeout=1800)
        except subprocess.TimeoutExpired:groups=terminate_own(process);code=124
    receipt=dict(mass=mass,injection=injection,toy=toy,exit_code=code,seconds=time.time()-start,result_v2=ready(result),terminated_own_process_groups=groups)
    atomic(directory/'EXTRA_runner_receipt.json',receipt)
    return receipt

def main():
    with (R/'EXTRA_runner.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        jobs=[(f'{m/10:.1f}',label,toy) for toy in range(50,1,-1) for m in range(120,104,-1) for label in ['rExpected','r0']]
        state=dict(status='running',pid=os.getpid(),workers=2,started_at=datetime.datetime.now().astimezone().isoformat(),target_jobs=len(jobs),full_experiment_target=1600,protocol=PROTOCOL,completed_jobs=0,extra_completed_results=0,runner_failures=0,skipped_already_v2=0,policy='IDs50..2, masses12..10.5; hesse68; never modify profile-control toy1; shared worker mutex; EXTRA receipts only')
        atomic(R/'EXTRA_runner_state.json',state);start=time.time();todo=iter(jobs)
        with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
            active={}
            def enqueue():
                for params in todo:
                    mass,label,toy=params
                    if ready(R/f'fits/m{mass}/{label}/toy{toy:03d}/result.json'):
                        state['completed_jobs']+=1;state['skipped_already_v2']+=1
                        continue
                    active[pool.submit(job,*params)]=params
                    return True
                return False
            enqueue();enqueue()
            while active:
                done,_=concurrent.futures.wait(active,return_when=concurrent.futures.FIRST_COMPLETED)
                for future in done:
                    params=active.pop(future)
                    try:
                        receipt=future.result()
                        if receipt.get('state')=='already_v2':state['skipped_already_v2']+=1
                        else:
                            state['extra_completed_results']+=int(receipt.get('result_v2',False))
                            state['runner_failures']+=int(receipt.get('exit_code',1)!=0 or not receipt.get('result_v2'))
                    except Exception as error:
                        receipt=dict(job=params,exception=repr(error));state['runner_failures']+=1
                    state['completed_jobs']+=1
                    state.update(seconds=time.time()-start,last_job=receipt,total_protocol_v2=count_protocol())
                    atomic(R/'EXTRA_runner_state.json',state)
                    print(json.dumps(state),flush=True)
                    if state['total_protocol_v2']<1600:enqueue()
            state.update(status='completed' if not state['runner_failures'] else 'completed_with_runner_failures',completed_at=datetime.datetime.now().astimezone().isoformat(),seconds=time.time()-start,total_protocol_v2=count_protocol())
            atomic(R/'EXTRA_runner_state.json',state);print(json.dumps(state),flush=True)
        return int(bool(state['runner_failures']))

if __name__=='__main__':raise SystemExit(main())
