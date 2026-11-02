"""Final queue of the preselected toy1 control profiles; points stay immutable."""
from pathlib import Path
import argparse, concurrent.futures as cf, fcntl, json, time
from toys_run import R, S, PY, subprocess_job, atomic, worker_budget


def main(args):
    base=R/'bias/profiles';base.mkdir(parents=True,exist_ok=True)
    with (base/'queue.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX)
        subprocess_job([PY,str(S/'toys_refresh_edm.py')],base/'refresh_early',600)
        tasks=[(f'{float(m):.1f}',label) for m in args.mass for label in ['r0','rExpected']]
        path=lambda t:base/f'm{t[0]}/{t[1]}/toy001/result.json'
        def job(task):
            mass,label=task;d=path(task).parent
            if path(task).exists():return dict(existing=True,task=task)
            point=R/f'bias/fits/m{mass}/{label}/toy001/result.json'
            try:p=json.loads(point.read_text())
            except Exception:p={}
            if not p.get('truth_source') or p.get('r') is None:
                atomic(path(task),dict(mass=float(mass),injection=label,toy=1,profile_valid=False,profile_outcome='skipped_unreadable_point',outcome_complete=True,exception='No readable point estimate to profile',primary_point_estimate_unchanged=True,coverage_evaluated=False))
                return dict(task=task,unavailable=True)
            if not p.get('fit_valid') or not p.get('CR_valid'):
                atomic(path(task),dict(mass=float(mass),injection=label,toy=1,profile_valid=False,profile_outcome='skipped_invalid_point_or_CR',outcome_complete=True,
                       point_result=str(point),point_r=p.get('r'),point_valid=p.get('fit_valid',False),CR_valid=p.get('CR_valid',False),point_failure_flags=p.get('failure_flags',[]),
                       reason='A converged positive point fit and CR fit are prerequisites; no selection on bias, Hesse covariance quality or interval size',
                       primary_point_estimate_unchanged=True,coverage_evaluated=False))
                return dict(task=task,skipped_invalid_point_or_CR=True)
            cmd=[PY,str(S/'toys_worker.py'),'--mass',mass,'--injection',label,'--toy','1','--finish-profile']
            q=subprocess_job(cmd,d,3000);q['task']=task
            if not path(task).exists():atomic(path(task),dict(mass=float(mass),injection=label,toy=1,profile_valid=False,profile_outcome='worker_timeout_or_failure',outcome_complete=True,exception='Profile worker failed before result',runner=q,primary_point_estimate_unchanged=True,coverage_evaluated=False))
            return q
        active={};start=time.time()
        with cf.ThreadPoolExecutor(max_workers=8) as pool:
            while True:
                for future in list(active):
                    if future.done():
                        task=active.pop(future)
                        try:print('PROFILE',future.result(),flush=True)
                        except Exception as error:atomic(path(task),dict(mass=float(task[0]),injection=task[1],toy=1,profile_valid=False,profile_outcome='scheduler_exception',outcome_complete=True,exception=repr(error),primary_point_estimate_unchanged=True,coverage_evaluated=False))
                pending=[t for t in tasks if not path(t).exists() and t not in active.values()]
                if not pending and not active:break
                budget=worker_budget(args)
                for task in pending:
                    if len(active)>=budget:break
                    active[pool.submit(job,task)]=task
                atomic(base/'progress.json',dict(status='running',completed=sum(path(t).exists() for t in tasks),target=len(tasks),active=[list(t) for t in active.values()],workers=budget,seconds=time.time()-start))
                time.sleep(3)
        records=[json.loads(path(t).read_text()) for t in tasks]
        atomic(base/'progress.json',dict(status='completed',completed=len(records),target=len(tasks),valid=sum(bool(q.get('profile_valid')) for q in records),lower_minimum_flags=sum(bool(q.get('profile_found_lower_minimum')) for q in records),seconds=time.time()-start,coverage_evaluated=False))
        print('ALL CONTROL PROFILES COMPLETED',len(records),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--mass',nargs='+',default=[f'{k/10:.1f}' for k in range(60,121)]);p.add_argument('--workers',type=int,default=4);p.add_argument('--worker-budget-file',type=Path);main(p.parse_args())
