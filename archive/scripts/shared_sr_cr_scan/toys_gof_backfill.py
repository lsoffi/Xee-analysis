"""Attach nominal-state diagnostics to earlier GoF scores without replacing them."""
from pathlib import Path
import argparse, concurrent.futures as cf, fcntl, json, math, shutil, subprocess, time, os

P=Path(__file__).resolve().parents[2]
R=P/'runs/20261031_shared_sr_cr_scan'
S=P/'scripts/shared_sr_cr_scan'
PY=str(P/'bin/python')
ENV=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',PYTHONDONTWRITEBYTECODE='1')


def audit_one(directory):
    from gof_diagnostics import DIAGNOSTIC_COMBINE,RTD_FLAG,read_nominal_diagnostics
    from toys_helpers import atomic_json,execute
    directory=Path(directory)
    with (directory/'worker.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX)
        path=directory/'result.json';q=json.loads(path.read_text())
        if q.get('nominal_diagnostics',{}).get('available'):
            return
        backup=directory/'before_nominal_diagnostics'
        backup.mkdir(exist_ok=True)
        if not (backup/'result.json').exists():shutil.copy2(path,backup/'result.json')
        diagnostic_dir=directory/'nominal_diagnostic';diagnostic_dir.mkdir(exist_ok=True)
        copy=diagnostic_dir/'workspace.root'
        if not copy.exists():shutil.copy2(directory/'workspace.root',copy)
        reused=directory.parent.name=='toy007' and directory.name=='m6.0'
        if reused:
            files=list((P/'work/shared_gof_diagnostic_runtime/v2/validation/toy007_6.0').glob('**/higgsCombine*.root'))
            if len(files)!=1:reused=False
        if not reused:
            command=json.loads((directory/'command.json').read_text())
            command[0]=str(DIAGNOSTIC_COMBINE)
            command=[str(copy) if x==str(directory/'workspace.root') else x for x in command]
            command+=['--X-rtd',RTD_FLAG]
            receipt=execute(command,diagnostic_dir,600)
            files=list(diagnostic_dir.glob('higgsCombine*.root'))
        else:receipt=dict(exit_code=0,seconds=0,reused_validated_run=True)
        if receipt['exit_code']==0 and len(files)==1:
            try:
                a=read_nominal_diagnostics(files[0],copy,expected_score=q.get('score'))
                a['reused_builder_validation']=reused
            except Exception as error:
                a=dict(available=False,attempted=True,outcome_complete=True,valid=False,flags=['diagnostic_parse_exception'],exception=repr(error))
        else:a=dict(available=False,attempted=True,outcome_complete=True,valid=False,flags=['diagnostic_execution_failed'],receipt=receipt)
        atomic_json(diagnostic_dir/'result.json',a)
        q.update(nominal_diagnostics=a,joint_EDM_available=a.get('available',False),
                 nominal_continuum_positivity_checked=a.get('available',False),nominal_status=a.get('status'),
                 nominal_edm=a.get('edm'),nominal_valid=a.get('valid',False),
                 diagnostic_runtime=str(DIAGNOSTIC_COMBINE),diagnostic_backfill=True,
                 raw_score_preserved=True,nominal_diagnostic_file=str(diagnostic_dir/'result.json'))
        q['valid']=bool(q.get('score') is not None and math.isfinite(q['score']) and q.get('CR_valid',False)
                        and not q.get('r_boundary',False) and not q.get('numerical_log_flags',[]) and a.get('valid',False))
        atomic_json(path,q)
        print('DIAGNOSTIC',str(directory.relative_to(R)),a.get('valid'),a.get('flags'),flush=True)


def main(args):
    if args.one:
        try:
            audit_one(args.one)
        except Exception as error:
            # Preparation failures are terminal audit outcomes, never missing scores.
            from toys_helpers import atomic_json
            directory=Path(args.one)
            with (directory/'worker.lock').open('a') as lock:
                fcntl.flock(lock,fcntl.LOCK_EX)
                path=directory/'result.json';q=json.loads(path.read_text())
                backup=directory/'before_nominal_diagnostics';backup.mkdir(exist_ok=True)
                if not (backup/'result.json').exists():shutil.copy2(path,backup/'result.json')
                a=dict(available=False,attempted=True,outcome_complete=True,valid=False,flags=['diagnostic_preparation_exception'],exception=repr(error))
                q.update(nominal_diagnostics=a,nominal_valid=False,valid=False,diagnostic_backfill=True,raw_score_preserved=True)
                atomic_json(path,q)
            print('AUDIT_FAILED',str(directory),repr(error),flush=True)
        return
    start=time.time();allcount=0
    while True:
        pending=[]
        for p in sorted((R/'gof').glob('toy*/m*/result.json')):
            q=json.loads(p.read_text())
            audit=q.get('nominal_diagnostics',{})
            if q.get('score') is not None and not audit.get('available') and not audit.get('outcome_complete') and not q.get('diagnostic_backfill'):
                pending.append(p.parent)
        if not pending:break
        def job(d):
            cmd=[PY,str(Path(__file__).resolve()),'--one',str(d)]
            target=d/'nominal_diagnostic';target.mkdir(exist_ok=True)
            with (target/'runner.log').open('w') as log:
                code=subprocess.run(cmd,cwd=P,env=ENV,stdout=log,stderr=subprocess.STDOUT,timeout=900).returncode
            return str(d),code
        with cf.ThreadPoolExecutor(max_workers=args.workers) as pool:
            for future in cf.as_completed([pool.submit(job,d) for d in pending]):
                directory,code=future.result();allcount+=1;print('BACKFILL',directory,code,flush=True)
        time.sleep(2)
    out=dict(status='completed',audited_records=allcount,seconds=time.time()-start,raw_scores_preserved=True)
    (R/'gof/nominal_backfill.json').write_text(json.dumps(out,indent=2));print(out,flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--one',type=Path);p.add_argument('--workers',type=int,default=1);main(p.parse_args())
