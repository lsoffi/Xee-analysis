"""Finalize missing-score QC outcomes from an already-running legacy scheduler.

This never fits, modifies a native score, or marks a scored legacy toy as audited.
"""
from pathlib import Path
import argparse, fcntl, json, os, shutil, time

P=Path(__file__).resolve().parents[2]
R=P/'runs/20261031_shared_sr_cr_scan'


def atomic(path,obj):
    tmp=path.with_suffix(path.suffix+'.qc.tmp');tmp.write_text(json.dumps(obj,indent=2));tmp.replace(path)


def finalize():
    changed=[]
    for path in sorted((R/'gof').glob('toy*/m*/result.json')):
        q=json.loads(path.read_text())
        a=q.get('nominal_diagnostics',{})
        if a.get('outcome_complete') or (q.get('score') is not None and not a.get('available')):continue
        with (path.parent/'worker.lock').open('a') as lock:
            try:fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
            except BlockingIOError:continue
            q=json.loads(path.read_text())
            a=q.get('nominal_diagnostics',{})
            if a.get('outcome_complete') or (q.get('score') is not None and not a.get('available')):continue
            backup=path.parent/'before_QC_completion';backup.mkdir(exist_ok=True)
            if not (backup/'result.json').exists():shutil.copy2(path,backup/'result.json')
            if a.get('available'):
                # Early v2 parser records predate the explicit completion field.
                a.update(attempted=True,outcome_complete=True);q['nominal_diagnostics']=a
                q['QC_completion_schema_updated']=True
            else:
                q.update(valid=False,nominal_valid=False,nominal_diagnostics=dict(available=False,attempted=False,outcome_complete=True,valid=False,
                         flags=['native_GoF_score_unavailable'],upstream_exception=q.get('exception'),upstream_failure_flags=q.get('failure_flags',[])),
                         missing_score_QC_finalized=True)
            atomic(path,q);changed.append(str(path.relative_to(R)))
    return changed


def main(args):
    count=0;start=time.time()
    while True:
        changed=finalize();count+=len(changed)
        if changed:print('FINALIZED',changed,flush=True)
        state=R/'toys_progress_main.json'
        finished=False
        if state.exists():
            q=json.loads(state.read_text());finished=q.get('status') in ['completed','incomplete','control_profiles_incomplete']
        atomic(R/'gof/QC_outcome_finalizer.json',dict(status='completed' if finished or not args.watch else 'watching',finalized=count,seconds=time.time()-start,updated_at=time.strftime('%Y-%m-%dT%H:%M:%S%z')))
        if finished or not args.watch:return
        time.sleep(15)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--watch',action='store_true');main(p.parse_args())
