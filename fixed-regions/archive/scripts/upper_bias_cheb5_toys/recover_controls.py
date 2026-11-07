"""Recover finished v2 control-point estimates from saved fits, never launch Combine."""
from pathlib import Path
import fcntl, hashlib, json, math, subprocess, time
from datetime import datetime
from zoneinfo import ZoneInfo
P=Path(__file__).resolve().parents[2]
R=P/'runs/20261030_upper_bias_cheb5_toys'
S=P/'scripts/upper_bias_cheb5_toys'

def receipts(directory):
    return {str(f.relative_to(directory)):hashlib.sha256(f.read_bytes()).hexdigest()
            for f in directory.glob('**/receipt.json')}

def anchor_valid(a):
    q=a.get('selected',{})
    d=q.get('deviance_with_aux')
    return bool(q.get('readable') and q.get('fit_valid') and q.get('receipt',{}).get('exit_code')==0
                and d is not None and math.isfinite(d))

entries=[]
for mass in [f'{i/10:.1f}' for i in range(105,121)]:
    for label in ['r0','rExpected']:
        d=R/f'fits/m{mass}/{label}/toy001';rp=d/'result.json'
        item=dict(mass=float(mass),injection=label,directory=str(d));entries.append(item)
        if not rp.exists():item['decision']='skip_incomplete';continue
        q=json.loads(rp.read_text())
        if q.get('protocol')!='two_fixed_anchors_v2':item['decision']='skip_legacy_protocol';continue
        if q.get('recovered_control_point') or q.get('profile_failed_point_retained'):
            item['decision']='skip_already_recovered_or_retained';continue
        if q.get('fit_valid'):
            item['decision']='skip_point_already_usable';continue
        mp=d/'multistart_v2.json'
        if not mp.exists():item['decision']='skip_no_anchor_audit';continue
        m=json.loads(mp.read_text())
        if not any(anchor_valid(a) for a in m.get('anchors',[])):
            item['decision']='skip_no_valid_anchor';continue
        lock=(d/'worker.lock').open('a')
        try:fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        except BlockingIOError:
            lock.close();item['decision']='skip_locked_worker';continue
        # Completed-v2 workers are idempotent. Release so the recovery worker owns its normal lock.
        fcntl.flock(lock,fcntl.LOCK_UN);lock.close()
        before=receipts(d);before_hash=hashlib.sha256(rp.read_bytes()).hexdigest()
        command=[str(P/'bin/python'),str(S/'toy_worker.py'),'--mass',mass,'--injection',label,
                 '--toy','1','--recover-control-point']
        item['command']=command;start=time.time()
        with (d/'control_point_recovery.log').open('w') as log:
            code=subprocess.run(command,cwd=P,stdout=log,stderr=subprocess.STDOUT,timeout=120).returncode
        result=json.loads(rp.read_text());after=receipts(d)
        item.update(decision='recovered' if result.get('recovered_control_point') and result.get('fit_valid') else 'recovery_not_successful',
                    exit_code=code,seconds=time.time()-start,
                    before_result_sha256=before_hash,after_result_sha256=hashlib.sha256(rp.read_bytes()).hexdigest(),
                    combine_receipts_unchanged=(before==after),before_receipt_count=len(before),after_receipt_count=len(after),
                    r=result.get('r'),fit_valid=result.get('fit_valid'),interval_valid=result.get('interval_valid'),
                    profile_failed_point_retained=result.get('profile_failed_point_retained'))
        assert before==after,'Recovery must not run or modify any Combine fit receipts'
        print(json.dumps(item),flush=True)
stamp=datetime.now(ZoneInfo('Europe/Rome')).isoformat()
audit=dict(time=stamp,procedure='NO-FIT control toy1 recovery; completed v2 only; existing valid anchor required',
           entries=entries,counts={key:sum(e['decision']==key for e in entries) for key in sorted({e['decision'] for e in entries})})
archive=R/'control_point_recovery_audits';archive.mkdir(exist_ok=True)
file=archive/(stamp.replace(':','-')+'.json');file.write_text(json.dumps(audit,indent=2))
(R/'control_point_recovery_audit.json').write_text(json.dumps(audit,indent=2))
print(json.dumps(audit['counts']),flush=True)
