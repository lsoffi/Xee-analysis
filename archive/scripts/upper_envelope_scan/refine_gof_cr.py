from pathlib import Path
import ROOT,json,subprocess,fcntl,shutil
from fit_helpers import readfit
ROOT.gSystem.Load('libHiggsAnalysisCombinedLimit')
P=Path.cwd();R=P/'runs/20261030_upper_envelope_scan'
for p in sorted(R.glob('gof/degree*/toy*/cr_result.json')):
    old=json.loads(p.read_text());D=p.parent
    if old['valid'] or (D/'cr_start2/receipt.json').exists():continue
    # Serialize with every normal worker, and retry solely on failed CR convergence.
    lock=(D/'worker.lock').open('a');fcntl.flock(lock,fcntl.LOCK_EX)
    old=json.loads(p.read_text())
    if old['valid']:lock.close();continue
    wd=D/'cr_start2';wd.mkdir(exist_ok=True)
    f=ROOT.TFile.Open(str(D/'cr_start0/workspace.root'));w=f.Get('w')
    for q in old['selected']['parameters']:
        if w.var(q['name']):w.var(q['name']).setVal(q['value'])
    w.writeToFile(str(wd/'workspace.root'));f.Close()
    cmd=json.loads((D/'cr_start0/command.json').read_text());cmd[cmd.index('-M')+2]=str(wd/'workspace.root');cmd[cmd.index('--cminDefaultMinimizerTolerance')+1]='1e-5';(wd/'command.json').write_text(json.dumps(cmd,indent=2))
    with (wd/'combine.log').open('w') as h:rc=subprocess.run(cmd,cwd=wd,stdout=h,stderr=subprocess.STDOUT,timeout=300).returncode
    q=readfit(wd);q['variant']=2;(wd/'receipt.json').write_text(json.dumps(dict(exit_code=rc,status=q['status'],covQual=q['covQual'],deviance=q['deviance'])))
    print(D,q['status'],q['covQual'],q['deviance'],flush=True)
    if rc==0 and q['status']==0:
        arc=D/'initial_failed_CR';arc.mkdir(exist_ok=True);shutil.copy2(p,arc/'cr_result.json')
        for x in D.glob('m*'):
            if x.is_dir():x.rename(arc/x.name)
        new=dict(old,selected=q,means=[v for v in q['parameters'] if v['name'].startswith('Ymean_')],attempts=old['attempts']+[q],valid=True,retry_reason='Both initial CR fits had nonzero status; same-data continuation at tighter tolerance')
        p.write_text(json.dumps(new,indent=2));lock.close()
        degree=D.parent.name.removeprefix('degree');toy=str(int(D.name.removeprefix('toy')))
        with (D/'refined_worker.log').open('w') as h:rc=subprocess.run([str(P/'bin/python'),str(P/'scripts/upper_envelope_scan/toy_worker.py'),'gof',degree,toy],cwd=P,stdout=h,stderr=subprocess.STDOUT,timeout=1200).returncode
        (D/'refined_worker_receipt.json').write_text(json.dumps(dict(exit_code=rc)))
    else:lock.close()
