"""One continuation for at most five failed degree3 CR fits, preserving all toys."""
from pathlib import Path
import argparse, fcntl, json, shutil, subprocess, time
import ROOT
from fit_helpers import readfit

ROOT.gSystem.Load('libHiggsAnalysisCombinedLimit')
P = Path.cwd()
R = P / 'runs/20261030_upper_envelope_scan'
BASE = R / 'gof/degree3'
MASSES = [f'{k / 10:.1f}' for k in range(105, 121)]


def sweep():
    failures = [p for p in sorted(BASE.glob('toy*/cr_result.json')) if not json.loads(p.read_text())['valid']]
    attempts = list(BASE.glob('toy*/cr_start2/receipt.json'))
    if len(failures) > 5:
        print('REFINEMENT CAP: more than5 failed CR datasets; preserve remaining nominal scores', len(failures), flush=True)
        return True
    for p in failures:
        D = p.parent
        if (D / 'cr_start2/receipt.json').exists():
            continue
        if len(attempts) >= 5:
            print('REFINEMENT CAP: five continuations already attempted', flush=True)
            return True
        if not all((D / f'm{m}/result.json').exists() for m in MASSES):
            continue
        lock = (D / 'worker.lock').open('a')
        fcntl.flock(lock, fcntl.LOCK_EX)
        old = json.loads(p.read_text())
        if old['valid'] or not all((D / f'm{m}/result.json').exists() for m in MASSES):
            lock.close()
            continue
        wd = D / 'cr_start2'
        wd.mkdir(exist_ok=True)
        f = ROOT.TFile.Open(str(D / 'cr_start0/workspace.root'))
        w = f.Get('w')
        for q in old['selected']['parameters']:
            if w.var(q['name']):
                w.var(q['name']).setVal(q['value'])
        w.writeToFile(str(wd / 'workspace.root'))
        f.Close()
        cmd = json.loads((D / 'cr_start0/command.json').read_text())
        cmd[cmd.index('-M') + 2] = str(wd / 'workspace.root')
        cmd[cmd.index('--cminDefaultMinimizerTolerance') + 1] = '1e-5'
        (wd / 'command.json').write_text(json.dumps(cmd, indent=2))
        start = time.time()
        with (wd / 'combine.log').open('w') as h:
            try:
                rc = subprocess.run(cmd, cwd=wd, stdout=h, stderr=subprocess.STDOUT, timeout=300).returncode
            except subprocess.TimeoutExpired:
                rc = 124
        if rc == 0 and (wd / 'fitDiagnostics_upper.root').exists():
            q = readfit(wd)
            q['variant'] = 2
            delta = q['deviance'] - old['selected']['deviance']
        else:
            q = dict(status=None, covQual=None, deviance=None)
            delta = None
        receipt = dict(exit_code=rc, seconds=time.time() - start, status=q['status'],
                       covQual=q['covQual'],deviance=q['deviance'],delta_deviance=delta)
        (wd / 'receipt.json').write_text(json.dumps(receipt, indent=2))
        attempts.append(wd / 'receipt.json')
        print('CR3 CONTINUATION', D.name, receipt, flush=True)
        if rc == 0 and q['status'] == 0:
            archive = D / 'initial_failed_CR'
            archive.mkdir(exist_ok=True)
            shutil.copy2(p, archive / 'cr_result.json')
            for m in MASSES:
                (D / f'm{m}').rename(archive / f'm{m}')
            new = dict(old, selected=q,
                       means=[v for v in q['parameters'] if v['name'].startswith('Ymean_')],
                       attempts=old['attempts'] + [q],valid=True,
                       retry_reason='Both initial CR fits had nonzero status; one same-data continuation at tolerance1e-5')
            p.write_text(json.dumps(new, indent=2))
            lock.close()
            toy = str(int(D.name.removeprefix('toy')))
            command = [str(P / 'bin/python'),str(P / 'scripts/upper_envelope_scan/toy_worker.py'),'gof','3',toy]
            (D / 'refined_worker_command.json').write_text(json.dumps(command, indent=2))
            with (D / 'refined_worker.log').open('w') as h:
                rc = subprocess.run(command,cwd=P,stdout=h,stderr=subprocess.STDOUT,timeout=1200).returncode
            (D / 'refined_worker_receipt.json').write_text(json.dumps(dict(exit_code=rc),indent=2))
            print('CR3 RECOMPUTED16 SCORES',D.name,rc,flush=True)
        else:
            lock.close()
    ready = len(list(BASE.glob('toy*/m*/result.json'))) == 1600
    pending = [p for p in BASE.glob('toy*/cr_result.json')
               if not json.loads(p.read_text())['valid'] and not (p.parent / 'cr_start2/receipt.json').exists()]
    return ready and not pending


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--watch',action='store_true')
    args = parser.parse_args()
    while True:
        done = sweep()
        if done or not args.watch:
            print('CR3 REFINEMENT FINISHED',flush=True)
            break
        time.sleep(15)
