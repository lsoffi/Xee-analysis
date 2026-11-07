"""Degree-3 saturated GoF: 100 shared datasets, 16 mass fits per dataset."""
from pathlib import Path
import argparse, concurrent.futures, fcntl, hashlib, json, math, os, subprocess, time

P = Path.cwd()
R = P / 'runs/20261030_upper_envelope_scan'
S = P / 'scripts/upper_envelope_scan'
PY = str(P / 'bin/python')
ENV = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1')
MASSES = ['10.8'] + [f'{k / 10:.1f}' for k in range(105, 121) if k != 108]


def generate():
    directory = R / 'gof/degree3/generation'
    directory.mkdir(parents=True, exist_ok=True)
    with (directory / 'degree3_generation.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        receipt = directory / 'generation.runner_receipt.json'
        source = R / 'degree3/m10.8/joint/workspace_postfit_b.root'
        ready = R / 'degree3/m10.8/results.json'
        announced = False
        while not (source.exists() and ready.exists()):
            if not announced:
                print('WAIT B-only generation workspace:', source, flush=True)
                announced = True
            time.sleep(5)
        source_hash = hashlib.sha256(source.read_bytes()).hexdigest()
        if receipt.exists():
            old = json.loads(receipt.read_text())
            if old['exit_code'] == 0 and (directory / 'truth.json').exists() and list(directory.glob('higgsCombine*.root')):
                if old['source_sha256'] != source_hash:
                    raise RuntimeError('Generation source changed after toys were produced; refusing implicit regeneration')
                return
        command = [PY, str(S / 'generator.py'), 'gof', '3', '10.8', '0']
        (directory / 'generation.runner_command.json').write_text(json.dumps(command, indent=2))
        start = time.time()
        with (directory / 'generation.runner.log').open('w') as log:
            try:
                code = subprocess.run(command, cwd=P, stdout=log, stderr=subprocess.STDOUT,
                                      timeout=600, env=ENV).returncode
            except subprocess.TimeoutExpired:
                code = 124
        receipt.write_text(json.dumps(dict(exit_code=code, seconds=time.time() - start,
                                           source=str(source), source_sha256=source_hash), indent=2))
        if code:
            raise RuntimeError(f'Degree3 GoF generation failed: {code}')
        print('GOF3 GENERATED 100 B-only toy datasets', flush=True)


def job(toy):
    directory = R / f'gof/degree3/toy{toy:03d}'
    directory.mkdir(parents=True, exist_ok=True)
    if all((directory / f'm{mass}/result.json').exists() for mass in MASSES):
        return 0
    command = [PY, str(S / 'toy_worker.py'), 'gof', '3', str(toy)] + MASSES
    (directory / 'worker.runner_command.json').write_text(json.dumps(command, indent=2))
    start = time.time()
    with (directory / 'worker.runner.log').open('w') as log:
        try:
            code = subprocess.run(command, cwd=P, stdout=log, stderr=subprocess.STDOUT,
                                  timeout=7200, env=ENV).returncode
        except subprocess.TimeoutExpired:
            code = 124
    count = sum((directory / f'm{mass}/result.json').exists() for mass in MASSES)
    receipt = dict(exit_code=code, seconds=time.time() - start, scores=count, target_scores=16)
    (directory / 'worker.runner_receipt.json').write_text(json.dumps(receipt, indent=2))
    print('GOF3', toy, receipt, flush=True)
    return code or int(count != 16)


def progress():
    base = R / 'gof/degree3'
    rows = [json.loads(p.read_text()) for p in base.glob('toy*/m*/result.json')]
    value = dict(target_scores=1600, scores=len(rows), target_datasets=100,
                 complete_datasets=sum(all((base / f'toy{toy:03d}/m{m}/result.json').exists() for m in MASSES)
                                       for toy in range(1, 101)),
                 finite_scores=sum(math.isfinite(q['score']) for q in rows),
                 CR_flagged_scores=sum(not q['CR_valid'] for q in rows),
                 r_boundary_scores=sum(q['r_boundary'] for q in rows),
                 interpretation='100 B-only datasets reused across16 masses; r-free saturated fits')
    (base / 'progress.json').write_text(json.dumps(value, indent=2))
    print('GOF3 PROGRESS', value, flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--workers', type=int, default=4)
    parser.add_argument('--toy-start', type=int, default=1)
    parser.add_argument('--toy-stop', type=int, default=100)
    args = parser.parse_args()
    if not 1 <= args.workers <= 8:
        parser.error('Workers must be between1 and8')
    if not 1 <= args.toy_start <= args.toy_stop <= 100:
        parser.error('Toy range must be a subset of1–100')
    generate()
    progress()
    failures = 0
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as pool:
        for future in concurrent.futures.as_completed([pool.submit(job, toy) for toy in range(args.toy_start, args.toy_stop + 1)]):
            try:
                failures += int(future.result() != 0)
            except Exception as error:
                failures += 1
                print('GOF3 EXCEPTION', repr(error), flush=True)
            progress()
    print('GOF3 FINISHED; failed toy jobs:', failures, flush=True)
    raise SystemExit(1 if failures else 0)
