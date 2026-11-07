"""Degree-3 Asimov closure scan; execution is explicit, never on import."""
from pathlib import Path
import argparse, concurrent.futures, fcntl, json, math, os, subprocess, time

P = Path.cwd()
R = P / 'runs/20261030_upper_envelope_scan'
S = P / 'scripts/upper_envelope_scan'
PY = str(P / 'bin/python')
ENV = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1')


def read_json_when_ready(path):
    announced = False
    while True:
        if path.exists():
            try:
                return json.loads(path.read_text())
            except json.JSONDecodeError:
                pass
        if not announced:
            print('WAIT observed input:', path, flush=True)
            announced = True
        time.sleep(5)


def execute(cmd, directory, name):
    directory.mkdir(parents=True, exist_ok=True)
    receipt = directory / (name + '.runner_receipt.json')
    if receipt.exists() and json.loads(receipt.read_text())['exit_code'] == 0:
        return 0
    (directory / (name + '.runner_command.json')).write_text(json.dumps(cmd, indent=2))
    start = time.time()
    with (directory / (name + '.runner.log')).open('w') as log:
        try:
            code = subprocess.run(cmd, cwd=P, stdout=log, stderr=subprocess.STDOUT,
                                  timeout=1800, env=ENV).returncode
        except subprocess.TimeoutExpired:
            code = 124
    receipt.write_text(json.dumps(dict(exit_code=code, seconds=time.time() - start), indent=2))
    return code


def job(mass):
    point = R / f'bias/degree3/m{mass}'
    point.mkdir(parents=True, exist_ok=True)
    # A mass-level lock protects generation as well as the generic worker's fit lock.
    with (point / 'degree3_runner.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        obs = read_json_when_ready(R / f'degree3/m{mass}/results.json')
        median = float(obs['limits']['expected_50'])
        if not math.isfinite(median) or median <= 0:
            raise ValueError(f'Invalid expected median at m={mass}: {median}')
        for label, value in [('r0', 0.0), ('rExpected', median), ('r10Expected', 10 * median)]:
            directory = point / label
            directory.mkdir(exist_ok=True)
            final = directory / f'm{mass}/result.json'
            if final.exists():
                continue
            meta = dict(mass=float(mass), fit_degree=3, generator_degree=3,
                        interpretation='same-degree Asimov closure, not alternative-background bias validation',
                        injection=label, rtrue=value, expected_median_at_generation=median,
                        expected_source=str(R / f'degree3/m{mass}/results.json'),
                        original_data_and_model_unchanged=True)
            (directory / 'degree3_interpretation.json').write_text(json.dumps(meta, indent=2))
            generation = [PY, str(S / 'generator.py'), 'bias', '3', mass, str(value), label]
            code = execute(generation, directory, 'generate')
            if code:
                print('BIAS3 GENERATION FAILED', mass, label, code, flush=True)
                return code
            truth_path = directory / 'generation/truth.json'
            truth = json.loads(truth_path.read_text())
            # Only clarify metadata in the new degree-3 output; preserve all numeric truth.
            truth.update(generator='Chebyshev3 fitted to data; same-degree closure for envelope3',
                         interpretation=meta['interpretation'])
            truth_path.write_text(json.dumps(truth, indent=2))
            fit = [PY, str(S / 'toy_worker.py'), 'bias', '3', mass, label]
            code = execute(fit, directory, 'fit')
            if code or not final.exists():
                print('BIAS3 FIT FAILED', mass, label, code, flush=True)
                return code or 1
            result = json.loads(final.read_text())
            print('BIAS3', mass, label, 'valid', result['valid'],
                  'bias/sigma', result['bias_sigma'], flush=True)
    return 0


def progress():
    base = R / 'bias/degree3'
    rows = [json.loads(p.read_text()) for p in base.glob('m*/*/m*/result.json')]
    value = dict(target=48, completed=len(rows), valid=sum(q['valid'] for q in rows),
                 interpretation='Cheb3 generator / degree3 fit: Asimov closure only')
    (base / 'progress.json').write_text(json.dumps(value, indent=2))
    print('BIAS3 PROGRESS', value, flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--workers', type=int, default=2)
    parser.add_argument('--masses', nargs='+', help='Optional subset; default all 16 masses')
    args = parser.parse_args()
    if not 1 <= args.workers <= 2:
        parser.error('Use one or two bias workers to preserve the shared CPU budget')
    masses = args.masses or ['10.8'] + [f'{k / 10:.1f}' for k in range(105, 121) if k != 108]
    (R / 'bias/degree3').mkdir(parents=True, exist_ok=True)
    failures = 0
    progress()
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as pool:
        for future in concurrent.futures.as_completed([pool.submit(job, m) for m in masses]):
            try:
                failures += int(future.result() != 0)
            except Exception as error:
                failures += 1
                print('BIAS3 EXCEPTION', repr(error), flush=True)
            progress()
    print('BIAS3 FINISHED; failed mass jobs:', failures, flush=True)
    raise SystemExit(1 if failures else 0)
