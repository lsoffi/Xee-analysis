"""Resumable bounded runner for the requested 1,600 toy experiment (50/cell)."""
from pathlib import Path
import argparse
import concurrent.futures
import datetime
import fcntl
import json
import os
import signal
import subprocess
import time

P = Path(__file__).resolve().parents[2]
R = P / 'runs/20261030_upper_bias_cheb5_toys'
S = P / 'scripts/upper_bias_cheb5_toys'
ENV = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1')


def atomic_json(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + '.tmp')
    tmp.write_text(json.dumps(value, indent=2, allow_nan=False))
    tmp.replace(path)


def run_job(job, timeout):
    mass, injection, toy = job
    directory = R / f'fits/m{mass}/{injection}/toy{toy:03d}'
    directory.mkdir(parents=True, exist_ok=True)
    if (directory / 'result.json').exists():
        existing = json.loads((directory / 'result.json').read_text())
        if existing.get('protocol') == 'two_fixed_anchors_v2':
            return dict(mass=mass, injection=injection, toy=toy, state='already_complete')
    generation = R / f'generation/m{mass}/{injection}'
    # A generated ensemble must be complete before any worker reads its ROOT file.
    while not (generation / 'truth.json').exists():
        time.sleep(2)
    command = [str(P / 'bin/python'), str(S / 'toy_worker.py'),
               '--mass', mass, '--injection', injection, '--toy', str(toy),
               '--interval', 'profile68' if toy == 1 else 'hesse68',
               '--fit-timeout', '600', '--upgrade-multistart']
    atomic_json(directory / 'runner_command.json', command)
    start = time.time()
    with (directory / 'runner.log').open('w') as log:
        process = subprocess.Popen(command, cwd=P, env=ENV, stdout=log,
                                   stderr=subprocess.STDOUT, start_new_session=True)
        try:
            code = process.wait(timeout=timeout)
        except subprocess.TimeoutExpired:
            os.killpg(process.pid, signal.SIGTERM)
            try:
                process.wait(timeout=10)
            except subprocess.TimeoutExpired:
                os.killpg(process.pid, signal.SIGKILL)
                process.wait()
            code = 124
    receipt = dict(mass=mass, injection=injection, toy=toy, exit_code=code,
                   seconds=time.time() - start, result_exists=(directory / 'result.json').exists())
    atomic_json(directory / 'runner_receipt.json', receipt)
    return receipt


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--workers', type=int, default=8)
    parser.add_argument('--timeout', type=int, default=1800)
    parser.add_argument('--toy-start', type=int, default=1)
    parser.add_argument('--toy-stop', type=int, default=50)
    parser.add_argument('--masses', nargs='+')
    args = parser.parse_args()
    if not 1 <= args.workers <= 8:
        parser.error('Use between one and eight workers.')
    if not 1 <= args.toy_start <= args.toy_stop <= 50:
        parser.error('User reduced the experiment to toy IDs 1–50.')
    R.mkdir(parents=True, exist_ok=True)
    with (R / 'runner.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        masses = args.masses or ['10.8'] + [f'{k / 10:.1f}' for k in range(105, 121) if k != 108]
        jobs = [(mass, injection, toy) for toy in range(args.toy_start, args.toy_stop + 1)
                for mass in masses for injection in ['r0', 'rExpected']]
        start = time.time()
        completed = 0
        failed = 0
        state = dict(status='running', started_at=datetime.datetime.now().astimezone().isoformat(),
                     pid=os.getpid(), workers=args.workers, target_jobs=len(jobs),
                     full_experiment_target=1600, toy_start=args.toy_start, toy_stop=args.toy_stop,
                     interval_policy='Profile controls toy1; all other toys Hesse diagnostic only',
                     masses=masses)
        atomic_json(R / 'runner_state.json', state)
        with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as pool:
            futures = {pool.submit(run_job, job, args.timeout): job for job in jobs}
            for future in concurrent.futures.as_completed(futures):
                try:
                    receipt = future.result()
                    failed += int(receipt.get('exit_code', 0) != 0 or not receipt.get('result_exists', True))
                except Exception as error:
                    failed += 1
                    receipt = dict(job=futures[future], exception=repr(error))
                completed += 1
                state.update(completed_jobs=completed, runner_failures=failed,
                             seconds=time.time() - start, last_job=receipt)
                atomic_json(R / 'runner_state.json', state)
                print(json.dumps(state), flush=True)
        state.update(status='completed' if not failed else 'completed_with_runner_failures',
                     completed_at=datetime.datetime.now().astimezone().isoformat())
        atomic_json(R / 'runner_state.json', state)
        return int(bool(failed))


if __name__ == '__main__':
    raise SystemExit(main())
