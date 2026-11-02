"""Optional extra GoF workers, consuming the already-generated degree3 toy sample.

Explicit launch only. By default processes toy IDs100..51 in reverse while the
main runner proceeds forwards. Generic toy_worker retains the per-toy file lock.
No toy generation and no change to models, data or completed scores occurs here.
"""
import argparse, concurrent.futures, fcntl, json
from run_gof_degree3 import R, job, progress


def extra_job(toy):
    directory = R / f'gof/degree3/toy{toy:03d}'
    directory.mkdir(parents=True, exist_ok=True)
    with (directory / 'worker.lock').open('a') as lock:
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            print('GOF3 EXTRA skips already-running toy', toy, flush=True)
            return 0  # The main runner owns every ID1..100 and will finish it.
        # Release before launching generic worker, which takes its own lock.
        fcntl.flock(lock, fcntl.LOCK_UN)
    return job(toy)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--workers', type=int, default=4)
    parser.add_argument('--toy-start', type=int, default=51)
    parser.add_argument('--toy-stop', type=int, default=100)
    args = parser.parse_args()
    if not 1 <= args.workers <= 4:
        parser.error('Extra workers must be between1 and4')
    if not 1 <= args.toy_start <= args.toy_stop <= 100:
        parser.error('Toy range must be within1–100')
    G = R / 'gof/degree3/generation'
    receipt = json.loads((G / 'generation.runner_receipt.json').read_text())
    truth = json.loads((G / 'truth.json').read_text())
    assert receipt['exit_code'] == 0 and truth['degree'] == 3 and truth['mode'] == 'gof' and truth['rtrue'] == 0
    assert (G / 'generator.root').exists() and len(list(G.glob('higgsCombine*.root'))) == 1
    # Deliberately reuse the saved generation snapshot, even if a later observed
    # extraction reserializes its source ROOT file. Never implicitly regenerate.
    print('GOF3 EXTRA REUSES existing100 toys; source SHA:', receipt['source_sha256'], flush=True)
    failures = 0
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures = [pool.submit(extra_job, toy) for toy in range(args.toy_stop, args.toy_start - 1, -1)]
        for future in concurrent.futures.as_completed(futures):
            try:
                failures += int(future.result() != 0)
            except Exception as error:
                failures += 1
                print('GOF3 EXTRA EXCEPTION', repr(error), flush=True)
            progress()
    print('GOF3 EXTRA FINISHED; failed toy jobs:', failures, flush=True)
    raise SystemExit(1 if failures else 0)
