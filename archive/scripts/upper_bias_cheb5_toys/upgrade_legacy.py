"""Uniformly apply the two fixed starting anchors to initial production toys."""
from pathlib import Path
import concurrent.futures
import json
import os
import subprocess
import time

P = Path(__file__).resolve().parents[2]
R = P / 'runs/20261030_upper_bias_cheb5_toys'
S = P / 'scripts/upper_bias_cheb5_toys'


def upgrade(path):
    row = json.loads(path.read_text())
    directory = path.parent
    command = [str(P / 'bin/python'), str(S / 'toy_worker.py'), '--mass', str(row['mass']),
               '--injection', row['injection'], '--toy', str(row['toy']), '--upgrade-multistart']
    (directory / 'upgrade_runner_command.json').write_text(json.dumps(command, indent=2))
    start = time.time()
    with (directory / 'upgrade_runner.log').open('w') as log:
        code = subprocess.run(command, cwd=P, stdout=log, stderr=subprocess.STDOUT,
                              env=dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1')).returncode
    receipt = dict(exit_code=code, seconds=time.time() - start)
    (directory / 'upgrade_runner_receipt.json').write_text(json.dumps(receipt, indent=2))
    print(str(directory.relative_to(R)), receipt, flush=True)


if __name__ == '__main__':
    attempted = set()
    while True:
        paths = list(R.glob('fits/*/*/toy*/result.json'))
        legacy = [p for p in paths if json.loads(p.read_text()).get('protocol') != 'two_fixed_anchors_v2']
        todo = [p for p in legacy if str(p) not in attempted]
        if todo:
            with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
                list(pool.map(upgrade, todo))
            attempted.update(str(p) for p in todo)
        state = json.loads((R / 'runner_state.json').read_text())
        if state['status'].startswith('completed') and not todo:
            break
        time.sleep(20)
