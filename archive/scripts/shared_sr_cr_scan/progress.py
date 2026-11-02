"""Run each incremental render in a fresh process, picking up report updates."""
from pathlib import Path
import argparse
import json
import subprocess
import sys
import time

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--watch', action='store_true')
    parser.add_argument('--interval', type=int, default=60)
    args = parser.parse_args()
    run = Path.cwd()/'runs/20261031_shared_sr_cr_scan'
    while True:
        result = subprocess.run([sys.executable, str(Path(__file__).with_name('report.py'))])
        try:
            complete = json.loads((run/'progress.json').read_text()).get('completed', False)
        except (OSError, json.JSONDecodeError):
            complete = False
        if not args.watch or complete:
            raise SystemExit(result.returncode)
        time.sleep(args.interval)
