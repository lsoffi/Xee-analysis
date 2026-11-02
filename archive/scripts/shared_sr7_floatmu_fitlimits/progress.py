"""Render once by default; --watch runs periodic fresh-process updates."""
from pathlib import Path
import argparse
import json
import subprocess
import sys
import time
P=Path(__file__).resolve().parents[2]
R=P/'runs/20261031_shared_sr7_floatmu_fitlimits'
if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--watch',action='store_true');parser.add_argument('--interval',type=int,default=60);args=parser.parse_args()
    while True:
        result=subprocess.run([sys.executable,str(Path(__file__).with_name('report.py'))],cwd=P)
        try:complete=json.loads((R/'progress.json').read_text()).get('watcher_done',False)
        except (OSError,json.JSONDecodeError):complete=False
        if not args.watch or (complete and result.returncode==0):raise SystemExit(result.returncode)
        time.sleep(max(10,args.interval))
