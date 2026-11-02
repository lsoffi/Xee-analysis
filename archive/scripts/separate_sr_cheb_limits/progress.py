"""One-shot by default; refresh only until completion, runner closure or user stop."""
from pathlib import Path
import argparse,json,subprocess,sys,time
P=Path(__file__).resolve().parents[2];R=P/'runs/20261031_separate_sr_cheb_limits'
if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--manifest',type=Path,default=R/'manifest.json');parser.add_argument('--watch',action='store_true');parser.add_argument('--interval',type=int,default=60);args=parser.parse_args()
    if args.watch and not args.manifest.exists():raise SystemExit('Set the authoritative point manifest before starting the watcher.')
    while True:
        if args.watch and (R/'stop_record.json').exists():raise SystemExit(0)
        result=subprocess.run([sys.executable,str(Path(__file__).with_name('report.py')),'--manifest',str(args.manifest)],cwd=P)
        try:state=json.loads((R/'progress.json').read_text())
        except (OSError,json.JSONDecodeError):state={}
        runner=json.loads((R/'remaining_runner.json').read_text()) if (R/'remaining_runner.json').exists() else {}
        if not args.watch or (result.returncode==0 and (state.get('watcher_done',False) or runner.get('status')=='completed')):raise SystemExit(result.returncode)
        for _ in range(max(10,args.interval)):
            if (R/'stop_record.json').exists():raise SystemExit(0)
            time.sleep(1)
