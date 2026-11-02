"""Numerical continuation only when native significance misses a known valid minimum."""
from pathlib import Path
import argparse,json,time
from run_observed import P,S,R,C,execute,common,extract,historical_expected,atomic

def run(mass):
    D=R/'observed'/('m'+mass);J=D/'joint';a=json.loads((D/'results.json').read_text())
    if 'Z' not in a or not a.get('fit',{}).get('s',{}).get('point_valid'):return
    if a.get('significance_consistent',True):return
    d=J/'significance_refined'
    if (d/'combine.receipt.json').exists():return
    bound=max(1.,10*abs(a['r']),5*historical_expected(mass))
    cmd=[C,'-M','Significance']+common(J/'workspace_postfit_s.root',mass,tolerance=.00001)+[
         '--redefineSignalPOIs','r','--rMin',str(-bound),'--rMax',str(bound),'--setParameters',f'r={a["r"]:.17g}']
    execute(cmd,d,timeout=600)
    b=extract(mass)
    atomic(d/'audit.json',dict(reason='Native result disagreed with a converged, positive, non-boundary FitDiagnostics minimum',
       initial_Z=a['Z'],FD_Z=a['Z_FD'],refined_Z=b.get('significance_refined_Z'),
       selected_Z=b['Z'],consistent=b['significance_consistent'],r_range=[-bound,bound],physical_model_changed=False))

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--watch',action='store_true');args=parser.parse_args()
    while True:
        files=sorted((R/'observed').glob('m*/results.json'))
        for p in files:run(p.parent.name[1:])
        if not args.watch:break
        finished=sum('Z' in json.loads(p.read_text()) for p in files)
        if finished==61:break
        time.sleep(15)
