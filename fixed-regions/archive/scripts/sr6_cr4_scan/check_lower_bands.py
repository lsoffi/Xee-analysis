from pathlib import Path
import json,subprocess,time,concurrent.futures
import numpy as np,uproot
P=Path.cwd();R=P/'runs/20261030_sr6_cr4_scan'
def job(arg):
 mass,key,quantile=arg;d=R/f'm{mass}/joint';a=json.loads((d.parent/'results.json').read_text());r=a['limits'][key];dest=d/'expected_profile'/f'check_{key}';dest.mkdir(parents=True,exist_ok=True);cmd=json.loads((d/'limits/command.json').read_text())+['--singlePoint',str(r)];cmd[cmd.index('--cminDefaultMinimizerStrategy')+1]='1';cmd[cmd.index('--cminDefaultMinimizerTolerance')+1]='0.001';(dest/'command.json').write_text(json.dumps(cmd));start=time.time()
 with (dest/'combine.log').open('w') as out:res=subprocess.run(cmd,cwd=dest,stdout=out,stderr=subprocess.STDOUT,timeout=600)
 arr=uproot.open(next(dest.glob('higgsCombine*.root')))['limit'].arrays(library='np');cls=float(arr['limit'][np.argmin(abs(arr['quantileExpected']-quantile))]);v=dict(mass=float(mass),key=key,r=r,quantile=quantile,CLs=cls,verified=abs(cls-.05)<.0005,seconds=time.time()-start);(dest/'result.json').write_text(json.dumps(v,indent=2));print(v,flush=True);return v
with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:rows=list(pool.map(job,[(m,k,q) for m in ['6.0','7.5','9.0'] for k,q in [('expected_2p5',.025),('expected_16',.16),('expected_50',.5)]]))
(R/'lower_bands_check.json').write_text(json.dumps(rows,indent=2))
