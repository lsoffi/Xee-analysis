from pathlib import Path
import json,subprocess,concurrent.futures,time
P=Path.cwd();R=P/'runs/20261029_gof_ftest_75';g=R/'degree2/toys_1000';pars=json.loads((P/'runs/20261029_upsilon_control/parameters.json').read_text())['joint']['b']['parameters'];coeff={x['name']:x['value'] for x in pars if x['name'].startswith('local_c')}
def job(shift):
 d=g/f'uniform_start{shift}';d.mkdir(exist_ok=True);c=json.loads((g/'fit_degree2/command.json').read_text());c[c.index('--setParameters')+1]='r=0,'+','.join(f'{k}={v+shift}' for k,v in coeff.items());c[c.index('--cminDefaultMinimizerStrategy')+1]='2';c[c.index('--cminDefaultMinimizerTolerance')+1]='0.00001';(d/'command.json').write_text(json.dumps(c,indent=2));t=time.time()
 with (d/'combine.log').open('w') as f:res=subprocess.run(c,cwd=d,stdout=f,stderr=subprocess.STDOUT)
 (d/'receipt.json').write_text(json.dumps(dict(exit_code=res.returncode,seconds=time.time()-t)));print(shift,res.returncode,round(time.time()-t,1),flush=True);assert res.returncode==0
with concurrent.futures.ThreadPoolExecutor(max_workers=2) as ex:list(ex.map(job,[-.05,.05]))
