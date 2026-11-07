from pathlib import Path
import json,subprocess,concurrent.futures
P=Path.cwd();R=P/'runs/20261029_cr_windows/baseline/degree2/gof/toys';pars=json.loads((P/'runs/20261029_upsilon_control/parameters.json').read_text())['control']['b']['parameters'];coeff={p['name']:p['value'] for p in pars if p['name'].startswith('local_c')}
def run(shift):
 d=R/f'start{shift}';d.mkdir(exist_ok=True);c=json.loads((R/'fit2/command.json').read_text());c[c.index('--cminDefaultMinimizerStrategy')+1]='2';c[c.index('--cminDefaultMinimizerTolerance')+1]='0.00001';c[c.index('--setParameters')+1]='r=0,'+','.join(f'{k}={v+shift}' for k,v in coeff.items());(d/'command.json').write_text(json.dumps(c,indent=2))
 with (d/'combine.log').open('w') as f:res=subprocess.run(c,cwd=d,stdout=f,stderr=subprocess.STDOUT)
 (d/'receipt.json').write_text(json.dumps({'exit_code':res.returncode}));assert res.returncode==0
with concurrent.futures.ThreadPoolExecutor(max_workers=2) as ex:list(ex.map(run,[-.05,.05]))
