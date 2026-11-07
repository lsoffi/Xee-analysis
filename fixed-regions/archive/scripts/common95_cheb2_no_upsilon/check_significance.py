from pathlib import Path
import subprocess,json,concurrent.futures,os
P=Path('/Users/imb/Documents/Codex/Analisi-statistica-Xee');R=P/'runs/20261029_common95_cheb2_no_upsilon'
def run(task):
 mass,stage=task;d=R/f'm{mass}';q=json.loads((d/'quality.json').read_text());cmd=json.loads((d/'significance/command.json').read_text());cmd[cmd.index('--cminDefaultMinimizerStrategy')+1]='2';cmd[cmd.index('--cminDefaultMinimizerTolerance')+1]='0.0001';cmd[cmd.index('-v')+1]='2';cmd+=['--setParameters',','.join(f"{v['name']}={v['value']:.17g}" for v in q[stage]['parameters'])]
 out=d/f'significance_seed_{stage}';out.mkdir(exist_ok=True);(out/'command.json').write_text(json.dumps(cmd,indent=2))
 with (out/'combine.log').open('w') as f:cp=subprocess.run(cmd,cwd=out,stdout=f,stderr=subprocess.STDOUT,timeout=180,env=dict(os.environ,OMP_NUM_THREADS='1'))
 import re
 log=(out/'combine.log').read_text();print(mass,stage,cp.returncode,re.findall(r'Significance: (.*)',log),flush=True)
with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:list(pool.map(run,[(m,s) for m in ['6.7','6.8','6.9'] for s in ['b','s']]))
