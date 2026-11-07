from pathlib import Path
import json,subprocess,concurrent.futures,time
P=Path.cwd();R=P/'runs/20261030_quick_no_upsilon'
def job(pair):
 a,label=pair;d=R/('m'+a['mass'])/label;cmd=json.loads((d/'fit/command.json').read_text());cmd[cmd.index('-d')+1]=str(d.parent/'generator.root');cmd+=['--freezeParameters',','.join(y['name'] for y in a['Y_scales'])];out=d/'reference';out.mkdir(exist_ok=True);(out/'command.json').write_text(json.dumps(cmd,indent=2));t=time.time()
 with (out/'combine.log').open('w') as f:r=subprocess.run(cmd,cwd=out,stdout=f,stderr=subprocess.STDOUT,timeout=300)
 (out/'receipt.json').write_text(json.dumps(dict(exit_code=r.returncode,seconds=time.time()-t)));print(a['mass'],label,flush=True)
with concurrent.futures.ThreadPoolExecutor(max_workers=2) as ex:list(ex.map(job,[(a,l) for a in json.loads((R/'cases.json').read_text()) for l in ['zero','injected']]))
