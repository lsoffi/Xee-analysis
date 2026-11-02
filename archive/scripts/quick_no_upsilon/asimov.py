from pathlib import Path
import subprocess,json,concurrent.futures
P=Path.cwd();R=P/'runs/20261030_quick_no_upsilon'
def run(cmd,d):
 d.mkdir(parents=True,exist_ok=True);(d/'command.json').write_text(json.dumps(cmd,indent=2))
 with (d/'combine.log').open('w') as f:subprocess.run(cmd,cwd=d,stdout=f,stderr=subprocess.STDOUT,timeout=300)
def job(pair):
 a,label=pair;base=R/('m'+a['mass'])/label;d=base/'asimov';cmd=json.loads((base/'generate/command.json').read_text());cmd[cmd.index('-t')+1]='-1';run(cmd,d/'generate');tf=next((d/'generate').glob('higgsCombine*.root'));cmd=json.loads((base/'fit/command.json').read_text());cmd[cmd.index('-t')+1]='-1';cmd[cmd.index('--toysFile')+1]=str(tf);cmd+=['--saveShapes','--saveNormalizations','--saveWorkspace'];run(cmd,d/'fit')
with concurrent.futures.ThreadPoolExecutor(max_workers=2) as ex:list(ex.map(job,[(a,l) for a in json.loads((R/'cases.json').read_text()) for l in ['zero','injected']]))
