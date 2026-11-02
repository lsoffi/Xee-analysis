from pathlib import Path
import json,subprocess,time,concurrent.futures,sys
P=Path.cwd();R=P/'runs/20261030_gof_thesis48';C=str(P/'bin/combine');tasks=json.loads((R/'tasks.json').read_text())
if '--pilot' in sys.argv:tasks=[t for t in tasks if t['era']=='2022' and t['mass']==7.5]
def execute(cmd,d):
 d.mkdir(exist_ok=True)
 if (d/'receipt.json').exists() and json.loads((d/'receipt.json').read_text())['exit_code']==0:return
 (d/'command.json').write_text(json.dumps(cmd,indent=2));start=time.time()
 with (d/'combine.log').open('w') as out:p=subprocess.run(cmd,cwd=d,stdout=out,stderr=subprocess.STDOUT,timeout=3600)
 (d/'receipt.json').write_text(json.dumps(dict(exit_code=p.returncode,seconds=time.time()-start),indent=2))
 if p.returncode:raise RuntimeError(str(d))
def task(t):
 d=Path(t['path']);m=str(t['mass']);seed=str(t['seed']);channels=json.loads((d/'audit.json').read_text())['channels'];opts=['-d',str(d/'workspace.root'),'-m',m,'--snapshotName','background_fit','--cminDefaultMinimizerStrategy','1','--cminDefaultMinimizerTolerance','0.001','--X-rtd','MINIMIZER_freezeDisassociatedParams','--cminRunAllDiscreteCombinations']
 execute([C,'-M','GoodnessOfFit','--algo','saturated','--fixedSignalStrength','0','-n','_observed']+opts,d/'observed')
 execute([C,'-M','GenerateOnly','--saveToys','--toysFrequentist','--expectSignal','0','--genBinnedChannels',','.join(channels),'-t','100','-s',seed,'-n','_toys']+opts,d/'generate')
 tf=next((d/'generate').glob('higgsCombine*.root'))
 execute([C,'-M','GoodnessOfFit','--algo','saturated','--fixedSignalStrength','0','--toysFile',str(tf),'--toysFrequentist','--bypassFrequentistFit','-t','100','-s',seed,'-n','_toys']+opts,d/'toys')
 print('DONE',t['era'],m,flush=True)
tasks.sort(key=lambda t:(t['era']!='2022',t['mass'],t['era']))
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
 for result in pool.map(task,tasks):pass
print('ALL DONE',len(tasks),flush=True)
