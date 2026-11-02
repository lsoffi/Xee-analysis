from pathlib import Path
import json,subprocess,time,concurrent.futures,sys
P=Path.cwd();SRC=P/'runs/20261030_sr6_cr4_scan';R=P/'runs/20261030_gof_profiled_combined';C=str(P/'bin/combine')
def execute(cmd,d):
 d.mkdir(parents=True,exist_ok=True)
 if (d/'receipt.json').exists() and json.loads((d/'receipt.json').read_text())['exit_code']==0:return
 (d/'command.json').write_text(json.dumps(cmd,indent=2));start=time.time()
 with (d/'combine.log').open('w') as out:p=subprocess.run(cmd,cwd=d,stdout=out,stderr=subprocess.STDOUT,timeout=14400)
 (d/'receipt.json').write_text(json.dumps(dict(exit_code=p.returncode,seconds=time.time()-start),indent=2))
 if p.returncode:raise RuntimeError(str(d))
def task(m,pilot=False):
 src=SRC/f'm{m:.1f}/joint/workspace.root';d=R/f'm{m:.1f}';common=[C,'-M','GoodnessOfFit',str(src),'--algo','saturated','--cminDefaultMinimizerStrategy','0','--X-rtd','MINIMIZER_freezeDisassociatedParams','--cminRunAllDiscreteCombinations','-m',f'{m:.1f}','--trackParameters','r'];seed=str(304800+round(m*10))
 execute(common+['-n','_observed'],d/'observed')
 n=10 if pilot else 100
 execute(common+['-t',str(n),'--toysFrequentist','--saveToys','--genBinnedChannels',','.join(f'{region}_ch{i}' for region in ['SR','CR'] for i in range(1,5)),'-s',seed,'-n','_toys'],d/('pilot10' if pilot else 'toys'))
 print('DONE',m,n,flush=True)
if '--pilot' in sys.argv:task(7.5,True)
else:
 masses=[7.5]+[k/10 for k in range(60,91) if k!=75]
 with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
  futures=[pool.submit(task,m) for m in masses]
  for future in concurrent.futures.as_completed(futures):
   future.result()
   env=__import__('os').environ.copy();env['MPLCONFIGDIR']=str(P/'work/mpl-cache')
   subprocess.run(['/Users/imb/Documents/Codex/2026-10-27/referenced-chatgpt-conversation-this-is-an/work/fit-env/bin/python',str(P/'scripts/gof_profiled_combined/report.py')],cwd=P,env=env,check=True)
 print('ALL DONE 31',flush=True)
