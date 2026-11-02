from pathlib import Path
import json,subprocess,time,concurrent.futures
P=Path.cwd();R=P/'runs/20261029_cr_windows';C=str(P/'bin/combine');N=500
opts=['-m','7.5','--snapshotName','control_fit','--setParameters','r=0','--freezeParameters','r','--cminDefaultMinimizerStrategy','1','--cminDefaultMinimizerTolerance','0.001','--X-rtd','MINIMIZER_freezeDisassociatedParams','-v','1']
def run(c,d):
 d.mkdir(parents=True,exist_ok=True);receipt=d/'receipt.json'
 if receipt.exists() and json.loads(receipt.read_text())['exit_code']==0:return
 (d/'command.json').write_text(json.dumps(c,indent=2));t=time.time()
 with (d/'combine.log').open('w') as f:res=subprocess.run(c,cwd=d,stdout=f,stderr=subprocess.STDOUT)
 receipt.write_text(json.dumps(dict(exit_code=res.returncode,seconds=time.time()-t)));assert res.returncode==0;print(str(d.relative_to(R)),round(time.time()-t,1),flush=True)
def job(item):
 tag,seedbase=item
 for n in [2,3]:
  root=R/tag/f'degree{n}'/'gof';w=str(root/'workspace.root');run([C,'-M','GoodnessOfFit','-d',w,'--algo','saturated','--fixedSignalStrength','0','-n','_observed']+opts,root/'observed')
  seed=seedbase+n;g=root/'toys';run([C,'-M','GenerateOnly','-d',w,'--saveToys','--toysFrequentist','--bypassFrequentistFit','--expectSignal','0','--genBinnedChannels','ch1,ch2,ch3,ch4','-t',str(N),'-s',str(seed),'-n','_toys']+opts,g)
  tf=str(g/f'higgsCombine_toys.GenerateOnly.mH7.5.{seed}.root')
  for fitted in ([2,3] if n==2 else [3]):
   run([C,'-M','GoodnessOfFit','-d',str(R/tag/f'degree{fitted}/gof/workspace.root'),'--algo','saturated','--fixedSignalStrength','0','--toysFile',tf,'--toysFrequentist','--bypassFrequentistFit','-t',str(N),'-s',str(seed),'-n','_gof']+opts,g/f'fit{fitted}')
with concurrent.futures.ThreadPoolExecutor(max_workers=2) as ex:list(ex.map(job,[(tag,91000+10*i) for i,tag in enumerate(['baseline','cr90_106','cr91_105','cr92_104'])]))
