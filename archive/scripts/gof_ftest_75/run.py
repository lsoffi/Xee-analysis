from pathlib import Path
import subprocess,json,time,sys
P=Path.cwd();R=P/'runs/20261029_gof_ftest_75';N=int(sys.argv[1]) if len(sys.argv)>1 else 10
C=str(P/'bin/combine');opts=['-m','7.5','--snapshotName','background_fit','--setParameters','r=0','--freezeParameters','r','--cminDefaultMinimizerStrategy','1','--cminDefaultMinimizerTolerance','0.001','--X-rtd','MINIMIZER_freezeDisassociatedParams','-v','1']
def run(cmd,d):
 d.mkdir(parents=True,exist_ok=True);(d/'command.json').write_text(json.dumps(cmd,indent=2));start=time.time()
 with (d/'combine.log').open('w') as f:res=subprocess.run(cmd,cwd=d,stdout=f,stderr=subprocess.STDOUT)
 (d/'receipt.json').write_text(json.dumps(dict(exit_code=res.returncode,seconds=time.time()-start)));print(str(d.relative_to(R)),res.returncode,round(time.time()-start,2),flush=True);assert res.returncode==0
for n in [2,3]:
 w=str(R/f'degree{n}/workspace.root')
 run([C,'-M','GoodnessOfFit','-d',w,'--algo','saturated','--fixedSignalStrength','0','-n','_observed']+opts,R/f'degree{n}/observed')
 seed=75000+n;g=R/f'degree{n}/toys_{N}';run([C,'-M','GenerateOnly','-d',w,'--saveToys','--toysFrequentist','--bypassFrequentistFit','--expectSignal','0','--genBinnedChannels',','.join(f'{reg}_ch{i}' for reg in ['SR','CR'] for i in range(1,5)),'-t',str(N),'-s',str(seed),'-n','_toys']+opts,g)
 tf=str(g/f'higgsCombine_toys.GenerateOnly.mH7.5.{seed}.root')
 for fitted in ([2,3] if n==2 else [3]):
  run([C,'-M','GoodnessOfFit','-d',str(R/f'degree{fitted}/workspace.root'),'--algo','saturated','--fixedSignalStrength','0','--toysFile',tf,'--toysFrequentist','--bypassFrequentistFit','-t',str(N),'-s',str(seed),'-n','_gof']+opts,g/f'fit_degree{fitted}')
