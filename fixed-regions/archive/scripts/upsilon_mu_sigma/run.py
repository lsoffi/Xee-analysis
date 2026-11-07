from pathlib import Path
import subprocess,json,time,concurrent.futures
P=Path.cwd();R=P/'runs/20261029_upsilon_mu_sigma';configs=json.loads((R/'configs.json').read_text())
def job(a):
 d=R/a['tag'];c=[str(P/'bin/combine'),'-M','MultiDimFit','-d',str(d/'workspace.root'),'--algo','singles','--redefineSignalPOIs','upsilon_mu,upsilon_sigma','--floatOtherPOIs','1','--noMCbonly','1','--saveFitResult','--saveWorkspace','--cl','0.68268949','--cminDefaultMinimizerStrategy','2','--cminDefaultMinimizerTolerance','0.001','-n','_shape','-v','1'];(d/'command.json').write_text(json.dumps(c,indent=2));start=time.time()
 with (d/'combine.log').open('w') as f:res=subprocess.run(c,cwd=d,stdout=f,stderr=subprocess.STDOUT)
 (d/'receipt.json').write_text(json.dumps(dict(exit_code=res.returncode,seconds=time.time()-start)));print(a['tag'],res.returncode,round(time.time()-start,2),flush=True)
with concurrent.futures.ThreadPoolExecutor(max_workers=2) as ex:list(ex.map(job,configs))
