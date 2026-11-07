import sys
MASS=sys.argv[1]
import pathlib,subprocess,json,concurrent.futures,time,os
P=pathlib.Path('/Users/imb/Documents/Codex/Analisi-statistica-Xee');R=P/'runs/20261029_local_gaussian_scan_7_8'/('m'+MASS);R.mkdir(parents=True,exist_ok=True)
def run(role):
 d=R/role;d.mkdir(exist_ok=True);method={'fitdiagnostics':'FitDiagnostics','significance':'Significance','limits':'AsymptoticLimits'}[role];cmd=[str(P/'bin/combine'),'-M',method,'-d',str(R/'workspace.root'),'-m',MASS,'-n','_local','--cminDefaultMinimizerStrategy','1','--cminDefaultMinimizerTolerance','0.001','--X-rtd','MINIMIZER_freezeDisassociatedParams','-v','1']
 if role=='fitdiagnostics':cmd+=['--saveShapes','--saveNormalizations','--keepFailures','--robustFit','1','--rMin','-5','--rMax','5']
 if role=='significance':cmd+=['--redefineSignalPOIs','r','--rMin','-5','--rMax','5']
 if role=='limits':cmd+=['--rMax','1']
 (d/'command.json').write_text(json.dumps(cmd,indent=2));t=time.time()
 with (d/'combine.log').open('w') as log:
  p=subprocess.run(cmd,cwd=d,stdout=log,stderr=subprocess.STDOUT,env=dict(os.environ,OMP_NUM_THREADS='1'),timeout=600)
 out=dict(role=role,exit_code=p.returncode,elapsed=time.time()-t);(d/'receipt.json').write_text(json.dumps(out,indent=2));print(out,flush=True)
with concurrent.futures.ThreadPoolExecutor(max_workers=3) as e:list(e.map(run,['fitdiagnostics','significance','limits']))
