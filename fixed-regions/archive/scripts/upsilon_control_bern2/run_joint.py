import sys,pathlib,subprocess,json,time,os,signal
MASS=sys.argv[1];P=pathlib.Path('/Users/imb/Documents/Codex/Analisi-statistica-Xee');R=P/'runs/20261029_upsilon_control_bern2'/'joint'
for role,method in [('fitdiagnostics','FitDiagnostics'),('limits','AsymptoticLimits'),('significance','Significance')]:
 d=R/role;d.mkdir(exist_ok=True);receipt=d/'receipt.json'
 if receipt.exists() and json.loads(receipt.read_text()).get('exit_code')==0:continue
 cmd=[str(P/'bin/combine'),'-M',method,'-d',str(R/'workspace.root'),'-m',MASS,'-n','_local','--cminDefaultMinimizerStrategy','1','--cminDefaultMinimizerTolerance','0.001','--X-rtd','MINIMIZER_freezeDisassociatedParams','-v','1']
 if role=='fitdiagnostics':cmd+=['--saveShapes','--saveNormalizations','--keepFailures','--robustFit','0','--minos','none','--rMin','-100','--rMax','100']
 if role=='significance':cmd+=['--redefineSignalPOIs','r','--rMin','-100','--rMax','100']
 if role=='limits':cmd+=['--rMax','100','--rRelAcc','0.001','--rAbsAcc','0.00001']
 (d/'command.json').write_text(json.dumps(cmd,indent=2));t=time.time()
 with (d/'combine.log').open('w') as f:
  proc=subprocess.Popen(cmd,cwd=d,stdout=f,stderr=subprocess.STDOUT,start_new_session=True,env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1'))
  try:rc=proc.wait(timeout=900)
  except subprocess.TimeoutExpired:os.killpg(proc.pid,signal.SIGKILL);proc.wait();rc=124
 result=dict(role=role,exit_code=rc,elapsed=time.time()-t);receipt.write_text(json.dumps(result,indent=2));print(result,flush=True)
