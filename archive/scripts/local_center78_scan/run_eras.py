import pathlib,json,subprocess,concurrent.futures,time,os
P=pathlib.Path('/Users/imb/Documents/Codex/Analisi-statistica-Xee');R=P/'runs/20261029_local_center78_scan';tasks=[]
for a in json.loads((R/'era_inputs.json').read_text()):
 for role in ['fitdiagnostics','significance']:tasks.append((a,role))
def run(task):
 a,role=task;d=pathlib.Path(a['path']).parent/role;d.mkdir(exist_ok=True);cmd=[str(P/'bin/combine'),'-M','FitDiagnostics' if role=='fitdiagnostics' else 'Significance','-d',a['path'],'-m','7.8','-n','_era','--cminDefaultMinimizerStrategy','1','--cminDefaultMinimizerTolerance','0.001','--X-rtd','MINIMIZER_freezeDisassociatedParams','--rMin','-30','--rMax','30','-v','1']
 if role=='fitdiagnostics':cmd+=['--saveShapes','--saveNormalizations','--keepFailures','--robustFit','1','--stepSize','0.0001','--setRobustFitStrategy','1']
 else:cmd+=['--redefineSignalPOIs','r']
 (d/'command.json').write_text(json.dumps(cmd,indent=2));t=time.time()
 with (d/'combine.log').open('w') as f:p=subprocess.run(cmd,cwd=d,stdout=f,stderr=subprocess.STDOUT,env=dict(os.environ,OMP_NUM_THREADS='1'),timeout=300)
 out=dict(**a,role=role,exit_code=p.returncode,seconds=time.time()-t);(d/'receipt.json').write_text(json.dumps(out,indent=2));print(out,flush=True)
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as e:list(e.map(run,tasks))
