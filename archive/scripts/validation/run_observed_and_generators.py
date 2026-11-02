import pathlib,json,subprocess,time,concurrent.futures
P=pathlib.Path('/Users/imb/Documents/Codex/Analisi-statistica-Xee');R=P/'runs/20261029_background_validation';C=str(P/'bin/combine')
base=['--rMin','-100','--rMax','100','--cminDefaultMinimizerStrategy','0','--X-rtd','MINIMIZER_freezeDisassociatedParams','--cminRunAllDiscreteCombinations']
def run(t):
 mass,model,role=t;w=R/role/f'm{mass}_{model}';w.mkdir(parents=True,exist_ok=True);input=R/f'workspaces/m{mass}_{model}.root'
 if role=='observed':args=['-M','FitDiagnostics','--saveShapes','--saveNormalizations','--keepFailures','--robustFit','1']
 else:args=['-M','MultiDimFit','--algo','none','--setParameters','r=0','--freezeParameters','r','--saveWorkspace','--saveFitResult']
 cmd=[C,*args,'-d',str(input),'-m',mass,'-n','_validation',*base,'-v','1'];(w/'command.json').write_text(json.dumps(cmd,indent=2));start=time.time()
 with (w/'combine.log').open('w') as f:
  try:result=subprocess.run(cmd,cwd=w,stdout=f,stderr=subprocess.STDOUT,timeout=1800);rc=result.returncode
  except subprocess.TimeoutExpired:rc=124
 receipt=dict(mass=mass,model=model,role=role,exit_code=rc,elapsed=time.time()-start);(w/'receipt.json').write_text(json.dumps(receipt,indent=2));print(json.dumps(receipt),flush=True)
with concurrent.futures.ThreadPoolExecutor(max_workers=6) as pool:list(pool.map(run,[(m,k,'observed') for m in ['7.5','9.5'] for k in ['original','degree5']]+[(m,k,'generators') for m in ['7.5','9.5'] for k in ['cheb7','bernstein']]))
