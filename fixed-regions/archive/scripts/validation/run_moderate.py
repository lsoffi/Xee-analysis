import pathlib,json,subprocess,time,concurrent.futures,threading
P=pathlib.Path('/Users/imb/Documents/Codex/Analisi-statistica-Xee');R=P/'runs/20261029_background_validation';C=str(P/'bin/combine');lock=threading.Lock();start=time.time();done=0

def execute(cmd,wd,timeout=1800):
 wd.mkdir(parents=True,exist_ok=True)
 if (wd/'receipt.json').exists() and json.loads((wd/'receipt.json').read_text())['exit_code']==0:return
 (wd/'command.json').write_text(json.dumps(cmd,indent=2));t=time.time()
 with (wd/'combine.log').open('w') as f:
  try:rc=subprocess.run(cmd,cwd=wd,stdout=f,stderr=subprocess.STDOUT,timeout=timeout).returncode
  except subprocess.TimeoutExpired:rc=124
 (wd/'receipt.json').write_text(json.dumps(dict(exit_code=rc,elapsed_seconds=time.time()-t),indent=2))
 if rc:raise RuntimeError(f'Failed {wd}: {rc}')

def run(task):
 global done
 mass,gen,inj,batch=task;mu=3.0;seed=401000+(['cheb7','bernstein'].index(gen)*20)+batch
 wd=R/f'toys/m{mass}/{gen}/{inj}/batch_{batch:02d}';g=wd/'generate';(wd).mkdir(parents=True,exist_ok=True)
 (wd/'scenario.json').write_text(json.dumps(dict(mass=float(mass),generator=gen,injection=mu,seed=seed,n_toys=10),indent=2))
 execute([C,'-M','GenerateOnly','-d',str(R/f'workspaces/generator_m{mass}_{gen}.root'),'-m',mass,'--snapshotName','validation_generator','--bypassFrequentistFit','--toysFrequentist','--expectSignal',str(mu),'--saveToys','-t','10','-s',str(seed),'-n','_toys','--genBinnedChannels','ch1,ch2,ch3,ch4'],g)
 tf=g/f'higgsCombine_toys.GenerateOnly.mH{mass}.{seed}.root'
 for fit in ['degree5','original']:
  inp=R/f'workspaces/m{mass}_{fit}'+pathlib.Path('') if False else R/f'workspaces/m{mass}_{fit}{".clean" if fit=="degree5" else ""}.root'
  execute([C,'-M','FitDiagnostics','-d',str(inp),'-m',mass,'--toysFile',str(tf),'--toysFrequentist','-t','10','-s',str(seed),'--rMin','-100','--rMax','100','--keepFailures','--cminDefaultMinimizerStrategy','0','--X-rtd','MINIMIZER_freezeDisassociatedParams','--cminRunAllDiscreteCombinations','-n','_toys','-v','1'],wd/fit)
 with lock:
  done+=1;status=dict(completed_batches=done,total_batches=20,expected_toy_fits=400,elapsed_seconds=time.time()-start);(R/'moderate_progress.json').write_text(json.dumps(status,indent=2));print(json.dumps(status),flush=True)

tasks=[(m,g,i,b) for b in range(10) for m in ['9.5'] for g in ['cheb7','bernstein'] for i in ['moderate']]
errors=[]
with concurrent.futures.ThreadPoolExecutor(max_workers=8) as pool:
 futures={pool.submit(run,t):t for t in tasks}
 for f in concurrent.futures.as_completed(futures):
  try:f.result()
  except Exception as e:errors.append(dict(task=futures[f],error=str(e)));print(str(e),flush=True)
(R/'moderate_execution_errors.json').write_text(json.dumps(errors,indent=2))
