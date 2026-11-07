import argparse,pathlib,json,subprocess,time,concurrent.futures,os,datetime,signal
P=pathlib.Path('/Users/imb/Documents/Codex/Analisi-statistica-Xee');R=P/'runs/20261029_alternative_background_scan';ap=argparse.ArgumentParser();ap.add_argument('--stage',choices=['seed','scan'],required=True);ap.add_argument('--pilot',action='store_true');ap.add_argument('--defer-active-pilot',action='store_true');ap.add_argument('--workers',type=int,default=8);args=ap.parse_args();tasks=[]
for a in json.loads((R/('audit_seed.json' if args.stage=='seed' else 'audit_scan.json')).read_text()):
 if args.stage=='scan' and args.pilot and a['mass_label'] not in ['6.0','7.5','9.5','11.0']:continue
 for role in (['seed'] if args.stage=='seed' else ['limits','significance','fitdiagnostics']):
  wd=R/(f"seeds/{a['family']}/{a['selection']}" if role=='seed' else f"fits/{a['family']}/m{a['mass_label']}/{role}");wd.mkdir(parents=True,exist_ok=True)
  method={'seed':'FitDiagnostics','fitdiagnostics':'FitDiagnostics','limits':'AsymptoticLimits','significance':'Significance'}[role];cmd=[str(P/'bin/combine'),'-M',method,'-d',a['path'],'-m',a['mass_label'],'-n','_seed' if role=='seed' else '_alternative','--cminDefaultMinimizerStrategy','0','--X-rtd','MINIMIZER_freezeDisassociatedParams','-v','1']
  if role in ['seed','fitdiagnostics']:cmd+=['--saveShapes','--saveNormalizations','--keepFailures','--rMin','-100','--rMax','100'];cmd+=['--skipSBFit'] if role=='seed' else ['--robustFit','1']
  elif role=='significance':cmd+=['--redefineSignalPOIs','r','--rMin','-100','--rMax','100']
  tasks.append(dict(mass=a['mass'],mass_label=a['mass_label'],family=a['family'],selection=a['selection'],role=role,input=a['path'],cwd=str(wd),command=cmd))
if args.defer_active_pilot:
 pending={t['cwd'] for t in json.loads((R/'commands_scan_pilot.json').read_text()) if not (pathlib.Path(t['cwd'])/'receipt.json').exists()};(R/'deferred_pilot_jobs.json').write_text(json.dumps(sorted(pending),indent=2));tasks=[t for t in tasks if t['cwd'] not in pending]
(R/('commands_'+args.stage+('_pilot' if args.pilot else '')+'.json')).write_text(json.dumps(tasks,indent=2))
def run(t):
 wd=pathlib.Path(t['cwd']);rp=wd/'receipt.json'
 if rp.exists():
  prior=json.loads(rp.read_text())
  if prior['exit_code']==0:return prior
 (wd/'command.json').write_text(json.dumps(t,indent=2));start=time.time();env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1');error=None
 with (wd/'combine.log').open('w') as f:
  proc=subprocess.Popen(t['command'],cwd=wd,env=env,stdout=f,stderr=subprocess.STDOUT,start_new_session=True)
  try:rc=proc.wait(timeout=900)
  except subprocess.TimeoutExpired:
   os.killpg(proc.pid,signal.SIGTERM)
   try:proc.wait(timeout=5)
   except subprocess.TimeoutExpired:os.killpg(proc.pid,signal.SIGKILL);proc.wait()
   rc=124;error='timeout after 900 seconds'
 receipt=dict(**t,exit_code=rc,elapsed_seconds=time.time()-start,error=error,finished_at=datetime.datetime.now(datetime.timezone.utc).isoformat());rp.write_text(json.dumps(receipt,indent=2));return receipt
results=[];start=time.time();name='progress_'+args.stage+('_pilot' if args.pilot else '')+'.json'
with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as pool:
 for future in concurrent.futures.as_completed([pool.submit(run,t) for t in tasks]):
  r=future.result();results.append(r);out=dict(completed=len(results),total=len(tasks),execution_failures=sum(r['exit_code']!=0 for r in results),elapsed_seconds=time.time()-start,last={k:r[k] for k in ['mass_label','family','role','exit_code']});tmp=R/(name+'.tmp');tmp.write_text(json.dumps(out,indent=2));tmp.replace(R/name);print(json.dumps(out),flush=True)
print('FINISHED',len(results),flush=True)
