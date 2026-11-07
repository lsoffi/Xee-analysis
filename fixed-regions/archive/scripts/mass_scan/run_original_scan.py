#!/usr/bin/env python3
import argparse,concurrent.futures,datetime,json,os,pathlib,subprocess,time
P=pathlib.Path('/Users/imb/Documents/Codex/Analisi-statistica-Xee');run=P/'runs/20261029_original_mass_scan';ap=argparse.ArgumentParser();ap.add_argument('--workers',type=int,default=8);ap.add_argument('--pilot',action='store_true');args=ap.parse_args()
inv=json.loads((run/'input_inventory.json').read_text());masses=sorted(set(x['mass_label'] for x in inv),key=float);lookup={(x['mass_label'],x['era']):x for x in inv};eras=['2022','2022EE','2023','2023BPix','allYears'];freeze=','.join('CMS_EXO25020_bkgScale_upsilon1s_'+e for e in eras[:-1]);tasks=[]
order=['9.0']+[m for m in masses if m!='9.0']
if args.pilot:order=['9.0']
for m in order:
 for role,era in [('limits','allYears'),('significance','allYears'),('fitdiagnostics','allYears')]+[('fitdiagnostics',e) for e in eras[:-1]]:
  wd=run/'fits'/('mass_'+m)/(role if era=='allYears' else role+'/'+era);wd.mkdir(parents=True,exist_ok=True);inp=lookup[m,era]['path'];method={'limits':'AsymptoticLimits','significance':'Significance','fitdiagnostics':'FitDiagnostics'}[role]
  cmd=[str(P/'bin/combine'),'-M',method,'-d',inp,'-m',m]
  if role=='fitdiagnostics':cmd+=['--saveNormalizations','--saveShapes','--keepFailures','--rMin','-30','--rMax','30','--cminDefaultMinimizerStrategy','0','--freezeParameters',freeze,'--robustFit','1','--X-rtd','MINIMIZER_freezeDisassociatedParams','--cminRunAllDiscreteCombinations','-n','_inclusive_'+era,'-v','3']
  elif role=='limits':cmd+=['--cminDefaultMinimizerStrategy','0','--X-rtd','MINIMIZER_freezeDisassociatedParams','-n','_original_allYears']
  else:cmd+=['--redefineSignalPOIs','r','--rMin','-50','--rMax','50','-n','.observed_significance_allYears','-v','1','--cminDefaultMinimizerStrategy','0','--X-rtd','MINIMIZER_freezeDisassociatedParams','--cminRunAllDiscreteCombinations']
  tasks.append(dict(mass=float(m),mass_label=m,role=role,era=era,input=inp,cwd=str(wd),command=cmd))
(run/'commands.json').write_text(json.dumps(tasks,indent=2))
def work(t):
 wd=pathlib.Path(t['cwd']);receipt=wd/'receipt.json'
 if receipt.exists():
  old=json.loads(receipt.read_text())
  if old.get('exit_code')==0:return old
 start=time.time();env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1');(wd/'command.json').write_text(json.dumps(t,indent=2));code=None;error=None
 try:
  with (wd/'combine.log').open('w') as out:
   proc=subprocess.run(t['command'],cwd=wd,env=env,stdout=out,stderr=subprocess.STDOUT,timeout=1800);code=proc.returncode
 except subprocess.TimeoutExpired:error='timeout after 1800 seconds';code=124
 except Exception as e:error=repr(e);code=-1
 r={**t,'exit_code':code,'elapsed_seconds':time.time()-start,'error':error,'finished_at':datetime.datetime.now(datetime.timezone.utc).isoformat()};receipt.write_text(json.dumps(r,indent=2));return r
results=[];start=time.time();progress_name='progress_pilot.json' if args.pilot else 'progress.json'
with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as pool:
 for future in concurrent.futures.as_completed([pool.submit(work,t) for t in tasks]):
  r=future.result();results.append(r);progress=dict(completed=len(results),total=len(tasks),failed=sum(a['exit_code']!=0 for a in results),elapsed_seconds=time.time()-start,last={k:r[k] for k in ['mass_label','role','era','exit_code']});tmp=run/(progress_name+'.tmp');tmp.write_text(json.dumps(progress,indent=2));tmp.replace(run/progress_name)
  with (run/'events.jsonl').open('a') as f:f.write(json.dumps(progress)+'\n')
  print(json.dumps(progress),flush=True)
(run/('receipts_pilot.json' if args.pilot else 'receipts.json')).write_text(json.dumps(results,indent=2));print('COMPLETED',len(results),'tasks',flush=True)
