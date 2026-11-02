from pathlib import Path
import subprocess,json,time,concurrent.futures,os,sys
P=Path.cwd();S=P/'scripts/sr6_cr4_scan';R=P/'runs/20261030_sr6_cr4_scan';C=str(P/'bin/combine');PY=str(P/'bin/python');PLOT='/Users/imb/Documents/Codex/2026-10-27/referenced-chatgpt-conversation-this-is-an/work/fit-env/bin/python'
def execute(cmd,d,log,project=False):
 d.mkdir(parents=True,exist_ok=True);start=time.time()
 with (d/log).open('w') as out:res=subprocess.run(cmd,cwd=P if project else d,stdout=out,stderr=subprocess.STDOUT,timeout=7200,env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1',MPLCONFIGDIR=str(P/'work/mpl-cache')))
 receipt=dict(exit_code=res.returncode,seconds=time.time()-start);(d/(log+'.receipt.json')).write_text(json.dumps(receipt));return receipt
def job(m):
 d=R/('m'+m);d.mkdir(exist_ok=True)
 if (d/'results.json').exists():return
 try:
  for script,args in [('prepare_signal.py',[m,'signal']),('prepare_control.py',[m,'control']),('build_joint.py',[m]),('seed.py',[m])]:
   rr=execute([PY,str(S/script)]+args,d,script+'.log',True)
   if rr['exit_code']:raise RuntimeError(script)
  # Expected uncertainty from prior conditional Asimov study sets a numerical scale only.
  reference=P/'runs/20261030_sr_upper_mass_bias/forward.jsonl';ref=[json.loads(x) for x in reference.read_text().splitlines()];a=next((x for x in ref if x['degree']==6 and abs(x['mass']-float(m))<1e-6 and x['rtrue']==0),None)
  scale=a['error'] if a else .2;rmax=max(2,15*scale)
  for role,method in [('fitdiagnostics','FitDiagnostics'),('limits','AsymptoticLimits'),('significance','Significance')]:
   wd=d/'joint'/role;cmd=[C,'-M',method,'-d',str(d/'joint/workspace.root'),'-m',m,'-n','_local','--cminDefaultMinimizerStrategy','1','--cminDefaultMinimizerTolerance','0.001','--X-rtd','MINIMIZER_freezeDisassociatedParams','--cminRunAllDiscreteCombinations','-v','1']
   if role=='fitdiagnostics':cmd+=['--saveShapes','--saveNormalizations','--keepFailures','--robustFit','0','--minos','none','--preFitValue','0.01','--rMin',str(-rmax),'--rMax',str(rmax)]
   if role=='significance':cmd+=['--redefineSignalPOIs','r','--rMin','-1','--rMax',str(rmax)]
   if role=='limits':cmd+=['--rMax',str(rmax),'--rRelAcc','0.001','--rAbsAcc','0.00001'];cmd[cmd.index('--cminDefaultMinimizerStrategy')+1]='2';cmd[cmd.index('--cminDefaultMinimizerTolerance')+1]='0.00001'
   wd.mkdir(parents=True,exist_ok=True);(wd/'command.json').write_text(json.dumps(cmd,indent=2));receipt=execute(cmd,wd,'combine.log');print(m,role,receipt,flush=True)
   if receipt['exit_code']:raise RuntimeError(role)
  for script in ['quality.py','report.py']:
   receipt=execute([PY if script=='quality.py' else PLOT,str(S/script),m],d,script+'.log',True)
   if receipt['exit_code']:raise RuntimeError(script)
  execute([PLOT,str(S/'summary.py')],R,'summary.log',True)
 except Exception as e:(d/'failure.json').write_text(json.dumps(dict(error=str(e))));print('FAIL',m,e,flush=True)
if __name__=='__main__':
 masses=sys.argv[1:] or [f'{k/10:.1f}' for k in range(60,91)]
 with concurrent.futures.ThreadPoolExecutor(max_workers=3) as ex:list(ex.map(job,masses))
 print('SCAN END',flush=True)
