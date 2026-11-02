from pathlib import Path
import subprocess,json,time,concurrent.futures,os
P=Path.cwd();S=P/'scripts/sr7_cr3_scan';R=P/'runs/20261030_sr7_cr3_scan';C=str(P/'bin/combine');PY=str(P/'bin/python');PLOT='/Users/imb/Documents/Codex/2026-10-27/referenced-chatgpt-conversation-this-is-an/work/fit-env/bin/python'
def execute(cmd,d,log):
 d.mkdir(parents=True,exist_ok=True);start=time.time()
 with (d/log).open('w') as out:res=subprocess.run(cmd,cwd=P if log.endswith('prepare.log') or log in ['build.log','quality.log','report.log'] else d,stdout=out,stderr=subprocess.STDOUT,timeout=1200,env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MPLCONFIGDIR=str(P/'work/mpl-cache')))
 return dict(exit_code=res.returncode,seconds=time.time()-start)
def job(m):
 d=R/('m'+m);d.mkdir(exist_ok=True)
 for region in ['signal','control']:
  rr=execute([PY,str(S/f'prepare_{region}.py'),m,region],d,region+'_prepare.log')
  if rr['exit_code']:return print('PREP FAIL',m,region,flush=True)
 rr=execute([PY,str(S/'build_joint.py'),m],d,'build.log')
 if rr['exit_code']:return print('BUILD FAIL',m,flush=True)
 rr=execute([PY,str(S/'seed.py'),m],d,'seed_prepare.log')
 if rr['exit_code']:return print('SEED FAIL',m,flush=True)
 ref=json.loads((P/'runs/20261029_cr91_dcb_scan'/('m'+m)/'results.json').read_text());rmax=max(2,8*ref['limits']['expected_50'])
 for role,method in [('fitdiagnostics','FitDiagnostics'),('significance','Significance'),('limits','AsymptoticLimits')]:
  wd=d/'joint'/role;wd.mkdir(exist_ok=True);cmd=[C,'-M',method,'-d',str(d/'joint/workspace.root'),'-m',m,'-n','_local','--cminDefaultMinimizerStrategy','1' if role!='limits' else '2','--cminDefaultMinimizerTolerance','0.001' if role!='limits' else '0.00001','--X-rtd','MINIMIZER_freezeDisassociatedParams','-v','1']
  if role=='fitdiagnostics':cmd+=['--saveShapes','--saveNormalizations','--keepFailures','--robustFit','0','--minos','none','--preFitValue','0.01','--rMin','-100','--rMax','100']
  if role=='significance':cmd+=['--redefineSignalPOIs','r','--rMin','-100','--rMax','100']
  if role=='limits':cmd+=['--rMax',str(rmax),'--rRelAcc','0.001','--rAbsAcc','0.00001']
  (wd/'command.json').write_text(json.dumps(cmd,indent=2));receipt=execute(cmd,wd,'combine.log');(wd/'receipt.json').write_text(json.dumps(receipt));print(m,role,receipt,flush=True)
 execute([PY,str(S/'quality.py'),m],d,'quality.log');execute([PLOT,str(S/'report.py'),m],d,'report.log')
with concurrent.futures.ThreadPoolExecutor(max_workers=2) as ex:list(ex.map(job,[f'{m/10:.1f}' for m in range(65,86)]))
