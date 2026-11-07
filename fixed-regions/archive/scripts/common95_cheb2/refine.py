import sys,pathlib,json,subprocess,os,time,signal,shutil,concurrent.futures
P=pathlib.Path('/Users/imb/Documents/Codex/Analisi-statistica-Xee');R=P/'runs/20261029_common95_cheb2';PY='/Users/imb/Documents/Codex/2026-10-27/referenced-chatgpt-conversation-this-is-an/work/fit-env/bin/python';ENV=dict(os.environ,MPLCONFIGDIR=str(P/'work/mpl-cache'),OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1')
def run(m):
 d=R/('m'+m);result=json.loads((d/'results.json').read_text());exp=result.get('limits',{}).get('expected_50',1);bound=max(1,min(1000,10*exp));archive=d/'initial_attempt';archive.mkdir(exist_ok=True)
 roles=[]
 if not result['fit_ok']:roles+=['fitdiagnostics','significance']
 if not result.get('limit_crossing',{}).get('verified'):roles+=['limits']
 for role in roles:
  rd=d/role;orig=json.loads((rd/'command.json').read_text());cmd=orig[:];cmd[cmd.index('--cminDefaultMinimizerStrategy')+1]='2';cmd[cmd.index('--cminDefaultMinimizerTolerance')+1]='0.0001'
  cmd[cmd.index('--rMax')+1]=str(bound)
  if '--rMin' in cmd:cmd[cmd.index('--rMin')+1]=str(-bound)
  if (archive/role).exists():continue
  shutil.move(str(rd),str(archive/role));rd.mkdir();(rd/'command.json').write_text(json.dumps(cmd,indent=2));start=time.time()
  with (rd/'combine.log').open('w') as log:
   proc=subprocess.Popen(cmd,cwd=rd,env=ENV,stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
   try:rc=proc.wait(timeout=600)
   except subprocess.TimeoutExpired:os.killpg(proc.pid,signal.SIGKILL);proc.wait();rc=124
  (rd/'receipt.json').write_text(json.dumps(dict(exit_code=rc,elapsed=time.time()-start,role=role,refined=True),indent=2))
 for step in ['quality','report']:
  if step=='quality':
   for f in ['quality.json','parameters.csv','results.json','fit.png','fit.pdf']:
    if (d/f).exists() and not (archive/f).exists():shutil.copy2(d/f,archive/f)
  with (d/(step+'.log')).open('w') as log:
   cp=subprocess.run([str(P/'bin/python') if step=='quality' else PY,str(P/f'scripts/common95_cheb2/{step}.py'),m],cwd=P,env=ENV,stdout=log,stderr=subprocess.STDOUT,timeout=120)
  if cp.returncode:print(m,'FAILED',step,flush=True);return
 r=json.loads((d/'results.json').read_text());print(m,'REFINED','fit_ok',r['fit_ok'],'CLs',r.get('limit_crossing',{}).get('verified'),flush=True)
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as e:list(e.map(run,sys.argv[1:]))
subprocess.run([PY,str(P/'scripts/common95_cheb2/summary.py')],env=ENV)
