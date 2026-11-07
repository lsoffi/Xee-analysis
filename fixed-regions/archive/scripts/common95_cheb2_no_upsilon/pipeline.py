import pathlib,subprocess,concurrent.futures,os,json,sys,time
P=pathlib.Path('/Users/imb/Documents/Codex/Analisi-statistica-Xee');R=P/'runs/20261029_common95_cheb2_no_upsilon';PY='/Users/imb/Documents/Codex/2026-10-27/referenced-chatgpt-conversation-this-is-an/work/fit-env/bin/python';ENV=dict(os.environ,MPLCONFIGDIR=str(P/'work/mpl-cache'),OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1')
def run(m):
 d=R/('m'+m);d.mkdir(parents=True,exist_ok=True)
 for step in ['prepare','run','quality','report']:
  if step=='prepare' and (d/'workspace.root').exists():continue
  py=PY if step=='report' else str(P/'bin/python') if step in ['prepare','quality'] else '/usr/bin/python3'
  with (d/(step+'.log')).open('w') as f:
   try:cp=subprocess.run([py,str(P/f'scripts/common95_cheb2_no_upsilon/{step}.py'),m],cwd=P,stdout=f,stderr=subprocess.STDOUT,env=ENV,timeout=2800)
   except Exception as e:return dict(mass=m,failed=step,error=str(e))
  if cp.returncode:print(m,step,'FAILED',flush=True);return dict(mass=m,failed=step)
 print(m,'COMPLETE',flush=True);return dict(mass=m,complete=True)
masses=sys.argv[1:] or [f'{m/10:.1f}' for m in range(65,101)];results=[]
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as e:
 futures={e.submit(run,m):m for m in masses}
 for f in concurrent.futures.as_completed(futures):
  results.append(f.result());(R/'pipeline_status.json').write_text(json.dumps(results,indent=2))
  subprocess.run([PY,str(P/'scripts/common95_cheb2_no_upsilon/summary.py')],env=ENV,stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)
print('FINISHED',json.dumps(results),flush=True)
