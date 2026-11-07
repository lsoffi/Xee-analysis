import pathlib,subprocess,concurrent.futures,os,json
P=pathlib.Path('/Users/imb/Documents/Codex/Analisi-statistica-Xee');R=P/'runs/20261029_local_poly2_scan_7_8';plotpython='/Users/imb/Documents/Codex/2026-10-27/referenced-chatgpt-conversation-this-is-an/work/fit-env/bin/python'
def run(task):
 model,m=task
 d=R/model/('m'+m);d.mkdir(parents=True,exist_ok=True)
 for step in ['prepare','run','quality','report']:
  py=plotpython if step=='report' else str(P/'bin/python') if step in ['prepare','quality'] else 'python3'
  with (d/(step+'.log')).open('w') as f:
   cp=subprocess.run([py,str(P/f'scripts/local_poly2_scan/{step}.py'),m,model],cwd=P,stdout=f,stderr=subprocess.STDOUT,env=dict(os.environ,MPLCONFIGDIR=str(P/'work/mpl-cache')),timeout=700)
  if cp.returncode:print(model,m,step,'FAILED',flush=True);return dict(model=model,mass=m,failed=step)
 print(model,m,'COMPLETE',flush=True);return dict(model=model,mass=m,complete=True)
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as e:results=list(e.map(run,[(model,f'{m/10:.1f}') for model in ['bernstein2','chebyshev2'] for m in range(70,81)]))
(R/'pipeline_status.json').write_text(json.dumps(results,indent=2))
