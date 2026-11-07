import pathlib,subprocess,concurrent.futures,os,json
P=pathlib.Path('/Users/imb/Documents/Codex/Analisi-statistica-Xee');R=P/'runs/20261029_local_dcb_scan_7_8';plotpython='/Users/imb/Documents/Codex/2026-10-27/referenced-chatgpt-conversation-this-is-an/work/fit-env/bin/python'
def run(m):
 d=R/('m'+m);d.mkdir(exist_ok=True)
 for step in ['prepare','run','quality','report']:
  py=plotpython if step=='report' else str(P/'bin/python') if step in ['prepare','quality'] else 'python3'
  with (d/(step+'.log')).open('w') as f:
   cp=subprocess.run([py,str(P/f'scripts/local_dcb_scan/{step}.py'),m],cwd=P,stdout=f,stderr=subprocess.STDOUT,env=dict(os.environ,MPLCONFIGDIR=str(P/'work/mpl-cache')),timeout=700)
  if cp.returncode:print(m,step,'FAILED',flush=True);return dict(mass=m,failed=step)
 print(m,'COMPLETE',flush=True);return dict(mass=m,complete=True)
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as e:results=list(e.map(run,[f'{m/10:.1f}' for m in range(70,81)]))
(R/'pipeline_status.json').write_text(json.dumps(results,indent=2))
