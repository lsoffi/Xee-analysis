import pathlib,json,subprocess,concurrent.futures,os
P=pathlib.Path('/Users/imb/Documents/Codex/Analisi-statistica-Xee');R=P/'runs/20261029_local_gaussian_scan_7_8';pp='/Users/imb/Documents/Codex/2026-10-27/referenced-chatgpt-conversation-this-is-an/work/fit-env/bin/python'
def run(m):
 d=R/('m'+m)/'fitdiagnostics';a=d/'initial_interval_attempt';a.mkdir(exist_ok=True)
 for p in list(d.iterdir()):
  if p.is_file():p.rename(a/p.name)
 cmd=json.loads((a/'command.json').read_text())+['--stepSize','0.001','--setRobustFitStrategy','1'];(d/'command.json').write_text(json.dumps(cmd,indent=2))
 with (d/'combine.log').open('w') as f:r=subprocess.run(cmd,cwd=d,stdout=f,stderr=subprocess.STDOUT,timeout=120)
 (d/'receipt.json').write_text(json.dumps(dict(role='fitdiagnostics',exit_code=r.returncode,refined_profile_step=True)))
 for s,py in [('quality',str(P/'bin/python')),('report',pp)]:
  with (d.parent/(s+'.log')).open('w') as f:subprocess.run([py,str(P/f'scripts/local_gaussian_scan/{s}.py'),m],stdout=f,stderr=subprocess.STDOUT,env=dict(os.environ,MPLCONFIGDIR=str(P/'work/mpl-cache')),check=True)
 print(m,'refined',flush=True)
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as e:list(e.map(run,['7.0','7.1','7.2','7.3']))
