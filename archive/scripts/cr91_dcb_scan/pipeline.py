from pathlib import Path
import subprocess,json,os,concurrent.futures,time
P=Path.cwd();R=P/'runs/20261029_cr91_dcb_scan';S=P/'scripts/cr91_dcb_scan';plot='/Users/imb/Documents/Codex/2026-10-27/referenced-chatgpt-conversation-this-is-an/work/fit-env/bin/python';masses=[f'{i/10:.1f}' for i in range(65,86)]
def run(c,log):
 with log.open('w') as f:subprocess.run(c,stdout=f,stderr=subprocess.STDOUT,check=True,env=dict(os.environ,MPLCONFIGDIR=str(P/'work/mpl-cache'),OMP_NUM_THREADS='1'))
def job(m):
 d=R/('m'+m);d.mkdir(exist_ok=True)
 for region in ['signal','control']:
  if not(d/region/'workspace.root').exists():run([str(P/'bin/python'),str(S/'prepare.py'),m,region],d/f'prepare_{region}.log')
 if not(d/'joint/workspace.root').exists():run([str(P/'bin/python'),str(S/'build_joint.py'),m],d/'build_joint.log')
 run(['python3',str(S/'run.py'),m],d/'run.log');run([str(P/'bin/python'),str(S/'quality.py'),m],d/'quality.log');run([plot,str(S/'report.py'),m],d/'report.log');print(m,'DONE',flush=True)
with concurrent.futures.ThreadPoolExecutor(max_workers=2) as ex:list(ex.map(job,masses))
