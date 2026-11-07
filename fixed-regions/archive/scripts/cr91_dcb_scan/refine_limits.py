from pathlib import Path
import json,subprocess,time,concurrent.futures,os
P=Path.cwd();R=P/'runs/20261029_cr91_dcb_scan';plot='/Users/imb/Documents/Codex/2026-10-27/referenced-chatgpt-conversation-this-is-an/work/fit-env/bin/python'
def job(p):
 row=json.loads(p.read_text());m=f'{row["mass"]:.1f}';J=p.parent/'joint';old=J/'limits';archive=J/'initial_limits';old.rename(archive);old.mkdir();c=json.loads((archive/'command.json').read_text());c[c.index('--cminDefaultMinimizerStrategy')+1]='2';c[c.index('--cminDefaultMinimizerTolerance')+1]='0.00001';c[c.index('--rMax')+1]=str(max(1,4*row['limits']['expected_50']));(old/'command.json').write_text(json.dumps(c,indent=2));t=time.time()
 with (old/'combine.log').open('w') as f:res=subprocess.run(c,cwd=old,stdout=f,stderr=subprocess.STDOUT)
 (old/'receipt.json').write_text(json.dumps(dict(exit_code=res.returncode,seconds=time.time()-t)));assert res.returncode==0
 with (p.parent/'report_refined.log').open('w') as f:subprocess.run([plot,str(P/'scripts/cr91_dcb_scan/report.py'),m,'numeric'],stdout=f,stderr=subprocess.STDOUT,check=True,env=dict(os.environ,MPLCONFIGDIR=str(P/'work/mpl-cache')))
 row=json.loads(p.read_text());print(m,row['limits']['observed'],row['numerical_ok'],flush=True)
with concurrent.futures.ThreadPoolExecutor(max_workers=2) as ex:list(ex.map(job,sorted(R.glob('m*/results.json'))))
