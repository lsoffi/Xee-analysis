from pathlib import Path
import json,subprocess,shutil,os
P=Path.cwd();R=P/'runs/20261030_sr7_cr3_scan';S=P/'scripts/sr7_cr3_scan';PLOT='/Users/imb/Documents/Codex/2026-10-27/referenced-chatgpt-conversation-this-is-an/work/fit-env/bin/python'
for result in sorted(R.glob('m*/results.json')):
 a=json.loads(result.read_text());d=result.parent/'joint'
 if a['significance_consistent'] or not a['fit_ok'] or (d/'significance_initial').exists():continue
 q=json.loads((d/'quality.json').read_text());seed=','.join(f"{p['name']}={p['value']}" for p in q['s']['parameters']);wd=d/'significance';old=d/'significance_initial';shutil.move(wd,old);wd.mkdir();cmd=json.loads((old/'command.json').read_text());cmd+=['--setParameters',seed];cmd[cmd.index('--cminDefaultMinimizerStrategy')+1]='0';(wd/'command.json').write_text(json.dumps(cmd,indent=2))
 with (wd/'combine.log').open('w') as f:res=subprocess.run(cmd,cwd=wd,stdout=f,stderr=subprocess.STDOUT,timeout=120)
 (wd/'receipt.json').write_text(json.dumps(dict(exit_code=res.returncode,reason='initialize from converged S+B fit to check local minimum')))
 with (result.parent/'report_refined.log').open('w') as f:subprocess.run([PLOT,str(S/'report.py'),f"{a['mass']:.1f}",'numeric'],cwd=P,stdout=f,stderr=subprocess.STDOUT,env=dict(os.environ,MPLCONFIGDIR=str(P/'work/mpl-cache')))
 b=json.loads(result.read_text());print(a['mass'],a['Z'],b['Z'],b['Z_from_FD'],b['numerical_ok'],flush=True)
