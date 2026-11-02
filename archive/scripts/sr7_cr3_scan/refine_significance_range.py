from pathlib import Path
import json,subprocess,shutil,os,concurrent.futures
P=Path.cwd();R=P/'runs/20261030_sr7_cr3_scan';S=P/'scripts/sr7_cr3_scan';PLOT='/Users/imb/Documents/Codex/2026-10-27/referenced-chatgpt-conversation-this-is-an/work/fit-env/bin/python'
def job(result):
 a=json.loads(result.read_text());d=result.parent/'joint';wd=d/'significance';old=d/'significance_wide';shutil.move(wd,old);wd.mkdir();cmd=json.loads((old/'command.json').read_text());cmd[cmd.index('--rMax')+1]=str(max(1,10*abs(a['r']),5*a['limits']['expected_50']));cmd[cmd.index('--rMin')+1]='-1';cmd[cmd.index('--cminDefaultMinimizerStrategy')+1]='2';cmd[cmd.index('--cminDefaultMinimizerTolerance')+1]='0.00001';(wd/'command.json').write_text(json.dumps(cmd,indent=2))
 with (wd/'combine.log').open('w') as out:res=subprocess.run(cmd,cwd=wd,stdout=out,stderr=subprocess.STDOUT,timeout=120)
 (wd/'receipt.json').write_text(json.dumps(dict(exit_code=res.returncode,reason='uniform scale-aware POI domain to resolve boundary minimization; bounds enclose unconstrained FD fit and expected limit')))
 with (result.parent/'report_sig_range.log').open('w') as out:subprocess.run([PLOT,str(S/'report.py'),f"{a['mass']:.1f}",'numeric'],cwd=P,stdout=out,stderr=subprocess.STDOUT,env=dict(os.environ,MPLCONFIGDIR=str(P/'work/mpl-cache')))
 b=json.loads(result.read_text());print(a['mass'],a['Z'],b['Z'],b['Z_from_FD'],b['numerical_ok'],flush=True)
with concurrent.futures.ThreadPoolExecutor(max_workers=2) as ex:list(ex.map(job,sorted(R.glob('m*/results.json'))))
