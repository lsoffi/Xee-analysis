from pathlib import Path
import subprocess,json,shutil,os
P=Path.cwd();R=P/'runs/20261030_fixed_sr_scan';S=P/'scripts/fixed_sr_scan';PY=str(P/'bin/python');PLOT='/Users/imb/Documents/Codex/2026-10-27/referenced-chatgpt-conversation-this-is-an/work/fit-env/bin/python'
for m,roles in [('7.2',['fitdiagnostics','significance','limits']),('7.8',['fitdiagnostics'])]:
 d=R/('m'+m)/'joint';q=json.loads((d/'quality.json').read_text());params={p['name']:p['value'] for p in q['b']['parameters']};params['r']=1e-6;seed=','.join(f'{k}={v}' for k,v in params.items())
 for role in roles:
  wd=d/role;archive=d/(role+'_refine2');shutil.move(str(wd),str(archive));wd.mkdir();cmd=json.loads((archive/'command.json').read_text());cmd[cmd.index('--cminDefaultMinimizerStrategy')+1]='0';cmd[cmd.index('--cminDefaultMinimizerTolerance')+1]='0.001';cmd+=['--setParameters',seed]
  if role=='limits':cmd+=['--rMin','-1']
  (wd/'command.json').write_text(json.dumps(cmd,indent=2))
  with (wd/'combine.log').open('w') as out:res=subprocess.run(cmd,cwd=wd,stdout=out,stderr=subprocess.STDOUT,timeout=300)
  print(m,role,res.returncode,flush=True)
 for py,script,log in [(PY,'quality.py','quality_refined.log'),(PLOT,'report.py','report_refined.log')]:
  with (d.parent/log).open('w') as f:subprocess.run([py,str(S/script),m],cwd=P,stdout=f,stderr=subprocess.STDOUT,env=dict(os.environ,MPLCONFIGDIR=str(P/'work/mpl-cache')))
