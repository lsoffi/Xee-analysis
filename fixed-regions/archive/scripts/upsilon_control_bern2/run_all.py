from pathlib import Path
import subprocess,json,os
P=Path.cwd();plotpython='/Users/imb/Documents/Codex/2026-10-27/referenced-chatgpt-conversation-this-is-an/work/fit-env/bin/python'
for n in [2,3]:
 S=P/f'scripts/upsilon_control_bern{n}';R=P/f'runs/20261029_upsilon_control_bern{n}'
 for script,args in [('prepare_signal.py',['7.5']),('prepare_control.py',['7.5']),('build_joint.py',[])]:
  with (R/(script+'.log')).open('w') as f:subprocess.run([str(P/'bin/python'),str(S/script)]+args,stdout=f,stderr=subprocess.STDOUT,check=True)
 d=R/'control/fitdiagnostics';d.mkdir(exist_ok=True)
 c=json.loads((P/'runs/20261029_upsilon_control_cheb3/control/fitdiagnostics/command.json').read_text());c[c.index('-d')+1]=str(R/'control/workspace.root');(d/'command.json').write_text(json.dumps(c,indent=2))
 with (d/'combine.log').open('w') as f:subprocess.run(c,cwd=d,stdout=f,stderr=subprocess.STDOUT,check=True)
 subprocess.run(['python3',str(S/'run_joint.py'),'7.5'],check=True)
 with (R/'extract.log').open('w') as f:subprocess.run([str(P/'bin/python'),str(S/'extract.py')],stdout=f,stderr=subprocess.STDOUT,check=True)
 subprocess.run([plotpython,str(S/'report.py')],env=dict(os.environ,MPLCONFIGDIR=str(P/'work/mpl-cache')),check=True)
