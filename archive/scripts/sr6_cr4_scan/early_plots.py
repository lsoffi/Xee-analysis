from pathlib import Path
import subprocess,time,json,os
P=Path.cwd();R=P/'runs/20261030_sr6_cr4_scan';S=P/'scripts/sr6_cr4_scan';PY=str(P/'bin/python');PLOT='/Users/imb/Documents/Codex/2026-10-27/referenced-chatgpt-conversation-this-is-an/work/fit-env/bin/python'
for iteration in range(180):
 for m in [f'{k/10:.1f}' for k in range(60,91)]:
  d=R/('m'+m)
  if (d/'fit_SR.png').exists() or not (d/'joint/fitdiagnostics/combine.log.receipt.json').exists():continue
  for script,interp in [('quality.py',PY),('plot_fit.py',PLOT)]:
   with (d/('early_'+script+'.log')).open('w') as out:subprocess.run([interp,str(S/script),m],cwd=P,stdout=out,stderr=subprocess.STDOUT,env=dict(os.environ,MPLCONFIGDIR=str(P/'work/mpl-cache')))
 subprocess.run([PLOT,str(S/'progress.py')],cwd=P,env=dict(os.environ,MPLCONFIGDIR=str(P/'work/mpl-cache')))
 if len(list(R.glob('m*/fit_SR.png')))==31:break
 time.sleep(20)
