from pathlib import Path
import subprocess,time,os,json
P=Path.cwd();S=P/'scripts/sr6_cr4_scan';R=P/'runs/20261030_sr6_cr4_scan';plot='/Users/imb/Documents/Codex/2026-10-27/referenced-chatgpt-conversation-this-is-an/work/fit-env/bin/python'
for it in range(180):
 subprocess.run([str(P/'bin/python'),str(S/'collect_bias.py')],cwd=P)
 subprocess.run([plot,str(S/'plot_bias.py')],cwd=P,env=dict(os.environ,MPLCONFIGDIR=str(P/'work/mpl-cache')))
 if len(json.loads((R/'bias_summary.json').read_text()))==124:break
 time.sleep(25)
