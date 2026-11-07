from pathlib import Path
import json,time
import numpy as np,uproot
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
P=Path.cwd();R=P/'runs/20261030_sr6_cr4_scan';masses=[f'{k/10:.1f}' for k in range(60,91)];rows=[]
for m in masses:
 d=R/('m'+m);st=dict(mass=float(m),fit=(d/'joint/quality.json').exists(),plot=(d/'fit_SR.png').exists(),result=(d/'results.json').exists());rows.append(st)
(R/'progress.json').write_text(json.dumps(dict(updated=time.strftime('%Y-%m-%d %H:%M:%S'),points=rows),indent=2))
# Keep progress usable before the first complete limit calculation.
if not any(x['result'] for x in rows):
 html='<meta charset="utf-8"><meta http-equiv="refresh" content="60"><style>body{font:17px system-ui;max-width:1200px;margin:30px auto}img{width:100%}</style><h1>Scan 6–9 GeV: envelope Chebyshev6/Bernstein6 nella SR</h1><p>CR 9,093333–10,506667, Chebyshev4, mu e normalizzazioni Υ libere e condivise. DCB segnale in SR e CR, unico r. Risultati parziali: controlli numerici in corso.</p>'
 for r in rows:
  if r['plot']:html+=f'<h2>Massa {r["mass"]:.1f} GeV</h2><img src="m{r["mass"]:.1f}/fit_SR.png"><a href="m{r["mass"]:.1f}/fit_CR.png">Fit CR</a>'
 html+='<p>Limiti, p-value, bias Asimov e GoF con 100 toy in elaborazione.</p>';(R/'index.html').write_text(html)
# Preliminary fit shapes are written as soon as FitDiagnostics is done, independently of limits.
print('plotted',sum(x['plot'] for x in rows),'complete',sum(x['result'] for x in rows),flush=True)

if any(x['result'] for x in rows):
 import subprocess,sys
 with (R/'summary_progress.log').open('w') as out:subprocess.run([sys.executable,str(P/'scripts/sr6_cr4_scan/summary.py')],cwd=P,stdout=out,stderr=subprocess.STDOUT)
