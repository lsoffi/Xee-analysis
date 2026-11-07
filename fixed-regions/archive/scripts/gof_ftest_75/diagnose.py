from pathlib import Path
import json,subprocess
P=Path.cwd();R=P/'runs/20261029_gof_ftest_75';g=R/'degree2/toys_1000'
for idx in [727,826,956]:
 d=g/f'retry/toy{idx}/fd';d.mkdir(exist_ok=True);c=json.loads((g/f'retry/toy{idx}/degree2/command.json').read_text());c[c.index('-M')+1]='FitDiagnostics'
 for opt in ['--algo','--fixedSignalStrength']:
  j=c.index(opt);del c[j:j+2]
 c+=['--skipSBFit','--minos','none','--robustFit','0','--saveShapes'];(d/'command.json').write_text(json.dumps(c,indent=2))
 with (d/'combine.log').open('w') as f:subprocess.run(c,cwd=d,stdout=f,stderr=subprocess.STDOUT,check=True)
