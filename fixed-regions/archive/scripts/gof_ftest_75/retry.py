from pathlib import Path
import json,subprocess
P=Path.cwd();R=P/'runs/20261029_gof_ftest_75'
ids=json.loads((R/'audit_toy_ids.json').read_text());g=R/'degree2/toys_1000'
for idx in ids:
 d=g/f'retry/toy{idx}';d.mkdir(parents=True,exist_ok=True);tf=d/'single.root';src=next(g.glob('higgsCombine*.root'));subprocess.run([str(P/'bin/python'),str(P/'scripts/validation/extract_single_toy.py'),str(src),str(idx),str(tf)],check=True,stdout=subprocess.DEVNULL)
 for n in [2,3]:
  dd=d/f'degree{n}';dd.mkdir(exist_ok=True);c=json.loads((g/f'fit_degree{n}/command.json').read_text());c[c.index('--toysFile')+1]=str(tf);c[c.index('-t')+1]='1';c[c.index('--cminDefaultMinimizerStrategy')+1]='2';c[c.index('--cminDefaultMinimizerTolerance')+1]='0.00001';(dd/'command.json').write_text(json.dumps(c,indent=2))
  with (dd/'combine.log').open('w') as f:subprocess.run(c,cwd=dd,stdout=f,stderr=subprocess.STDOUT,check=True)
 print(idx,flush=True)
