from pathlib import Path
import json,subprocess
P=Path.cwd();R=P/'runs/20261029_gof_ftest_75';g=R/'degree2/toys_1000';pars=json.loads((P/'runs/20261029_upsilon_control/parameters.json').read_text())['joint']['b']['parameters'];coeff={x['name']:x['value'] for x in pars if x['name'].startswith('local_c')}
for idx in [727,826,956]:
 for shift in [-.02,.02,-.05,.05]:
  d=g/f'retry/toy{idx}/start{shift}';d.mkdir(exist_ok=True);c=json.loads((g/f'retry/toy{idx}/degree2/command.json').read_text());c[c.index('--setParameters')+1]='r=0,'+','.join(f'{k}={v+shift}' for k,v in coeff.items());c+=['--trackParameters','rgx{local_.*}'];(d/'command.json').write_text(json.dumps(c,indent=2))
  with (d/'combine.log').open('w') as f:subprocess.run(c,cwd=d,stdout=f,stderr=subprocess.STDOUT,check=True)
