import sys,json
from pathlib import Path
src,seed,g=sys.argv[1:]
choices=[]
for sub in ('results','retry_results'):
 p=Path(src)/sub/seed/'result.json'
 if p.exists():
  o=json.loads(p.read_text());assert o['seed']==int(seed) and o['generator']==g
  b=o['branches']['L3']
  if b['ok']: choices.append(b)
if not choices:sys.exit(1)
Path('L3.json').write_text(json.dumps(min(choices,key=lambda b:b['score'])))
