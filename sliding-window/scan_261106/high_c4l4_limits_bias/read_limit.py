import sys,json,math,shutil,ROOT
from pathlib import Path
m,out=sys.argv[1:]
files=list(Path.cwd().glob('higgsCombine*.AsymptoticLimits.*.root'));assert len(files)==1
f=ROOT.TFile.Open(str(files[0]));assert f and not f.IsZombie()
t=f.Get('limit');assert t
rows=[dict(quantile=float(r.quantileExpected),limit=float(r.limit)) for r in t]
assert len(rows)==6 and all(math.isfinite(r['limit']) and r['limit']>=0 for r in rows)
for q in (-1,.025,.16,.5,.84,.975):assert sum(abs(r['quantile']-q)<.001 for r in rows)==1
expected=sorted((r for r in rows if r['quantile']>=0),key=lambda r:r['quantile'])
assert all(a['limit']<=b['limit'] for a,b in zip(expected,expected[1:]))
p=Path(out)/f'm{m}.tmp';p.write_text(json.dumps(dict(mass=float(m),quantiles=rows)));p.rename(p.with_suffix('.json'))
shutil.copy2(files[0],Path(out)/f'm{m}.root')
