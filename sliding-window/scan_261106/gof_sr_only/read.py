import sys,json,math,ROOT
from pathlib import Path
m,kind,out=sys.argv[1:]
files=list(Path.cwd().glob('higgsCombine*.GoodnessOfFit.*.root'));assert len(files)==1
f=ROOT.TFile.Open(str(files[0]));assert f and not f.IsZombie()
t=f.Get('limit');assert t
values=[float(r.limit) for r in t]
assert len(values)==(500 if kind=='toys' else 1),len(values)
assert all(math.isfinite(v) and v>=0 for v in values)
p=Path(out)/f'{m}_{kind}.tmp';p.write_text(json.dumps(dict(mass=float(m),kind=kind,statistics=values)));p.rename(p.with_suffix('.json'))
import shutil
shutil.copy2(files[0],Path(out)/f'{m}_{kind}.root')
