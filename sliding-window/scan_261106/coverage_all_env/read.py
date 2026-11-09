import sys,json,math,ROOT
from pathlib import Path
seed,gen,inj,out=sys.argv[1:]
f=ROOT.TFile.Open(str(next(Path.cwd().glob('higgsCombine*.AsymptoticLimits.*.root'))))
t=f.Get('limit');values=[float(row.limit) for row in t if row.quantileExpected<0]
assert len(values)==1 and math.isfinite(values[0]) and 0<=values[0]<=1,values
p=Path(out)/(seed+'.tmp');p.write_text(json.dumps(dict(seed=int(seed),generator=gen,injection=float(inj),CLs=values[0],excluded=values[0]<0.05)));p.rename(p.with_suffix('.json'))
