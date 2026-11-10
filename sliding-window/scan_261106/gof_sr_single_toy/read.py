import ROOT as R,sys,json,math,shutil
from pathlib import Path
m,seed,out=sys.argv[1:];fs=list(Path.cwd().glob('higgsCombine*.GoodnessOfFit.*.root'));assert len(fs)==1
f=R.TFile.Open(str(fs[0]));assert f and not f.IsZombie();t=f.Get('limit');assert t
v=[float(x.limit) for x in t];assert len(v)==1 and math.isfinite(v[0]) and v[0]>=0
p=Path(out)/f'{m}_{seed}.json';tmp=p.with_suffix('.tmp');tmp.write_text(json.dumps({'mass':float(m),'seed':int(seed),'statistic':v[0]}));tmp.replace(p)
shutil.copy2(fs[0],p.with_suffix('.root'))
