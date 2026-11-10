import ROOT as R,sys,json,math
from pathlib import Path
m,kind,point,start,tag=sys.argv[1:];base=Path(__file__).resolve().parent;p=base/'results'/tag
fs=list(p.glob('higgsCombine*.root'));assert len(fs)==1
f=R.TFile.Open(str(fs[0]));assert f and not f.IsZombie();t=f.Get('limit');assert t
names=[x.GetName() for x in t.GetListOfBranches()];rows=[]
for row in t:
 d={n:float(getattr(row,n)) for n in names};rows.append(d)
assert rows and all(math.isfinite(x) for d in rows for x in d.values())
out=dict(mass=float(m),kind=kind,r_test=float(point),start=float(start),tag=tag,rows=rows)
ff=p/'multidimfit.root'
if ff.exists():
 rf=R.TFile.Open(str(ff));fit=rf.Get('fit_mdf')
 if fit:out['saved_fit']={'status':fit.status(),'covQual':fit.covQual(),'edm':fit.edm(),'minNll':fit.minNll()}
# Saved fit_mdf may describe the unconditional fit; do not label it as conditional convergence.
(p/'result.json').write_text(json.dumps(out,indent=2))
