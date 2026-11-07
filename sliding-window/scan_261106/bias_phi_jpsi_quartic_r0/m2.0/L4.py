import sys,json,math
from pathlib import Path
import numpy as np,ROOT
ROOT.gSystem.Load('libHiggsAnalysisCombinedLimit')
run,toyfile,g,seed=sys.argv[1:];run=Path(run)
tf=ROOT.TFile.Open(toyfile);data=tf.Get('toys/toy_1');snap=tf.Get('toys/toy_1_snapshot')
aux={p.GetName():p.getVal() for p in snap};assert len(aux)==10
obs={}
for i in range(data.numEntries()):
 row=data.get(i);key=row.getCatLabel('CMS_channel');ib=int(row.getRealValue('CMS_th1x'));obs[key,ib]=obs.get((key,ib),0)+data.weight()
assert len(obs)>100 and sum(obs.values())==data.sumEntries()
out={}
for name,degree in [('L4', 4)]:
 f=ROOT.TFile.Open(str(next(Path(name).glob('fitDiagnostics*.root'))));fit=f.Get('fit_s')
 fw=ROOT.TFile.Open(str(run/(name+'_2.0.root')));w=fw.Get('w');md=json.loads((run/(name+'.json')).read_text())
 if not fit:
  out[name]=dict(r=None,D=None,penalty=4*degree,score=None,status=-999,cov=-1,edm=None,ok=False)
  f.Close();fw.Close();continue
 pars=fit.floatParsFinal()
 for p in pars:w.var(p.GetName()).setVal(p.getVal())
 for k,v in aux.items():w.var(k).setVal(v)
 dev=0;ntot=0;valid=True
 for ch in md['channels']:
  y=np.array([obs.get((ch['label'],i),0) for i in range(len(ch['counts']))]);ntot+=sum(y)
  v=np.zeros(len(y))
  for proc in ch['processes']:
   nodes=np.array([w.arg(k).getVal() for k in proc['bin_nodes']]);v+=nodes*(w.arg(proc['coefficient']).getVal()/sum(nodes))
  if np.min(v)<=0:valid=False;break
  idx=y>0;dev+=2*(sum(v-y)+sum(y[idx]*np.log(y[idx]/v[idx])))
 if valid:assert ntot==data.sumEntries()
 for k in md['constraints']:
  u,mu,s=list(w.pdf(k).servers());dev+=((u.getVal()-mu.getVal())/s.getVal())**2
 r=pars.find('r').getVal();penalty=4*degree
 out[name]=dict(r=r,D=dev,penalty=penalty,score=dev+penalty,status=fit.status(),cov=fit.covQual(),edm=fit.edm(),ok=valid and fit.status()==0 and abs(r)<79.999 and math.isfinite(dev))
 f.Close();fw.Close()
print('SCORE',json.dumps(dict(seed=int(seed),generator=g,branches=out)))
