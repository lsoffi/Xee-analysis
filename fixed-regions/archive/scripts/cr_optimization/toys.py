from study import *
from scipy.special import betainc
from concurrent.futures import ProcessPoolExecutor
import os
CAND=[(232,338),(240,330),(240,323),(247,323),(255,315)]
N=200

def bernstein_generator(datasets=None,save=True):
 lo=6+232/75;hi=6+338/75;fits=[]
 for era in range(4):
  A,n,mask=matrices(lo,hi,4,era);n=n if datasets is None else np.asarray(datasets[era]);ix=np.flatnonzero(mask);u=(edges[np.r_[ix,ix[-1]+1]]-lo)/(hi-lo);B=np.array([np.diff(betainc(j+1,5-j,u))/5 for j in range(5)]).T;A=np.column_stack([B,A[:,-2:]])
  def cost(*p):
   mu=A@p
   return 2*np.sum(mu-n+xlogy(n,n/mu)) if min(mu)>0 else 1e30
  def grad(*p):return 2*A.T@(1-n/(A@p))
  m=Minuit(cost,*([max(1,n.sum()-A[:,-2:].sum())]*5+[1.,1.]),grad=grad);m.errordef=1
  for j in range(5):m.limits[j]=(1e-6,None)
  for j in [5,6]:m.limits[j]=(0,10)
  m.tol=1e-5;m.strategy=2;m.migrad();m.hesse();assert m.valid
  fits.append(dict(p=list(m.values),mu=(A@np.array(m.values)).tolist(),D=m.fval))
 g=dict(lo=lo,hi=hi,degree=4,family='positive Bernstein4',fits=fits,D=sum(f['D'] for f in fits));g['p_asym']=float(chi2.sf(g['D'],4*(106-7)))
 if save:(R/'bernstein_generator.json').write_text(json.dumps(g))
 return g

def task(args):
 typ,idx=args;rng=np.random.default_rng(987120+idx+(0 if typ=='bias' else 100000));scan=json.loads((R/'scan.json').read_text());result=[]
 if typ=='bias':
  g=json.loads((R/'bernstein_generator.json').read_text());data=[rng.poisson(f['mu']) for f in g['fits']]
  truth=[]
  for a,b in [(6,9.093333333333334),(10.506666666666668,12)]:
   mask=(edges[:-1]>=a-1e-8)&(edges[1:]<=b+1e-8);truth.append(sum(np.array(ER[i]['templates'])[:,mask].sum(axis=1)@np.array(g['fits'][i]['p'][-2:]) for i in range(4)))
  for li,ui in CAND:
   lo,hi=6+li/75,6+ui/75;d=[x[li-232:ui-232] for x in data]
   for degree in [1,2,3]:
    obs=next(r for r in scan if abs(r['lo']-lo)<1e-7 and abs(r['hi']-hi)<1e-7 and r['degree']==degree);r=combined(lo,hi,degree,d,[f['p'] for f in obs['fits']]);result.append(dict(kind=typ,toy=idx,li=li,ui=ui,degree=degree,valid=r['valid'],prediction=r['prediction'],error=r['error'],truth=truth,D=r['D']))
 else:
  # idx enumerates a toy index; generate independently for each null degree and window.
  for li,ui in CAND:
   lo,hi=6+li/75,6+ui/75
   for degree in [1,2,3]:
    obs=next(r for r in scan if abs(r['lo']-lo)<1e-7 and abs(r['hi']-hi)<1e-7 and r['degree']==degree);d=[rng.poisson(f['mu']) for f in obs['fits']];r=combined(lo,hi,degree,d,[f['p'] for f in obs['fits']]);rec=dict(kind=typ,toy=idx,li=li,ui=ui,degree=degree,valid=r['valid'],D=r['D'])
    if degree<3:
     seeds=[f['p'][:degree+1]+[0.]+f['p'][-2:] for f in r['fits']];r2=combined(lo,hi,degree+1,d,seeds);rec.update(delta=r['D']-r2['D'],valid_next=r2['valid'])
    result.append(rec)
 return result
if __name__=='__main__':
 g=bernstein_generator();print('Bernstein4 p asym',g['p_asym'],flush=True)
 jobs=[(typ,i) for typ in ['gof','bias'] for i in range(N)]
 with ProcessPoolExecutor(max_workers=3) as pool, (R/'toys.jsonl').open('w') as out:
  for j,res in enumerate(pool.map(task,jobs)):
   for row in res:out.write(json.dumps(row)+'\n')
   out.flush()
   if j%10==9:print('completed',j+1,'/',len(jobs),flush=True)
 print('DONE',flush=True)
