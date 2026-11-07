from pathlib import Path
import json, numpy as np, time
from numpy.polynomial import chebyshev as C
from iminuit import Minuit
from scipy.special import xlogy
from scipy.stats import chi2
P=Path.cwd();R=P/'runs/20261030_cr_optimization';I=json.loads((R/'input.json').read_text());edges=np.array(I['edges']);ER=I['eras']
def matrices(lo,hi,deg,era):
 mask=(edges[:-1]>=lo-1e-8)&(edges[1:]<=hi+1e-8); ix=np.flatnonzero(mask); a,b=edges[ix[0]],edges[ix[-1]+1];te=2*(edges[np.r_[ix,ix[-1]+1]]-a)/(b-a)-1
 bases=[]
 for j in range(deg+1):
  z=np.zeros(j+1);z[-1]=1;integ=C.chebint(z);bases.append(np.diff(C.chebval(te,integ))/2)
 B=np.array(bases).T;Y=np.array(ER[era]['templates'])[:,mask].T
 return np.column_stack([B,Y]),np.array(ER[era]['counts'])[mask],mask

def fit(lo,hi,deg,era,data=None,seed=None):
 A,n,mask=matrices(lo,hi,deg,era);n=n if data is None else np.asarray(data);N=n.sum();T=C.chebvander(np.linspace(-1,1,129),deg);nt=deg+1
 def fun(*p):
  p=np.array(p);mu=A@p;pol=T@p[:nt]
  if np.min(mu)<=0 or np.min(pol)<=0:return 1e30+1e12*np.sum(np.minimum(pol,0)**2)+1e12*np.sum(np.minimum(mu,0)**2)
  return 2*np.sum(mu-n+xlogy(n,n/mu))
 def grad(*p):
  p=np.array(p);mu=A@p
  if np.min(mu)<=0 or np.min(T@p[:nt])<=0:return np.zeros(len(p))
  return 2*A.T@(1-n/mu)
 initial=np.r_[max(1,N-A[:,-2:].sum()),np.zeros(deg),1.,1.]
 starts=[initial] if seed is None else [np.array(seed),initial]
 best=None
 for start in starts:
  m=Minuit(fun,*start,grad=grad);m.errordef=1;m.limits[0]=(1e-5,None);m.limits[nt]=(0,10);m.limits[nt+1]=(0,10);m.errors=tuple([max(1,N*.01)]*nt+[.1,.1]);m.strategy=1;m.tol=1e-5;m.migrad(ncall=5000);m.hesse()
  if best is None or m.fval<best.fval:best=m
 m=best;p=np.array(m.values);cov=np.array(m.covariance) if m.covariance is not None else np.zeros((len(p),len(p)))
 roots=C.chebroots(C.chebder(p[:nt]));rr=[float(z.real) for z in roots if abs(z.imag)<1e-7 and -1<z.real<1];pos=min(C.chebval(np.r_[-1,1,rr],p[:nt]));mu=A@p
 gradients=[]
 for a,b in [(6,9.093333333333334),(10.506666666666668,12)]:
  sel=(edges[:-1]>=a-1e-8)&(edges[1:]<=b+1e-8);g=np.r_[np.zeros(nt),np.array(ER[era]['templates'])[:,sel].sum(axis=1)];gradients.append(g)
 G=np.array(gradients);pred=G@p;pcov=G@cov@G.T
 corr=cov/np.sqrt(np.maximum(np.outer(np.diag(cov),np.diag(cov)),1e-200))
 return dict(p=p.tolist(),cov=cov.tolist(),D=float(m.fval),valid=bool(m.valid and m.fmin.has_posdef_covar and pos>0),accurate=bool(m.fmin.has_accurate_covar),edm=float(m.fmin.edm),minpoly=float(pos),mu=mu.tolist(),prediction=pred.tolist(),prediction_cov=pcov.tolist(),boundary=[bool(p[-2]<1e-4 or p[-2]>9.9999),bool(p[-1]<1e-4 or p[-1]>9.9999)],max_Y_cont_corr=float(np.max(abs(corr[-2:,:nt]))),Ycorr=float(corr[-1,-2]),fraction=np.array(ER[era]['fractions'])[:,mask].sum(axis=1).tolist())

def combined(lo,hi,deg,datasets=None,seeds=None):
 fits=[fit(lo,hi,deg,i,None if datasets is None else datasets[i],None if seeds is None else seeds[i]) for i in range(4)];pred=np.sum([f['prediction'] for f in fits],axis=0);cov=np.sum([f['prediction_cov'] for f in fits],axis=0);D=sum(f['D'] for f in fits);nb=len(fits[0]['mu']);return dict(lo=lo,hi=hi,width=hi-lo,degree=deg,fits=fits,D=D,dof=4*(nb-deg-3),p_asym=float(chi2.sf(D,4*(nb-deg-3))),valid=all(f['valid'] for f in fits),prediction=pred.tolist(),error=np.sqrt(np.diag(cov)).tolist(),boundary_count=sum(sum(f['boundary']) for f in fits),maxcorr=max(f['max_Y_cont_corr'] for f in fits))

def scan():
 out=[]
 for li in [232,240,247,255]:
  for ui in [338,330,323,315]:
   for deg in [1,2,3]:
    lo=6+li/75;hi=6+ui/75;rec=combined(lo,hi,deg);out.append(rec);(R/'scan.json').write_text(json.dumps(out));print(lo,hi,deg,rec['valid'],round(rec['p_asym'],4),np.round(rec['error'],1),flush=True)
 return out
if __name__=='__main__':scan()
