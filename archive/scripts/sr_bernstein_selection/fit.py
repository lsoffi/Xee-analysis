from pathlib import Path
import numpy as np,json,sys
from scipy.special import erf,xlogy,betainc
from math import comb
from numpy.polynomial import Polynomial
from numpy.polynomial import chebyshev as C
from iminuit import Minuit
from concurrent.futures import ProcessPoolExecutor
P=Path.cwd();R=P/'runs/20261030_sr_bernstein_selection';I=json.loads((P/'runs/20261030_cr_optimization/input.json').read_text());YS=json.loads((P/'runs/20261029_upsilon_peak_windows/parameters.json').read_text())['upsilon'];CR=json.loads((P/'runs/20261030_cr_degree4_mu/fits.json').read_text());edges=np.array(I['edges'])[:339];eras=[e['era'] for e in I['eras']];NS=232;NC=106

def primitive(x,mu,sig,q):
 t=(np.asarray(x)-mu)/sig;v=np.sqrt(np.pi/2)*erf(t/np.sqrt(2))
 for side,sgn in [('L',-1),('R',1)]:
  a=abs(q['alpha'+side]);n=q['n'+side];z=sgn*t;sel=z>a;B=n/a-a;A=(n/a)**n*np.exp(-a*a/2)
  tail=A*np.log((B+z[sel])/(B+a)) if abs(n-1)<1e-8 else A*((B+z[sel])**(1-n)-(B+a)**(1-n))/(1-n)
  v[sel]=sgn*(np.sqrt(np.pi/2)*erf(a/np.sqrt(2))+tail)
 return v*sig

def fractions(e,mu,sig,q):
 v=primitive(np.r_[6,e,12],mu,sig,q);return np.diff(v[1:-1])/(v[-1]-v[0])
def basis(e,d):
 t=2*(e-e[0])/(e[-1]-e[0])-1;out=[]
 for j in range(d+1):
  z=np.zeros(j+1);z[-1]=1;out.append(np.diff(C.chebval(t,C.chebint(z)))/2)
 return np.array(out).T

def setup(deg,ei):
 nt=deg+1;B=np.zeros((338,nt+5));u=(edges[:NS+1]-6)/(edges[NS]-6);B[:NS,:nt]=np.array([np.diff(betainc(j+1,deg-j+1,u))/(deg+1) for j in range(deg+1)]).T;B[NS:,nt:]=basis(edges[NS:],4);n=np.array(I['eras'][ei]['counts'])[:338];qs=[next(q for q in YS if q['era']==eras[ei] and q['process']==p) for p in ['upsilon1s','upsilon2s']];rates=np.sum(I['eras'][ei]['templates'],axis=1);return B,n,qs,rates

def model(p,deg,ei,setupvals=None):
 B,n,qs,rates=setup(deg,ei) if setupvals is None else setupvals;nc=deg+6;cont=B@p[:nc];ys=[p[nc+j]*rates[j]*fractions(edges,p[nc+2+j],qs[j]['sigmaL'],qs[j]) for j in range(2)];return cont,ys

def exact_min(c):
 z=C.chebroots(C.chebder(c));rr=[v.real for v in z if abs(v.imag)<1e-7 and -1<v.real<1];return float(min(C.chebval(np.r_[-1,1,rr],c)))

def fit(args):
 deg,ei,data,seed=args;sv=setup(deg,ei);B,n,qs,rates=sv;n=n if data is None else np.array(data);nt=deg+1;nc=nt+5;u=np.linspace(0,1,257);grids=[np.array([comb(deg,j)*u**j*(1-u)**(deg-j) for j in range(deg+1)]).T,C.chebvander(2*u-1,4)];cr=next(f for f in CR if f['idx']==1 and f['era']==eras[ei]);cp=np.array(cr['p']);p0=np.r_[max(1,n[:NS].sum()-sum(rates[j]*np.sum(I['eras'][ei]['fractions'][j][:NS])*cp[5+j] for j in range(2))),np.zeros(deg),cp[:5],cp[5:7],cp[[7,9]]]
 # Seed SR from polynomial least squares to the observed remainder, then preserve positive support.
 yy=sum(p0[nc+j]*rates[j]*fractions(edges,p0[nc+2+j],qs[j]['sigmaL'],qs[j]) for j in range(2));co=np.linalg.lstsq(B[:NS,:nt],n[:NS]-yy[:NS],rcond=None)[0]
 p0[:nt]=np.maximum(co,1.0)
 def fun(*pars):
  p=np.array(pars);v1=grids[0]@p[:nt];v2=grids[1]@p[nt:nc]
  if min(v1)<=0 or min(v2)<=0:return 1e25+1e12*(np.minimum(v1,0)@np.minimum(v1,0)+np.minimum(v2,0)@np.minimum(v2,0))
  cont,ys=model(p,deg,ei,sv);mu=cont+ys[0]+ys[1]
  if min(mu)<=0:return 1e25
  return 2*np.sum(mu-n+xlogy(n,n/mu))
 alt=p0.copy();alt[-1]=qs[1]['mean'];alt[nc+1]=.5;starts=[p0,alt] if seed is None else [np.array(seed),p0,alt];best=None;attempts=[]
 for start in starts:
  m=Minuit(fun,*start);m.errordef=1
  for j in range(nt):m.limits[j]=(1e-5,None)
  m.limits[nt]=(1e-5,None)
  for j in range(2):m.limits[nc+j]=(0,10);m.limits[nc+2+j]=(qs[j]['mean']-.2,qs[j]['mean']+.2)
  m.errors=tuple([max(1,n[:NS].sum()*.005)]*nt+[max(1,n[NS:].sum()*.01)]*5+[.05,.1,.005,.005]);m.strategy=2;m.tol=.0001;m.migrad(ncall=12000);m.hesse();attempts.append(dict(D=m.fval,valid=bool(m.valid)))
  if best is None or m.fval<best.fval:best=m
 m=best;p=np.array(m.values);cov=np.array(m.covariance) if m.covariance is not None else np.zeros((len(p),len(p)));cont,ys=model(p,deg,ei,sv);mu=cont+ys[0]+ys[1];mins=[float(min(p[:nt])),exact_min(p[nt:nc])];dev=2*(mu-n+xlogy(n,n/mu));bound=[bool(min(p[j]-m.limits[j][0],m.limits[j][1]-p[j])<.001*(m.limits[j][1]-m.limits[j][0])) for j in range(nc,nc+4)]
 return dict(degree=deg,era=eras[ei],p=p.tolist(),cov=cov.tolist(),errors=list(m.errors),D=m.fval,D_SR=float(dev[:NS].sum()),D_CR=float(dev[NS:].sum()),valid=bool(m.valid and m.fmin.has_posdef_covar and min(mins)>0),accurate=bool(m.fmin.has_accurate_covar),minimum=mins,boundary=bound,mu=mu.tolist(),continuum=cont.tolist(),Y=[y.tolist() for y in ys],counts=n.tolist(),attempts=attempts,edm=m.fmin.edm)
if __name__=='__main__':
 original=json.loads((P/'runs/20261030_sr_degree_selection/observed.json').read_text());out=[]
 for deg in range(3,7):
  jobs=[]
  for ei,e in enumerate(eras):
   f=next(f for f in original if f['degree']==deg and f['era']==e);poly=C.Chebyshev(f['p'][:deg+1]).convert(kind=Polynomial)(Polynomial([-1,2]));coeff=[sum(poly.coef[k]*comb(j,k)/comb(deg,k) for k in range(j+1)) for j in range(deg+1)];seed=np.r_[np.maximum(coeff,1.),f['p'][deg+1:]];jobs.append((deg,ei,None,seed))
  with ProcessPoolExecutor(max_workers=3) as pool:
   for f in pool.map(fit,jobs):out.append(f);(R/'observed.json').write_text(json.dumps(out));print(deg,f['era'],f['valid'],f['D'],flush=True)
 print('DONE',flush=True)
