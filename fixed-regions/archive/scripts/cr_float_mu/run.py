from pathlib import Path
import sys,json,numpy as np
from scipy.special import erf,xlogy
from numpy.polynomial import chebyshev as C
from iminuit import Minuit
from concurrent.futures import ProcessPoolExecutor
P=Path.cwd();R=P/'runs/20261030_cr_float_mu';old=P/'runs/20261030_cr_optimization';I=json.loads((old/'input.json').read_text());S=json.loads((old/'scan.json').read_text());Y=json.loads((P/'runs/20261029_upsilon_peak_windows/parameters.json').read_text())['upsilon'];edges=np.array(I['edges']);eras=[a['era'] for a in I['eras']]
def primitive(x,mu,sig,q):
 t=(np.asarray(x)-mu)/sig;v=np.sqrt(np.pi/2)*erf(t/np.sqrt(2))
 for side,sgn in [('L',-1),('R',1)]:
  a=abs(q['alpha'+side]);n=q['n'+side];z=sgn*t;sel=z>a;B=n/a-a;A=(n/a)**n*np.exp(-a*a/2)
  tail=A*np.log((B+z[sel])/(B+a)) if abs(n-1)<1e-8 else A*((B+z[sel])**(1-n)-(B+a)**(1-n))/(1-n)
  v[sel]=sgn*(np.sqrt(np.pi/2)*erf(a/np.sqrt(2))+tail)
 return v*sig

def fractions(e,mu,sig,q):
 v=primitive(np.r_[6,e,12],mu,sig,q);return np.diff(v[1:-1])/(v[-1]-v[0])

def fit(args):
 idx,ei,data,seed_override=args;cfg=S[idx];deg=cfg['degree'];nt=deg+1;lo,hi=cfg['lo'],cfg['hi'];mask=(edges[:-1]>=lo-1e-8)&(edges[1:]<=hi+1e-8);ix=np.flatnonzero(mask);e=edges[np.r_[ix,ix[-1]+1]];t=2*(e-lo)/(hi-lo)-1;basis=[]
 for j in range(nt):
  c=np.zeros(j+1);c[-1]=1;basis.append(np.diff(C.chebval(t,C.chebint(c)))/2)
 B=np.array(basis).T;grid=C.chebvander(np.linspace(-1,1,129),deg);n=np.array(I['eras'][ei]['counts'])[mask] if data is None else np.array(data);N=n.sum();qs=[next(q for q in Y if q['era']==eras[ei] and q['process']==pr) for pr in ['upsilon1s','upsilon2s']];rates=np.sum(I['eras'][ei]['templates'],axis=1)
 def model(p,ret=False):
  cont=B@p[:nt];ys=[p[nt+j]*rates[j]*fractions(e,p[nt+2+2*j],p[nt+3+2*j],qs[j]) for j in range(2)]
  return (cont,ys) if ret else cont+ys[0]+ys[1]
 def fun(*args):
  p=np.array(args);pol=grid@p[:nt]
  if min(pol)<=0:return 1e25+1e12*np.sum(np.minimum(pol,0)**2)
  mu=model(p)
  if min(mu)<=0:return 1e25
  return 2*np.sum(mu-n+xlogy(n,n/mu))
 init=cfg['fits'][ei]['p']+[v for q in qs for v in [q['mean'],q['sigmaL']]];alt=init.copy();alt[nt]=max(.3,init[nt]);alt[nt+1]=.5;alt[-3]*=.85;alt[-1]*=1.15
 starts=[init,alt] if seed_override is None else [seed_override,init];best=None;attempts=[]
 for st in starts:
  m=Minuit(fun,*st);m.errordef=1;m.limits[0]=(1e-5,None)
  for j in range(2):m.limits[nt+j]=(0,10);m.limits[nt+2+2*j]=(qs[j]['mean']-.2,qs[j]['mean']+.2);m.limits[nt+3+2*j]=(.02,.4)
  
  for j in range(2):m.values[nt+3+2*j]=qs[j]['sigmaL'];m.fixed[nt+3+2*j]=True
  m.errors=tuple([max(1,N*.01)]*nt+[.05,.1,.005,.005,.005,.005]);m.strategy=1;m.tol=.0001;m.migrad(ncall=5000);m.hesse();attempts.append(dict(D=m.fval,valid=bool(m.valid)))
  if best is None or m.fval<best.fval:best=m
 m=best;p=np.array(m.values);cov=np.zeros((len(p),len(p)));free=[j for j in range(len(p)) if not m.fixed[j]]
 if m.covariance is not None:
  for ja in free:
   for jb in free:cov[ja,jb]=m.covariance[ja,jb];roots=C.chebroots(C.chebder(p[:nt]));rr=[z.real for z in roots if abs(z.imag)<1e-7 and -1<z.real<1];minpol=float(min(C.chebval(np.r_[-1,1,rr],p[:nt])));cont,ys=model(p,True)
 def predictions(p):
  out=[]
  for a,b in [(6,9.093333333333334),(10.506666666666668,12)]:out.append(sum(p[nt+j]*rates[j]*fractions(np.array([a,b]),p[nt+2+2*j],p[nt+3+2*j],qs[j])[0] for j in range(2)))
  return np.array(out)
 jac=[]
 for j in range(len(p)):
  h=1e-5*max(1,abs(p[j]));pp=p.copy();pm=p.copy();pp[j]+=h;pm[j]-=h;jac.append((predictions(pp)-predictions(pm))/(2*h))
 jac=np.array(jac).T;pcov=jac@cov@jac.T;bounds=[]
 for j in range(nt,len(p)):
  a,b=m.limits[j];bounds.append(bool(not m.fixed[j] and min(p[j]-a,b-p[j])<1e-3*(b-a)))
 return dict(idx=idx,era=eras[ei],lo=lo,hi=hi,degree=deg,p=p.tolist(),errors=list(m.errors),cov=cov.tolist(),D=m.fval,valid=bool(m.valid and m.fmin.has_posdef_covar and minpol>0),accurate=bool(m.fmin.has_accurate_covar),boundary=bounds,attempts=attempts,minpoly=minpol,counts=n.tolist(),x=((e[:-1]+e[1:])/2).tolist(),continuum=cont.tolist(),Y=[y.tolist() for y in ys],mu=model(p).tolist(),prediction=predictions(p).tolist(),prediction_error=np.sqrt(np.maximum(0,np.diag(pcov))).tolist(),nominal_shapes=[[q['mean'],q['sigmaL']] for q in qs])
if __name__=='__main__':
 # Check analytic DCB bin integrals against the ROOT export for all eight nominal shapes.
 maxdiff=0
 for ei,era in enumerate(eras):
  for j,pr in enumerate(['upsilon1s','upsilon2s']):
   q=next(q for q in Y if q['era']==era and q['process']==pr);v=fractions(edges,q['mean'],q['sigmaL'],q);maxdiff=max(maxdiff,float(np.max(abs(v-np.array(I['eras'][ei]['fractions'][j])))))
 assert maxdiff<1e-8; (R/'validation.json').write_text(json.dumps(dict(max_root_bin_fraction_difference=maxdiff,mu_bounds='nominal +/-0.2 GeV independently for each resonance and era',sigma='fixed at original MC values',tails='original fixed',likelihood='bin-integrated Poisson',seeds=2)))
 res=[]
 with ProcessPoolExecutor(max_workers=3) as pool:
  for r in pool.map(fit,[(idx,ei,None,None) for idx in range(len(S)) for ei in range(4)]):
   res.append(r);(R/'fits.json').write_text(json.dumps(res));print(len(res),'/192',r['idx'],r['era'],r['valid'],r['boundary'],flush=True)
 print('DONE',flush=True)
