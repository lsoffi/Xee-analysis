from run import *
from scipy.special import betainc
from scipy.stats import chi2
lo=6+232/75;hi=6+338/75;fits=[];mask=(edges[:-1]>=lo-1e-8)&(edges[1:]<=hi+1e-8);ix=np.flatnonzero(mask);u=(edges[np.r_[ix,ix[-1]+1]]-lo)/(hi-lo);B=np.array([np.diff(betainc(j+1,6-j,u))/6 for j in range(6)]).T
for ei in range(4):
 n=np.array(I['eras'][ei]['counts'])[mask];Ytemp=np.array(I['eras'][ei]['templates'])[:,mask].T;A=np.column_stack([B,Ytemp])
 def fun(*p):
  mu=A@p
  return 2*np.sum(mu-n+xlogy(n,n/mu)) if min(mu)>0 else 1e30
 def grad(*p):return 2*A.T@(1-n/(A@p))
 m=Minuit(fun,*([max(1,n.sum()-Ytemp.sum())]*6+[1.,1.]),grad=grad);m.errordef=1
 for j in range(6):m.limits[j]=(1e-6,None)
 for j in [6,7]:m.limits[j]=(0,10)
 m.strategy=2;m.tol=1e-5;m.migrad();m.hesse();assert m.valid
 fits.append(dict(p=list(m.values),mu=(A@np.array(m.values)).tolist(),D=m.fval))
g=dict(lo=lo,hi=hi,degree=5,fits=fits,D=sum(f['D'] for f in fits),family='positive Bernstein5 fixed original Upsilon shapes');g['p_asym']=float(chi2.sf(g['D'],4*(106-8)));(R/'generator5.json').write_text(json.dumps(g));print(g['D'],g['p_asym'])
