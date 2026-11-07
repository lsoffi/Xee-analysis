from pathlib import Path
import sys
sys.path.insert(0,str(Path.cwd()/'scripts/sr_degree_selection'))
from fit import *
R=P/'runs/20261030_sr_bernstein7_bias'
SIG=json.loads((P/'runs/20261029_upsilon_peak_windows/parameters.json').read_text())['signal'];RATES=json.loads((R/'signal_rates.json').read_text());GEN=json.loads((R/'generators.json').read_text());ANCH=json.loads((R/'asimov_anchors.json').read_text());LIM=json.loads((P/'runs/20261030_sr7_cr3_scan/summary.json').read_text())

TANCH=[json.loads(l) for l in (R/'toy_anchors.jsonl').read_text().splitlines()]
def poisson_deviance(mu,n):
 mask=n>0;z=(mu[mask]-n[mask])/n[mask];term=z-np.log1p(z);small=abs(z)<1e-4;q=z[small];term[small]=q*q*(.5+q*(-1/3+q*(.25-q/5)));return 2*(np.sum(n[mask]*term)+np.sum(mu[~mask]))
def runone(arg):
 gen,deg,mass,inject,toyindex,verify=arg;nc=deg+6;npar=nc+4;nt=deg+1;svs=[setup(deg,ei) for ei in range(4)];anchors=[next(f for f in TANCH if f['toy']==toyindex and f['degree']==deg and f['era']==e) for e in eras];generators=[next(f for f in GEN if f['degree']==gen and f['era']==e) for e in eras];rates=next(f['rates'] for f in RATES if abs(f['mass']-mass)<1e-7);signals=[]
 for ei,e in enumerate(eras):
  q=next(q for q in SIG if q['era']==e and abs(q['mass']-mass)<1e-7);signals.append(rates[ei]*fractions(edges,q['mean'],q['sigmaL'],q))
 rtrue=next(f['limits']['expected_50'] for f in LIM if abs(f['mass']-mass)<1e-7) if inject else 0.;data=[np.array(g['mu'])+rtrue*s for g,s in zip(generators,signals)];rng=np.random.default_rng(765421+toyindex+(int(round(mass*1000)) if inject else 0));data=[rng.poisson(x) for x in data];grid1=C.chebvander(np.linspace(-1,1,257),deg);grid2=C.chebvander(np.linspace(-1,1,257),4)
 def fun(*args):
  p=np.array(args);v=0.
  for ei in range(4):
   a=p[ei*npar:(ei+1)*npar];c1=grid1@a[:nt];c2=grid2@a[nt:nc]
   if min(c1)<=0 or min(c2)<=0:return 1e25+1e12*(np.minimum(c1,0)@np.minimum(c1,0)+np.minimum(c2,0)@np.minimum(c2,0))
   cont,ys=model(a,deg,ei,svs[ei]);mu=cont+ys[0]+ys[1]+p[-1]*signals[ei]
   if min(mu)<=0:return 1e25+1e12*np.minimum(mu,0)@np.minimum(mu,0)
   v+=poisson_deviance(mu,data[ei])
  return v
 def grad(*args):
  p=np.array(args);v=np.zeros(len(p))
  for ei in range(4):
   off=ei*npar;a=p[off:off+npar];B,n,qs,ratesy=svs[ei];c1=grid1@a[:nt];c2=grid2@a[nt:nc]
   if min(c1)<=0 or min(c2)<=0:
    v[off:off+nt]=2e12*grid1.T@np.minimum(c1,0);v[off+nt:off+nc]=2e12*grid2.T@np.minimum(c2,0);return v
   shapes=[ratesy[j]*fractions(edges,a[nc+2+j],qs[j]['sigmaL'],qs[j]) for j in range(2)];mu=B@a[:nc]+a[nc]*shapes[0]+a[nc+1]*shapes[1]+p[-1]*signals[ei];factor=2*(1-data[ei]/np.maximum(mu,1e-12));v[off:off+nc]=B.T@factor
   for j in range(2):
    v[off+nc+j]=shapes[j]@factor;h=1e-5;der=ratesy[j]*a[nc+j]*(fractions(edges,a[nc+2+j]+h,qs[j]['sigmaL'],qs[j])-fractions(edges,a[nc+2+j]-h,qs[j]['sigmaL'],qs[j]))/(2*h);v[off+nc+2+j]=der@factor
   v[-1]+=signals[ei]@factor
  return v
 start=np.r_[np.concatenate([f['p'] for f in anchors]),rtrue];m=Minuit(fun,*start,grad=grad);m.errordef=1;m.precision=1e-12
 for ei in range(4):
  off=ei*npar;B,n,qs,ratesy=svs[ei];m.limits[off]=(1e-5,None);m.limits[off+nt]=(1e-5,None)
  for j in range(2):m.limits[off+nc+j]=(0,10);m.limits[off+nc+2+j]=(qs[j]['mean']-.2,qs[j]['mean']+.2)
  for j in range(npar):m.errors[off+j]=max(1e-5,anchors[ei]['errors'][j])
 m.limits[-1]=(-max(5,10*rtrue),max(5,10*rtrue));m.errors[-1]=max(.01,rtrue*.2);m.strategy=0;m.tol=.001;m.migrad(ncall=10000);p=np.array(m.values);mins=[]
 # Exact observed likelihood Hessian: linear yields/polynomials plus the two nonlinear means.
 info=np.zeros((len(p),len(p)))
 for ei in range(4):
  off=ei*npar;a=p[off:off+npar];B,n,qs,ratesy=svs[ei];base=[];ders=[];seconds=[]
  for j in range(2):
   h=1e-4;f0=ratesy[j]*fractions(edges,a[nc+2+j],qs[j]['sigmaL'],qs[j]);fp=ratesy[j]*fractions(edges,a[nc+2+j]+h,qs[j]['sigmaL'],qs[j]);fm=ratesy[j]*fractions(edges,a[nc+2+j]-h,qs[j]['sigmaL'],qs[j]);base.append(f0);ders.append((fp-fm)/(2*h));seconds.append((fp-2*f0+fm)/(h*h))
  mu=B@a[:nc]+a[nc]*base[0]+a[nc+1]*base[1]+p[-1]*signals[ei];jac=np.column_stack([B,base[0],base[1],a[nc]*ders[0],a[nc+1]*ders[1],signals[ei]]);local=jac.T@((data[ei]/mu**2)[:,None]*jac);factor=1-data[ei]/mu
  for j in range(2):
   cross=factor@ders[j];local[nc+j,nc+2+j]+=cross;local[nc+2+j,nc+j]+=cross;local[nc+2+j,nc+2+j]+=a[nc+j]*(factor@seconds[j])
  ids=list(range(off,off+npar))+[len(p)-1];info[np.ix_(ids,ids)]+=local
 scale=np.sqrt(np.maximum(abs(np.diag(info)),1e-20));normalized=info/np.outer(scale,scale);eig=np.linalg.eigvalsh(normalized);hessian_valid=bool(np.all(np.isfinite(normalized)) and min(eig)>1e-10);cov=np.linalg.inv(normalized)/np.outer(scale,scale) if hessian_valid else np.full(info.shape,np.nan);rerror=float(np.sqrt(cov[-1,-1])) if hessian_valid else float(m.errors[-1]);verification=None
 if verify:
  m.hesse();verification=dict(analytic_r_error=rerror,minuit_r_error=m.errors[-1],relative_difference=abs(rerror/m.errors[-1]-1),analytic_hessian_valid=hessian_valid,minuit_valid=bool(m.valid))

 for ei in range(4):a=p[ei*npar:(ei+1)*npar];mins.extend([exact_min(a[:nt]),exact_min(a[nt:nc])])
 valid=bool(not m.fmin.is_above_max_edm and not m.fmin.has_reached_call_limit and np.isfinite(m.fval) and min(mins)>0 and abs(p[-1])<4.999);return dict(generator=gen,degree=deg,mass=mass,injection='expected' if inject else 'zero',rtrue=rtrue,rhat=p[-1],error=rerror,bias_sigma=(p[-1]-rtrue)/rerror if rerror>0 else None,toy=toyindex,error_valid=bool(hessian_valid),minuit_covariance_valid=bool(m.fmin.has_posdef_covar),above_max_edm=bool(m.fmin.is_above_max_edm),call_limit=bool(m.fmin.has_reached_call_limit),hessian_verification=verification,minimum_scaled_hessian_eigenvalue=float(min(eig)),data_hash=__import__('hashlib').sha256(np.array(data).tobytes()).hexdigest(),valid=valid,D=m.fval,edm=m.fmin.edm,minimum=mins,parameters=p.tolist(),errors=list(m.errors),conditional_signal_shape=True)


if __name__=='__main__':
 import time
 initial=[json.loads(l) for l in (R/'stable_toys.jsonl').read_text().splitlines()];out=[]
 for d,mass in [(3,6.5),(5,7.5),(6,6.5),(6,7.5),(6,8.5),(7,7.5)]:
  start=time.perf_counter();q=runone(('bern7',d,mass,False,0,False));q['seconds']=time.perf_counter()-start;ref=next(x for x in initial if x['degree']==d and x['mass']==mass and x['toy']==0);q['delta_D']=q['D']-ref['D'];q['delta_r_sigma']=(q['rhat']-ref['rhat'])/ref['error'];out.append(q);print(d,mass,q['valid'],q['seconds'],q['delta_D'],q['delta_r_sigma'],flush=True)
 (R/'strategy0_check.json').write_text(json.dumps(out))
