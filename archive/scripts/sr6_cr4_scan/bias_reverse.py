from pathlib import Path
import sys
sys.path.insert(0,str(Path.cwd()/'scripts/sr_bernstein_selection'))
from fit import *
R=P/'runs/20261030_sr_cheb7_bern_bias'
DEST=P/'runs/20261030_sr6_cr4_scan'
INPUT=json.loads((DEST/'signals.json').read_text());SIG=[dict(q,mass=z['mass']) for z in INPUT for q in z['eras']];RATES=[dict(mass=z['mass'],rates=[q['full_rate_r1'] for q in z['eras']]) for z in INPUT]
GEN=json.loads((R/'generators.json').read_text());ANCH=json.loads((R/'asimov_anchors.json').read_text())

def poisson_deviance(mu,n):
 z=(mu-n)/np.maximum(n,1e-300);t=np.where(abs(z)<1e-4,z*z*(.5+z*(-1/3+z*(.25-z/5))),z-np.log1p(z));return 2*np.sum(np.where(n>0,n*t,mu))

def runone(arg):
 gen,deg,mass,inject=arg;nc=deg+6;npar=nc+4;nt=deg+1;svs=[setup(deg,ei) for ei in range(4)];anchors=[next(f for f in ANCH if f['generator_degree']==gen and f['degree']==deg and f['era']==e) for e in eras];generators=[next(f for f in GEN if f['degree']==gen and f['era']==e) for e in eras];rates=next(f['rates'] for f in RATES if abs(f['mass']-mass)<1e-7);signals=[]
 for ei,e in enumerate(eras):
  q=next(q for q in SIG if q['era']==e and abs(q['mass']-mass)<1e-7);signals.append(rates[ei]*fractions(edges,q['mean'],q['sigmaL'],q))
 rtrue=1. if inject else 0.;data=[np.array(g['mu'])+rtrue*s for g,s in zip(generators,signals)];u=np.linspace(0,1,257);grid1=np.array([comb(deg,j)*u**j*(1-u)**(deg-j) for j in range(deg+1)]).T;grid2=C.chebvander(np.linspace(-1,1,257),4)
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
  off=ei*npar;B,n,qs,ratesy=svs[ei];
  for j in range(nt):m.limits[off+j]=(1e-5,None)
  m.limits[off+nt]=(1e-5,None)
  for j in range(2):m.limits[off+nc+j]=(0,10);m.limits[off+nc+2+j]=(qs[j]['mean']-.2,qs[j]['mean']+.2)
  for j in range(npar):m.errors[off+j]=max(1e-5,anchors[ei]['errors'][j])
 m.limits[-1]=(-max(5,10*rtrue),max(5,10*rtrue));m.errors[-1]=max(.01,rtrue*.2);m.strategy=2;m.tol=.001;m.migrad(ncall=30000);m.hesse();p=np.array(m.values);mins=[]
 for ei in range(4):a=p[ei*npar:(ei+1)*npar];mins.extend([float(min(a[:nt])),exact_min(a[nt:nc])])
 valid=bool(m.valid and m.fmin.has_posdef_covar and min(mins)>0);correlations={}
 if m.covariance is not None:
  for ei,e in enumerate(eras):
   for j,name in enumerate(['normY1','normY2','muY1','muY2']):
    idx=ei*npar+nc+j;den=np.sqrt(max(0,m.covariance[-1,-1]*m.covariance[idx,idx]));correlations[e+'_'+name]=float(m.covariance[-1,idx]/den) if den>0 else None
 return dict(generator=gen,degree=deg,mass=mass,injection='one' if inject else 'zero',rtrue=rtrue,rhat=p[-1],error=m.errors[-1],bias_sigma=(p[-1]-rtrue)/m.errors[-1] if m.errors[-1]>0 else None,valid=valid,D=m.fval,edm=m.fmin.edm,minimum=mins,parameters=p.tolist(),errors=list(m.errors),conditional_signal_shape=True,r_bounds=list(m.limits[-1]),r_at_boundary=bool(min(p[-1]-m.limits[-1][0],m.limits[-1][1]-p[-1])<.001*(m.limits[-1][1]-m.limits[-1][0])),r_Y_correlations=correlations,signal_CR_fraction=[float(sum(x[NS:])/rates[i]) for i,x in enumerate(signals)])
if __name__=='__main__':
 jobs=[(gen,deg,k/10,inj) for gen in ['cheb7'] for deg in [6] for k in range(60,91) for inj in [False,True]]
 with ProcessPoolExecutor(max_workers=2) as pool,(DEST/'reverse.jsonl').open('w') as out:
  for j,r in enumerate(pool.map(runone,jobs)):
   out.write(json.dumps(r)+'\n');out.flush()
   if j%14==13:print(j+1,'/',len(jobs),'latest',r['generator'],r['degree'],r['mass'],r['valid'],r['bias_sigma'],flush=True)
 print('DONE',flush=True)
