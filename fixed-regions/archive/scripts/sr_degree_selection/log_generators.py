from fit import *
from numpy.polynomial.legendre import leggauss
from scipy.stats import chi2
G=json.loads((R/'generators.json').read_text())
def one(arg):
 d,ei=arg;B,n,qs,rates=setup(8,ei);original=next(g for g in G if g['degree']==8 and g['era']==eras[ei]);grid,w=leggauss(12);qpoints=(edges[:NS,None]+edges[1:NS+1,None])/2+(edges[1:NS+1,None]-edges[:NS,None])*grid/2;qt=2*(qpoints-6)/(edges[NS]-6)-1;Q=C.chebvander(qt,d)[:,:,1:];coef=np.linalg.lstsq(C.chebvander(2*((edges[:NS]+edges[1:NS+1])/2-6)/(edges[NS]-6)-1,d),np.log(np.maximum(np.array(original['continuum'])[:NS],1e-6)),rcond=None)[0][1:];nt=d+1;nc=nt+5;init=np.r_[np.sum(original['continuum'][:NS]),coef,original['p'][9:14],original['p'][14:]];BC=basis(edges[NS:],4);gridC=C.chebvander(np.linspace(-1,1,257),4)
 def model2(p):
  z=np.einsum('ijk,k->ij',Q,p[1:nt]);z-=np.max(z);yy=np.exp(z)@w;yy/=yy.sum();cont=np.r_[p[0]*yy,BC@p[nt:nc]];ys=[p[nc+j]*rates[j]*fractions(edges,p[nc+2+j],qs[j]['sigmaL'],qs[j]) for j in range(2)];return cont,ys
 def fun(*args):
  p=np.array(args);v=gridC@p[nt:nc]
  if min(v)<=0:return 1e25+1e12*np.minimum(v,0)@np.minimum(v,0)
  cont,ys=model2(p);mu=cont+ys[0]+ys[1]
  if min(mu)<=0:return 1e25
  return 2*np.sum(mu-n+xlogy(n,n/mu))
 m=Minuit(fun,*init);m.errordef=1;m.limits[0]=(1,None);m.limits[nt]=(1e-5,None)
 for j in range(1,nt):m.limits[j]=(-10,10)
 for j in range(2):m.limits[nc+j]=(0,10);m.limits[nc+2+j]=(qs[j]['mean']-.2,qs[j]['mean']+.2)
 m.strategy=2;m.tol=.0001;m.migrad(ncall=15000);m.hesse();p=np.array(m.values);cont,ys=model2(p);mu=cont+ys[0]+ys[1];return dict(degree=f'exp{d}',era=eras[ei],D=m.fval,valid=bool(m.valid and m.fmin.has_posdef_covar and exact_min(p[nt:nc])>0),p=p.tolist(),mu=mu.tolist(),continuum=cont.tolist(),Y=[y.tolist() for y in ys],family='exponential of Chebyshev in SR, Chebyshev4 CR; shared floating Y mu/norm',p_asym_era=float(chi2.sf(m.fval,338-len(p))))
if __name__=='__main__':
 out=[]
 with ProcessPoolExecutor(max_workers=3) as pool:
  for f in pool.map(one,[(d,ei) for d in [5,6] for ei in range(4)]):out.append(f);print(f['degree'],f['era'],f['D'],f['valid'],flush=True)
 (R/'log_generators.json').write_text(json.dumps(out));(R/'generators_extended.json').write_text(json.dumps(G+out))
 anchors=json.loads((R/'asimov_anchors.json').read_text())
 for gen in [7,8,'exp5','exp6']:
  jobs=[]
  for ei,e in enumerate(eras):
   ob=next(f for f in G if f['degree']==8 and f['era']==e);truth=next(f for f in G+out if f['degree']==gen and f['era']==e);jobs.append((8,ei,truth['mu'],ob['p']))
  with ProcessPoolExecutor(max_workers=3) as pool:
   for f in pool.map(fit,jobs):f['generator_degree']=gen;anchors.append(f);print('anchor8',gen,f['era'],f['valid'],flush=True)
 (R/'asimov_anchors_extended.json').write_text(json.dumps(anchors))
