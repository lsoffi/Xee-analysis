from fit import *
SIG=json.loads((P/'runs/20261029_upsilon_peak_windows/parameters.json').read_text())['signal'];RATES=json.loads((R/'signal_rates.json').read_text());GEN=json.loads((R/'generators.json').read_text());ANCH=json.loads((R/'asimov_anchors.json').read_text());LIM=json.loads((P/'runs/20261030_sr7_cr3_scan/summary.json').read_text())

def runone(arg):
 gen,deg,mass,inject=arg;nc=deg+6;npar=nc+4;nt=deg+1;svs=[setup(deg,ei) for ei in range(4)];anchors=[next(f for f in ANCH if f['generator_degree']==gen and f['degree']==deg and f['era']==e) for e in eras];generators=[next(f for f in GEN if f['degree']==gen and f['era']==e) for e in eras];rates=next(f['rates'] for f in RATES if abs(f['mass']-mass)<1e-7);signals=[]
 for ei,e in enumerate(eras):
  q=next(q for q in SIG if q['era']==e and abs(q['mass']-mass)<1e-7);signals.append(rates[ei]*fractions(edges,q['mean'],q['sigmaL'],q))
 rtrue=next(f['limits']['expected_50'] for f in LIM if abs(f['mass']-mass)<1e-7) if inject else 0.;data=[np.array(g['mu'])+rtrue*s for g,s in zip(generators,signals)];grid1=C.chebvander(np.linspace(-1,1,257),deg);grid2=C.chebvander(np.linspace(-1,1,257),4)
 def fun(*args):
  p=np.array(args);v=0.
  for ei in range(4):
   a=p[ei*npar:(ei+1)*npar];c1=grid1@a[:nt];c2=grid2@a[nt:nc]
   if min(c1)<=0 or min(c2)<=0:return 1e25+1e12*(np.minimum(c1,0)@np.minimum(c1,0)+np.minimum(c2,0)@np.minimum(c2,0))
   cont,ys=model(a,deg,ei,svs[ei]);mu=cont+ys[0]+ys[1]+p[-1]*signals[ei]
   if min(mu)<=0:return 1e25+1e12*np.minimum(mu,0)@np.minimum(mu,0)
   v+=2*np.sum(mu-data[ei]+xlogy(data[ei],data[ei]/mu))
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
 start=np.r_[np.concatenate([f['p'] for f in anchors]),rtrue];start[-1]+=.01;analytic=grad(*start);numeric=[]
 for j in range(len(start)):
  h=1e-6*max(1,abs(start[j]));pp=start.copy();pm=start.copy();pp[j]+=h;pm[j]-=h;numeric.append((fun(*pp)-fun(*pm))/(2*h))
 numeric=np.array(numeric);err=np.max(abs(analytic-numeric)/np.maximum(1,abs(numeric)));assert err<.001,(err,analytic,numeric)
 return dict(max_scaled_gradient_difference=float(err),parameters=len(start),signal_fraction_SR_CR=[float(sum(s))/rates[i] for i,s in enumerate(signals)])
 m=Minuit(fun,*start,grad=grad);m.errordef=1
 for ei in range(4):
  off=ei*npar;B,n,qs,ratesy=svs[ei];m.limits[off]=(1e-5,None);m.limits[off+nt]=(1e-5,None)
  for j in range(2):m.limits[off+nc+j]=(0,10);m.limits[off+nc+2+j]=(qs[j]['mean']-.2,qs[j]['mean']+.2)
  for j in range(npar):m.errors[off+j]=max(1e-5,anchors[ei]['errors'][j])
 m.limits[-1]=(-max(5,10*rtrue),max(5,10*rtrue));m.errors[-1]=max(.01,rtrue*.2);m.strategy=1;m.tol=.0001;m.migrad(ncall=20000);m.hesse();p=np.array(m.values);mins=[]
 for ei in range(4):a=p[ei*npar:(ei+1)*npar];mins.extend([exact_min(a[:nt]),exact_min(a[nt:nc])])
 valid=bool(m.valid and m.fmin.has_posdef_covar and min(mins)>0);return dict(generator=gen,degree=deg,mass=mass,injection='expected' if inject else 'zero',rtrue=rtrue,rhat=p[-1],error=m.errors[-1],bias_sigma=(p[-1]-rtrue)/m.errors[-1] if m.errors[-1]>0 else None,valid=valid,D=m.fval,edm=m.fmin.edm,minimum=mins,parameters=p.tolist(),errors=list(m.errors),conditional_signal_shape=True)

if __name__=='__main__':
 q=runone((8,6,7.5,False));(R/'gradient_check.json').write_text(json.dumps(q,indent=2));print(q)
