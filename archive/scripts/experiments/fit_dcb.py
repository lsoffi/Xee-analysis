import os
os.environ['MPLCONFIGDIR']='work/mplconfig'
import json, re
from pathlib import Path
import numpy as np
from scipy.special import erf
from scipy.stats import chi2
from iminuit import Minuit
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
D=np.genfromtxt('work/data_obs_resonant.csv',delimiter=',',names=True)
meta=Path('work/data_metadata.txt').read_text()
obs=[n for n in D.dtype.names if n not in ('entry','weight','weight_error_lo','weight_error_hi')]
assert len(obs)==1,obs
name=obs[0]; x=D[name]; w=D['weight']
match=re.search(r'observable '+name+r' min (\S+) max (\S+) bins (\d+) unit (.*)',meta)
low,high=float(match[1]),float(match[2]); unit=match[4].strip()
assert np.all(w>=0), 'negative weights require a different statistical interpretation'
isbinned='class RooDataHist' in meta
N=w.sum()
def shape(t,al,nl,ar,nr):
 out=np.exp(-.5*t*t); l=t < -al; r=t > ar
 out[l]=np.exp(-.5*al*al-nl*np.log1p((-al-t[l])*al/nl))
 out[r]=np.exp(-.5*ar*ar-nr*np.log1p((t[r]-ar)*ar/nr))
 return out
def integ(a,b,al,nl,ar,nr):
 a,b=np.broadcast_arrays(np.asarray(a,float),np.asarray(b,float)); out=np.zeros_like(a)
 def tail(u,v,alpha,n):
  # Integral in positive distance from transition.
  p=1-n; A=1+u*alpha/n; B=1+v*alpha/n
  if abs(p)<1e-6: val=np.log(B/A)
  else: val=A**p*np.expm1(p*np.log(B/A))/p
  return np.exp(-alpha*alpha/2)*n/alpha*val
 l=a < -al; z=np.minimum(b,-al);out[l]+=tail(-al-z[l],-al-a[l],al,nl)
 r=b > ar; z=np.maximum(a,ar);out[r]+=tail(z[r]-ar,b[r]-ar,ar,nr)
 aa=np.maximum(a,-al);bb=np.minimum(b,ar);c=bb>aa
 out[c]+=np.sqrt(np.pi/2)*(erf(bb[c]/np.sqrt(2))-erf(aa[c]/np.sqrt(2)))
 return out
if isbinned:
 centers=np.sort(np.unique(x)); spacing=np.diff(centers)
 assert np.allclose(spacing,spacing[0]),'variable bins need explicit bin edges'
 bw=spacing[0]; edges=np.r_[centers-bw/2,centers[-1]+bw/2]; hist=np.array([w[x==c].sum() for c in centers])
else:
 edges=np.linspace(low,high,61);hist,_=np.histogram(x,bins=edges,weights=w)
def probs(edges,mu,sigma,alphaL,nL,alphaR,nR):
 z=(np.asarray(edges)-mu)/sigma
 return integ(z[:-1],z[1:],alphaL,nL,alphaR,nR)/integ((low-mu)/sigma,(high-mu)/sigma,alphaL,nL,alphaR,nR)
def nll(mu,sigma,alphaL,nL,alphaR,nR):
 norm=sigma*float(integ((low-mu)/sigma,(high-mu)/sigma,alphaL,nL,alphaR,nR))
 if norm<=0:return 1e100
 if isbinned:
  p=probs(edges,mu,sigma,alphaL,nL,alphaR,nR);return float(-np.sum(hist*np.log(np.maximum(p,1e-300))))
 t=(x-mu)/sigma;s=shape(t,alphaL,nL,alphaR,nR)
 return float(-np.sum(w*(np.log(np.maximum(s,1e-300))-np.log(norm))))
qs=np.interp([.16,.5,.84],np.cumsum(w[np.argsort(x)])/N,np.sort(x))
width=max((qs[2]-qs[0])/2,(high-low)/1000)
fits=[]
for alpha in [.7,1.5,2.5]:
 for n in [2.,8.]:
  m=Minuit(nll,mu=qs[1],sigma=width*.7,alphaL=alpha,nL=n,alphaR=alpha,nR=n)
  m.errordef=Minuit.LIKELIHOOD;m.limits['mu']=(low,high);m.limits['sigma']=((high-low)*1e-5,(high-low))
  for k in ['alphaL','alphaR']:m.limits[k]=(.05,15)
  for k in ['nL','nR']:m.limits[k]=(.05,200)
  m.migrad(ncall=15000);fits.append(m)
m=min(fits,key=lambda a:a.fval);m.hesse()
pars=list(m.parameters);values=np.array(list(m.values));cov=np.array(m.covariance) if m.covariance is not None else np.full((6,6),np.nan)
weighted=not np.allclose(w,1)
if weighted and not isbinned:
 def lp(p):
  mu,sig,al,nl,ar,nr=p
  return np.log(np.maximum(shape((x-mu)/sig,al,nl,ar,nr),1e-300))-np.log(sig*float(integ((low-mu)/sig,(high-mu)/sig,al,nl,ar,nr)))
 scores=[]
 for i,v in enumerate(values):
  h=max(abs(v)*1e-5,1e-6);a=values.copy();b=values.copy();a[i]+=h;b[i]-=h;scores.append((lp(a)-lp(b))/(2*h))
 S=np.array(scores);J=(S*w)@(S*w).T;cov=cov@J@cov
err=np.sqrt(np.maximum(np.diag(cov),0));expected=N*probs(edges,*values)
mask=expected>=5;pearson=float(np.sum((hist[mask]-expected[mask])**2/expected[mask]));ndof=int(mask.sum()-7)
fitinfo={'dataset':'data_obs_resonant','observable':name,'unit':unit,'fit_range':[low,high],'entries':len(x),'sum_weights':float(N),'weighted':weighted,'binned_input':isbinned,'fit_method':'binned multinomial likelihood' if isbinned else 'unbinned weighted maximum likelihood' if weighted else 'unbinned maximum likelihood','parameters':{p:{'value':float(v),'error':float(e)} for p,v,e in zip(pars,values,err)},'covariance':cov.tolist(),'valid':bool(m.valid),'accurate_covariance':bool(m.fmin.has_accurate_covar),'parameters_at_limit':bool(m.fmin.has_parameters_at_limit),'edm':float(m.fmin.edm),'nll':float(m.fval),'pearson_chi2_bins_expected_ge5':pearson,'ndof_approx':ndof,'pvalue_approx':float(chi2.sf(pearson,ndof)) if ndof>0 and not weighted else None,'multistart_nll':[float(f.fval) for f in fits]}
Path('outputs/dcb_fit_results.json').write_text(json.dumps(fitinfo,indent=2))
np.savetxt('outputs/dcb_covariance.csv',cov,delimiter=',',header=','.join(pars))
fig,(ax,pax)=plt.subplots(2,1,figsize=(8,7),sharex=True,gridspec_kw={'height_ratios':[3,1]},layout='constrained')
cent=(edges[1:]+edges[:-1])/2
v,_=np.histogram(x,bins=edges,weights=w*w) if not isbinned else (hist,None)
ax.errorbar(cent,hist,yerr=np.sqrt(v),fmt='o',ms=3,color='black',label='data_obs_resonant')
ax.stairs(expected,edges,color='#c43b35',lw=1.8,label='Double Crystal Ball (bin integrals)')
ax.set_ylabel('Events / bin' if not weighted else 'Sum of weights / bin');ax.legend(loc='upper left');ax.set_title('Double Crystal Ball fit — data_obs_resonant')
summary='\n'.join(f'{p} = {v:.5g} ± {e:.2g}' for p,v,e in zip(pars,values,err));ax.text(.98,.96,summary,transform=ax.transAxes,va='top',ha='right',fontsize=9,bbox={'facecolor':'white','alpha':.9,'edgecolor':'.8'})
pull=(hist-expected)/np.sqrt(np.maximum(expected,1e-12));pax.axhline(0,color='.4');pax.axhspan(-2,2,color='.9');pax.plot(cent,pull,'o',ms=3,color='black');pax.set_ylabel('Pearson pull');pax.set_xlabel(name+(f' [{unit}]' if unit else ''));pax.set_ylim(-max(4,np.max(abs(pull))*1.1),max(4,np.max(abs(pull))*1.1))
fig.savefig('outputs/dcb_fit.png',dpi=180);fig.savefig('outputs/dcb_fit.pdf')
print(json.dumps(fitinfo,indent=2));print(m.fmin)
