exec(open('work/fit_dcb.py').read().split('qs=np.interp')[0])
from scipy.integrate import quad
# Keep native histogram uncertainties. Bins with exactly zero variance do not
# constrain a Gaussian chi-square; their number is explicitly reported.
variance=D['weight_error_lo']**2; use=variance>0
assert np.allclose(D['weight_error_lo'],D['weight_error_hi'])
# Validate the analytic tail integrals independently.
for pp in [(1.3,2.,1.8,3.),(.5,.8,.7,1.2)]:
 numerical=quad(lambda q:shape(np.array([q]),*pp)[0],-4,5,points=[-pp[0],pp[2]],epsabs=1e-9)[0]
 assert np.allclose(integ(-4,5,*pp),numerical,rtol=1e-7)
def cost(mu,sigma,alphaL,nL,alphaR,nR,yield_):
 pred=yield_*probs(edges,mu,sigma,alphaL,nL,alphaR,nR)
 return float(np.sum((hist[use]-pred[use])**2/variance[use]))
fits=[]
for sig in [.08,.18,.35]:
 for al in [.5,1.5,3.]:
  m=Minuit(cost,mu=9.43,sigma=sig,alphaL=al,nL=3,alphaR=al,nR=5,yield_=N)
  m.errordef=Minuit.LEAST_SQUARES;m.limits['mu']=(6,12);m.limits['sigma']=(.005,2)
  for k in ['alphaL','alphaR']:m.limits[k]=(.1,10)
  for k in ['nL','nR']:m.limits[k]=(1.01,100)
  m.limits['yield_']=(1,10*N);m.strategy=1;m.tol=.01;m.migrad(ncall=20000);fits.append(m)
m=min(fits,key=lambda v:v.fval);m.strategy=2;m.migrad(ncall=30000);m.hesse()
P=list(m.parameters);val=np.array(list(m.values));errors=np.array(list(m.errors));cov=np.array(m.covariance)
expected=m.values['yield_']*probs(edges,*val[:6]);ndof=int(use.sum()-len(P));pvalue=float(chi2.sf(m.fval,ndof))
res={'dataset':'data_obs_resonant','class':'RooDataHist','observable':name,'unit':unit or 'not specified in input','fit_range':[low,high],'native_bins':len(hist),'bins_in_chi2':int(use.sum()),'zero_variance_bins_excluded':int((~use).sum()),'sum_weights':float(N),'sum_weights_squared':float(variance.sum()),'method':'chi-square with native SumW2 errors and bin-integrated double-sided Crystal Ball; floating normalization','parameters':{k:{'value':float(v),'error_hesse':float(e)} for k,v,e in zip(P,val,errors)},'covariance':cov.tolist(),'chi2':float(m.fval),'ndof':ndof,'pvalue_asymptotic':pvalue,'valid_minimum':bool(m.valid),'accurate_covariance':bool(m.fmin.has_accurate_covar),'forced_positive_definite_covariance':bool(m.fmin.has_made_posdef_covar),'parameters_at_limit':bool(m.fmin.has_parameters_at_limit),'edm':float(m.fmin.edm),'multistart_chi2':[float(f.fval) for f in fits]}
# Profile errors: robust checks for weakly constrained shape parameters.
try:
 m.minos('mu','sigma',ncall=20000)
 res['minos']={k:{'lower':float(e.lower),'upper':float(e.upper),'valid':bool(e.is_valid)} for k,e in m.merrors.items()}
except Exception as e:res['minos_error']=str(e)
Path('outputs/dcb_fit_results.json').write_text(json.dumps(res,indent=2));np.savetxt('outputs/dcb_covariance.csv',cov,delimiter=',',header=','.join(P))
fig,(ax,pa)=plt.subplots(2,1,figsize=(9,7),sharex=True,gridspec_kw={'height_ratios':[3,1]},layout='constrained')
cent=(edges[1:]+edges[:-1])/2
ax.errorbar(cent[use],hist[use],yerr=np.sqrt(variance[use]),fmt='.',ms=4,lw=.7,color='.25',label='Data (native SumW2 errors)')
ax.stairs(expected,edges,color='#be2e3b',lw=1.7,label='Double Crystal Ball, bin integrals')
ax.set_title('data_obs_resonant — weighted binned fit');ax.set_ylabel('Sum of weights / bin');ax.legend(loc='upper left',fontsize=9)
ss='\n'.join(f'{p} = {v:.5g} ± {e:.2g}' for p,v,e in zip(P[:6],val[:6],errors[:6]));ss+=f'\nχ²/ndf = {m.fval:.1f}/{ndof}'
ax.text(.98,.96,ss,ha='right',va='top',transform=ax.transAxes,fontsize=9,bbox={'facecolor':'white','alpha':.9,'edgecolor':'.8'})
pull=(hist[use]-expected[use])/np.sqrt(variance[use]);pa.axhspan(-2,2,color='.92');pa.axhline(0,color='.4',lw=.8);pa.plot(cent[use],pull,'.',color='.25',ms=4);pa.set_ylabel('(data−fit)/error');pa.set_xlabel('mass (unit not specified in file)');pa.set_ylim(-max(4,np.max(abs(pull))*1.15),max(4,np.max(abs(pull))*1.15))
fig.savefig('outputs/dcb_fit.png',dpi=180);fig.savefig('outputs/dcb_fit.pdf')
np.savetxt('outputs/dcb_binned_fit.csv',np.c_[cent,hist,np.sqrt(variance),expected,use.astype(int)],delimiter=',',header='mass,sum_weights,sumw2_error,dcb_expected,used_in_chi2')
print(json.dumps(res,indent=2));print(m.fmin)
