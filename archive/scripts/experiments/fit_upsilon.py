exec(open('work/fit_dcb.py').read().split('qs=np.interp')[0])
from iminuit.cost import ExtendedBinnedNLL
original_edges=edges.copy()
# Rebin complete adjacent native bins; preserve sum(w) and sum(w^2).
y=w.reshape(-1,2).sum(axis=1)
v=(D['weight_error_lo']**2).reshape(-1,2).sum(axis=1)
e=original_edges[::2]
sel=(e[:-1]>=8.8-1e-8)&(e[1:]<=9.8+0.014)
idx=np.flatnonzero(sel);edges=e[idx[0]:idx[-1]+2];hist=y[sel];variance=v[sel]
low,high=float(edges[0]),float(edges[-1]);N=hist.sum()
def cdf(z,mu,sigma,alphaL,nL,alphaR,nR,yield_):
 t=(np.asarray(z)-mu)/sigma
 return yield_*integ((low-mu)/sigma,t,alphaL,nL,alphaR,nR)/integ((low-mu)/sigma,(high-mu)/sigma,alphaL,nL,alphaR,nR)
cost=ExtendedBinnedNLL(np.c_[hist,variance],edges,cdf)
fits=[]
for sigma in [.10,.18,.25]:
 for alpha in [.7,1.5,2.5]:
  m=Minuit(cost,mu=9.43,sigma=sigma,alphaL=alpha,nL=4,alphaR=alpha,nR=4,yield_=N)
  m.limits['mu']=(9.1,9.65);m.limits['sigma']=(.01,.5)
  for p in ['alphaL','alphaR']:m.limits[p]=(.1,10)
  for p in ['nL','nR']:m.limits[p]=(1.01,100)
  m.limits['yield_']=(.05*N,5*N);m.tol=.01;m.migrad(ncall=30000);fits.append(m)
m=min(fits,key=lambda a:a.fval);m.strategy=2;m.migrad(ncall=50000);m.hesse()
vals=np.array(list(m.values));errs=np.array(list(m.errors));names=list(m.parameters)
try:m.minos('mu','sigma',ncall=30000)
except Exception as exc:print('MINOS warning',exc)
expected=np.diff(cdf(edges,*vals));good=variance>0
pull=np.full(len(hist),np.nan);pull[good]=(hist[good]-expected[good])/np.sqrt(variance[good])
res={'dataset':'data_obs_resonant','signal_only':True,'model':'double-sided Crystal Ball, no background','rebin':2,'native_bin_width':float(np.diff(original_edges)[0]),'rebinned_bin_width':float(np.diff(edges)[0]),'nominal_range':[8.8,9.8],'actual_range_aligned_to_bins':[low,high],'bins':len(hist),'empty_bins':int(np.count_nonzero(hist==0)),'sum_weights_in_window':float(N),'sum_weights_squared_in_window':float(variance.sum()),'effective_entries':float(N*N/variance.sum()),'method':'iminuit ExtendedBinnedNLL weighted Bohm-Zech approximation using sum(w) and sum(w^2), exact DCB bin integrals','parameters':{p:{'value':float(a),'error_hesse':float(b),'limits':list(m.limits[p])} for p,a,b in zip(names,vals,errs)},'minos':{k:{'lower':float(a.lower),'upper':float(a.upper),'valid':bool(a.is_valid)} for k,a in m.merrors.items()},'covariance':np.array(m.covariance).tolist(),'valid_minimum':bool(m.valid),'accurate_covariance':bool(m.fmin.has_accurate_covar),'forced_positive_definite_covariance':bool(m.fmin.has_made_posdef_covar),'parameters_at_limit':bool(m.fmin.has_parameters_at_limit),'edm':float(m.fmin.edm),'weighted_deviance':float(m.fval),'ndof_nominal':int(len(hist)-m.nfit),'asymptotic_pvalue_nominal':float(chi2.sf(m.fval,len(hist)-m.nfit)),'chi2_sumw2_diagnostic':float(np.nansum(pull*pull)),'multistart_deviance':[float(f.fval) for f in fits]}
Path('outputs/upsilon_dcb_rebin2_results.json').write_text(json.dumps(res,indent=2))
np.savetxt('outputs/upsilon_dcb_rebin2_bins.csv',np.c_[edges[:-1],edges[1:],hist,variance,expected,pull],delimiter=',',header='low,high,sumw,sumw2,expected,pull_sumw2')
fig,(ax,pa)=plt.subplots(2,1,figsize=(8.5,7),sharex=True,gridspec_kw={'height_ratios':[3,1]},layout='constrained')
cent=(edges[:-1]+edges[1:])/2
ax.errorbar(cent,hist,yerr=np.sqrt(variance),fmt='o',ms=3.5,color='.2',lw=.8,label='data_obs_resonant, rebin 2 (SumW2)')
ax.stairs(expected,edges,color='#b72b3d',lw=2,label='Double Crystal Ball; weighted likelihood')
ax.set_title('Υ(1S) region — signal only');ax.set_ylabel('Sum of weights / bin');ax.legend(loc='upper left',fontsize=9)
ss=f'μ = {vals[0]:.4f} ({m.merrors["mu"].lower:+.4f}, {m.merrors["mu"].upper:+.4f})\nσ = {vals[1]:.4f} ({m.merrors["sigma"].lower:+.4f}, {m.merrors["sigma"].upper:+.4f})\nMINOS errors; tails at limits\nWeighted deviance / nominal ndf\n= {m.fval:.2f} / {len(hist)-m.nfit}'
ax.text(.98,.96,ss,transform=ax.transAxes,ha='right',va='top',fontsize=9,bbox={'facecolor':'white','edgecolor':'.8','alpha':.9})
pa.axhspan(-2,2,color='.92');pa.axhline(0,color='.4',lw=.8);pa.plot(cent,pull,'o',ms=3.5,color='.2');pa.set_ylabel('(data−fit)/σ(data)');pa.set_ylim(-max(3.5,np.nanmax(abs(pull))*1.2),max(3.5,np.nanmax(abs(pull))*1.2));pa.set_xlabel('mass (unit not specified in file)')
fig.savefig('outputs/upsilon_dcb_rebin2.png',dpi=180);fig.savefig('outputs/upsilon_dcb_rebin2.pdf')
print(json.dumps(res,indent=2));print(m.fmin)
