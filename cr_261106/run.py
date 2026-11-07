"""Bounded comparison of simple CR models, exact binned Poisson likelihood."""
import argparse
import csv
import hashlib
import json
import math
import time
from array import array
from pathlib import Path
import numpy as np
import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.gInterpreter.Declare((Path(__file__).parent/'kernel.h').read_text())
def objective(cr,era,degree,float_width,tail_count=0):
    rr=[r for r in DATA_ROWS if r['era']==era and float(r['bin_low_GeV'])>=cr['low']-1e-9 and float(r['bin_high_GeV'])<=cr['high']+1e-9]
    o=ROOT.ThreeCR.Objective()
    for r in rr:o.edges.push_back(float(r['bin_low_GeV']));o.counts.push_back(float(r['events']))
    o.edges.push_back(float(rr[-1]['bin_high_GeV']))
    assert abs(o.edges[0]-cr['low'])<1e-9 and abs(o.edges.back()-cr['high'])<1e-9
    o.degree=degree;o.floatWidth=float_width;o.tailCount=tail_count
    pp=[]
    for i,(_,m) in enumerate(cr['resonances']):
        p=shape(era,m);pp.append(p)
        s=ROOT.ThreeCR.Shape()
        for k,v in [('mean',p['mean']),('sigma',p['sigmaL']),('al',p['alphaL']),('ar',p['alphaR']),('nl',p['nL']),('nr',p['nR'])]:setattr(s,k,v)
        o.shape[i]=s
    o.prepare()
    return o,pp


def fit(o,cr,start=None,variant=0,hesse=True,fast=False):
    total=sum(o.counts)
    n=o.npar()
    initial=[math.log(max(1,total*f)) for f in cr['fractions']]+[0.,0.]
    if o.floatWidth:initial += [0.,0.]
    initial += [0.]*o.tailCount
    initial += [0.]*(n-len(initial))
    if o.degree==-2:
        initial[o.nbase():]=[1.25,math.log(.3),math.log(.6)]
    if start is not None:initial=list(start)
    if variant:
        initial[3]+=(.25 if variant==1 else -.25)
        initial[4]+=(-.25 if variant==1 else .25)
        if o.floatWidth:
            initial[5]+=math.log(.8 if variant==1 else 1.2)
            initial[6]+=math.log(1.2 if variant==1 else .8)
        for i in range(o.nbase(),n):initial[i]+=(.15 if variant==1 else -.15)
    names=['ln_N_continuum','ln_N_'+cr['resonances'][0][0],'ln_N_'+cr['resonances'][1][0],
           'mean_shift_sigma_'+cr['resonances'][0][0],'mean_shift_sigma_'+cr['resonances'][1][0]]
    limits=[(math.log(1e-5),math.log(3*total))]*3+[(-2.,2.)]*2
    if cr.get('ordered_means'):
        midpoint=(o.shape[0].mean+o.shape[1].mean)/2
        limits[3]=((o.edges.front()-o.shape[0].mean)/o.shape[0].sigma,(midpoint-o.shape[0].mean)/o.shape[0].sigma)
        limits[4]=((midpoint-o.shape[1].mean)/o.shape[1].sigma,(o.edges.back()-o.shape[1].mean)/o.shape[1].sigma)
    if o.floatWidth:
        names+=['ln_sigma_scale_'+cr['resonances'][0][0],'ln_sigma_scale_'+cr['resonances'][1][0]]
        limits += [(-math.log(4),math.log(4))]*2
    if o.tailCount:
        tailnames=['ln_alphaL_scale','ln_alphaR_scale','ln_nL_minus1_scale','ln_nR_minus1_scale']
        if o.fractionalN:
            names += [n+'_'+res[0] for res in cr['resonances'] for n in ['ln_alphaL_scale','ln_alphaR_scale','ln_nL_scale','ln_nR_scale']]
            af,nf=cr['alpha_fraction'],cr['n_fraction']
            limits += ([(math.log(1-af),math.log(1+af))]*2+[(math.log(1-nf),math.log(1+nf))]*2)*2
        elif o.independentAlpha:
            names += [n+'_'+res[0] for res in cr['resonances'] for n in tailnames[:2]]
            fraction=cr['alpha_fraction']
            limits += [(math.log(1-fraction),math.log(1+fraction))]*4
        elif o.tailCount==8:
            names += [n+'_'+res[0] for res in cr['resonances'] for n in tailnames]
            if cr.get('wide_tails'):
                for k in range(2):
                    s=o.shape[k]
                    limits += [(math.log(.05/s.al),math.log(5/s.al)),(math.log(.05/s.ar),math.log(5/s.ar)),
                               (math.log(.01/(s.nl-1)),math.log(99/(s.nl-1))),
                               (math.log(.01/(s.nr-1)),math.log(99/(s.nr-1)))]
                    if variant:
                        ti=(7 if o.floatWidth else 5)+4*k
                        initial[ti+2]=math.log(((2. if variant==1 else 40.)-1)/(s.nl-1))
                        initial[ti+3]=math.log(((2. if variant==1 else 40.)-1)/(s.nr-1))
            else:
                limits += ([(-math.log(3),math.log(3))]*2+[(-math.log(5),math.log(5))]*2)*2
        else:
            names += tailnames[:o.tailCount]
            limits += [(-math.log(3),math.log(3))]*min(2,o.tailCount)
            if o.tailCount==4:limits += [(-math.log(5),math.log(5))]*2
    if o.degree==-4:
        names+=['quadratic_slope_left','quadratic_slope_right'];limits +=[(0,1000)]*2
    elif o.degree==-3:
        names+=['polynomial_linear','polynomial_quadratic'];limits +=[(-100,100)]*2
    elif o.degree==-2:
        names+=['alpha_scaled','ln_beta_scaled','ln_gamma_scaled']
        limits += [(1.,5.),(math.log(.005),math.log(10.)),(math.log(.01),math.log(100.))]
    elif o.degree<0:
        names+=['exponential_slope'];limits +=[(-12,-1e-6) if cr.get('decreasing_exponential') else (-12,12)]
    else:
        names += [f'ln_b{i}_over_b0' for i in range(1,o.degree+1)]
        limits += [(-8,8)]*o.degree
    m=ROOT.Math.Factory.CreateMinimizer('Minuit2','Migrad')
    m.SetMaxFunctionCalls(50000 if not fast else 10000)
    m.SetMaxIterations(10000)
    m.SetTolerance(1e-5 if not fast else 1e-3)
    m.SetPrecision(1e-12)
    m.SetStrategy(2 if not fast else 1)
    m.SetPrintLevel(-1)
    callback=lambda p:o.evaluate(p)
    fn=ROOT.Math.Functor(callback,n)
    m.SetFunction(fn)
    for i,(name,value,(lo,hi)) in enumerate(zip(names,initial,limits)):
        value=min(hi-1e-5,max(lo+1e-5,value))
        m.SetLimitedVariable(i,name,value,.01,lo,hi)
    ok=bool(m.Minimize())
    migrad_status=int(m.Status())
    if hesse:m.Hesse()
    values=[float(m.X()[i]) for i in range(n)]
    errors=[float(m.Errors()[i]) for i in range(n)]
    covariance=[[float(m.CovMatrix(i,j)) for j in range(n)] for i in range(n)] if hesse else None
    components=[list(c) for c in o.components(array('d',values))]
    expected=[sum(c[i] for c in components) for i in range(len(o.counts))]
    D=2*o.evaluate(array('d',values))
    ndf=len(o.counts)-n
    flags=[names[i] for i,(lo,hi) in enumerate(limits) if min(values[i]-lo,hi-values[i])<.001*(hi-lo)]
    valid=ok and migrad_status==0 and math.isfinite(D) and m.Edm()<1e-3
    return dict(ok=ok,migrad_status=migrad_status,status=int(m.Status()),covQual=int(m.CovMatrixStatus()),
                edm=float(m.Edm()),calls=int(m.NCalls()),D=D,ndf=ndf,
                p_chi2_approx=float(ROOT.Math.chisquared_cdf_c(D,ndf)),
                numerical_valid=valid,covariance_valid=int(m.CovMatrixStatus())==3,
                parameter_names=names,values=values,errors=errors,covariance=covariance,boundary_flags=flags,
                expected=expected,components=components,
                poisson_deviance_residuals=[math.copysign(math.sqrt(max(0,2*(mu-obs+(obs*math.log(obs/mu) if obs else 0)))),obs-mu) for obs,mu in zip(o.counts,expected)])



import glob
configs=json.loads(Path('config.json').read_text())
for q in configs:
 era=q['era']; cr=q['cr']; group=cr['name']
 cr['low']=round(cr['low']*200)/200;cr['high']=round(cr['high']*200)/200
 f=ROOT.TFile.Open(glob.glob('/eos/user/n/npalmeri/shared/datasets_261106/261106/'+era+'/*.root')[0]);d=f.Get('w').data('data_obs')
 DATA_ROWS=[]
 for i in range(d.numEntries()):
  mass=d.get(i).getRealValue('mass')
  if cr['low']<mass<cr['high']:DATA_ROWS.append(dict(era=era,bin_low_GeV=round(mass-.0025,6),bin_high_GeV=round(mass+.0025,6),events=d.weight()))
 def shape(e,m):return q['pp'][[v[1] for v in cr['resonances']].index(m)]
 o,pp=objective(cr,era,q['degree'],q['float_width'],q['tail_count'])
 o.fractionalN=cr.get('fractional_n',False)
 r=fit(o,cr,q['seed'])
 record=dict(cr=cr,era=era,degree=q['degree'],float_width=q['float_width'],tail_count=q['tail_count'],signal_parameters=pp,fit=r,edges=list(o.edges),counts=list(o.counts))
 path=Path(group,era);path.mkdir(parents=True,exist_ok=True);(path/'display_fit.json').write_text(json.dumps(record))
 print('FIT',group,era,r['status'],r['covQual'],r['D'],r['ndf'],flush=True)
 f.Close()
