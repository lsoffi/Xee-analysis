"""Read native C0 SR+CR fit snapshots; never minimize a likelihood."""
from pathlib import Path
import json, math, os
import ROOT
import numpy as np

ERAS=['2022','2022EE','2023','2023BPix']
REGIONS=[('SR','mass_SR'),('CR','mass_CR')]


def dump(path,data):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    def clean(x):
        if isinstance(x,dict):return {str(k):clean(v) for k,v in x.items()}
        if isinstance(x,(list,tuple)):return [clean(v) for v in x]
        if isinstance(x,(float,np.floating)):return float(x) if math.isfinite(x) else None
        if isinstance(x,np.integer):return int(x)
        return x
    tmp=path.with_suffix(path.suffix+'.tmp.'+str(os.getpid()));tmp.write_text(json.dumps(clean(data),indent=2,allow_nan=False));tmp.replace(path)


def model_side(w):
    return 'low' if w.var('mass_SR').getMax()<=w.var('mass_CR').getMin()+1e-8 else 'high'


def apply_fit(w,fit):
    for a in list(fit.floatParsFinal())+list(fit.constPars()):
        if a.InheritsFrom('RooAbsCategory'):
            if w.cat(a.GetName()):w.cat(a.GetName()).setIndex(a.getIndex())
        elif w.var(a.GetName()):w.var(a.GetName()).setVal(a.getVal())


def polynomial_minimum(w):
    out={};degree=1
    while w.var(f'local_c{degree+1}_{ERAS[0]}'):degree+=1
    for era in ERAS:
        for reg,deg,suffix in [('SR',degree,''),('CR',4,'_CR')]:
            p=np.polynomial.Chebyshev([1.]+[w.var(f'local_c{k}_{era}{suffix}').getVal() for k in range(1,deg+1)])
            roots=p.deriv().roots();points=[-1.,1.]+[float(q.real) for q in roots if abs(q.imag)<1e-9 and -1<q.real<1]
            out[f'{reg}_{era}']=float(min(p(q) for q in points))
    return out


def add_pdf(w,stage,channel):
    mc=w.obj('ModelConfig' if stage=='s' else 'ModelConfig_bonly');cp=mc.GetPdf().getPdf(channel)
    if cp.InheritsFrom('RooAddPdf'):return cp
    return next(p for p in cp.getComponents() if p.InheritsFrom('RooAddPdf') and p.canBeExtended())


def inspect_model(w,stage,data_name='data_obs'):
    mc=w.obj('ModelConfig' if stage=='s' else 'ModelConfig_bonly');data=w.data(data_name)
    bins={};deviance=0.;minimum=math.inf
    for region,xname in REGIONS:
        x=w.var(xname);obs=ROOT.RooArgSet(x);width=(x.getMax()-x.getMin())/x.getBins()
        for i,era in enumerate(ERAS,1):
            channel=f'{region}_ch{i}';pdf=add_pdf(w,stage,channel);expected=pdf.expectedEvents(obs)
            sub=data.reduce(ROOT.RooFit.Cut(f'CMS_channel==CMS_channel::{channel}'))
            curves={p.GetName():[] for p in pdf.pdfList()};yields=[a.getVal() for a in pdf.coefList()]
            centers=[];counts=[];total=[];ns=0.;term=0.;valid=True
            for j in range(sub.numEntries()):
                row=sub.get(j);xx=row.getRealValue(xname);n=sub.weight();x.setVal(xx)
                mu=pdf.getVal(obs)*expected*width;minimum=min(minimum,mu)
                centers.append(xx);counts.append(n);total.append(mu);ns+=n
                if mu<=0 or not math.isfinite(mu):valid=False
                elif n:term+=n*math.log(n/mu)
                for p,y in zip(pdf.pdfList(),yields):curves[p.GetName()].append(p.getVal(obs)*y*width)
            d=2*(expected-ns+term) if valid else None
            deviance=deviance+d if d is not None else math.inf
            bins[channel]=dict(region=region,era=era,mass=centers,data=counts,total=total,components=curves,
                expected_events=expected,observed_events=ns,bin_width=width,deviance=d,
                component_yields={p.GetName():y for p,y in zip(pdf.pdfList(),yields)})
    aux=0.;globs=list(mc.GetGlobalObservables());assert len(globs)==10
    for g in globs:
        name=g.GetName();assert name.endswith('_In') and w.var(name[:-3]),name
        aux+=(w.var(name[:-3]).getVal()-g.getVal())**2
    return dict(exact_deviance=deviance+aux,data_deviance=deviance,aux_penalty=aux,auxiliary_constraints=10,minimum_bin=minimum),bins


def continuity(w,stage):
    side=model_side(w);edge=w.var('mass_SR').getMax() if side=='low' else w.var('mass_SR').getMin();out={}
    for i,era in enumerate(ERAS,1):
        regions={}
        for region,xname in REGIONS:
            x=w.var(xname);old=x.getVal();x.setVal(edge);obs=ROOT.RooArgSet(x);pdf=add_pdf(w,stage,f'{region}_ch{i}')
            components={}
            for p,y in zip(pdf.pdfList(),pdf.coefList()):
                name=p.GetName();proc='dy' if name.startswith('local_continuum_') else next(z for z in ['Zd','upsilon1s','upsilon2s'] if name.startswith('shape_'+z+'_'))
                components[proc]=float(p.getVal(obs)*y.getVal())
            regions[region]=dict(total=float(pdf.getVal(obs)*pdf.expectedEvents(obs)),components=components)
            x.setVal(old)
        values={k:regions['SR']['components'].get(k,0.)-regions['CR']['components'].get(k,0.) for k in set(regions['SR']['components'])|set(regions['CR']['components'])}
        values['total']=regions['SR']['total']-regions['CR']['total']
        relative={k:abs(v)/max(abs(regions['SR']['total'] if k=='total' else regions['SR']['components'].get(k,0.)),abs(regions['CR']['total'] if k=='total' else regions['CR']['components'].get(k,0.)),1e-12) for k,v in values.items()}
        out[era]=dict(**regions,difference_SR_minus_CR=values,relative_difference=relative,
            SR_continuum_yield=float(w.function('local_yield_'+era).getVal()),CR_continuum_yield=float(w.var('local_yield_'+era+'_CR').getVal()))
    return dict(side=side,boundary=edge,units='events/GeV',eras=out,maximum_relative_difference=max(v for e in out.values() for v in e['relative_difference'].values()))


def read_stage(ws_path,fd_path,stage,output=None):
    source=ROOT.TFile.Open(str(ws_path));w=source.Get('w') if source and not source.IsZombie() else None
    if not w:return dict(readable=False,source=str(fd_path),failure='missing_workspace'),None
    f=ROOT.TFile.Open(str(fd_path))
    if not f or f.IsZombie() or not f.Get('fit_'+stage):
        if f:f.Close()
        source.Close();return dict(readable=False,source=str(fd_path),failure='missing_fit_'+stage),None
    try:
        fit=f.Get('fit_'+stage)
        assert not any(p.InheritsFrom('RooMultiPdf') for p in w.allPdfs())
        pre=[dict(name=v.GetName(),value=float(w.var(v.GetName()).getVal()),error=float(w.var(v.GetName()).getError())) for v in fit.floatParsFinal() if w.var(v.GetName())]
        premap={v['name']:v['value'] for v in pre};apply_fit(w,fit)
        if stage=='b':w.var('r').setVal(0)
        params=[dict(name=v.GetName(),value=v.getVal(),error=v.getError(),low=v.getMin(),high=v.getMax(),prefit=premap.get(v.GetName())) for v in fit.floatParsFinal()]
        metric,bins=inspect_model(w,stage);mins=polynomial_minimum(w);r=fit.floatParsFinal().find('r')
        rv=r.getVal() if r else 0.;error=r.getError() if r else None
        boundary=bool(r and min(abs(rv-r.getMin()),abs(rv-r.getMax()))<1e-4*max(1.,r.getMax()-r.getMin()))
        means={v.GetName():v.getVal() for v in w.allVars() if v.GetName().startswith('Ymean_')}
        assert len(means)==8 and all(not w.var(n).isConstant() and fit.floatParsFinal().find(n) for n in means),'All eight Upsilon means must float'
        a=dict(readable=True,source=str(fd_path),status=fit.status(),covQual=fit.covQual(),edm=fit.edm(),parameters=params,prefit_parameters=pre,
            prefit_values_are_initializations_not_constraints=True,r=rv,error=error,r_boundary=boundary,continuum_minima=mins,minimum=min(mins.values()),
            floating_means=means,mean_boundary=[n for n in means if min(abs(w.var(n).getVal()-w.var(n).getMin()),abs(w.var(n).getVal()-w.var(n).getMax()))<1e-4*(w.var(n).getMax()-w.var(n).getMin())],**metric)
        a['edm_valid']=math.isfinite(a['edm']) and 0<=a['edm']<.01
        a['point_valid']=a['status']==0 and a['edm_valid'] and a['minimum']>0 and a['minimum_bin']>0 and math.isfinite(a['exact_deviance']) and not boundary
        a['numerical_ok']=a['point_valid'] and a['covQual']==3
        a['continuity']=continuity(w,stage);a['continuity_passed']=a['continuity']['maximum_relative_difference']<1e-8
        if output:
            output=Path(output);state=dict(variables={v.GetName():[float(v.getVal()),bool(v.isConstant())] for v in w.allVars()},categories={v.GetName():v.getIndex() for v in w.allCats()})
            previous=output.with_suffix('.state.json');same=output.exists() and previous.exists() and json.loads(previous.read_text())==state
            if not same:w.saveSnapshot('postfit_'+stage,w.allVars());w.writeToFile(str(output));dump(previous,state)
        return a,bins
    finally:
        f.Close();source.Close()
