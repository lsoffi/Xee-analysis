"""ROOT extraction for the common-continuum SRL+CR+SRH scan."""
from pathlib import Path
import itertools, json, math, re, os
import ROOT
import numpy as np

ERAS = ['2022', '2022EE', '2023', '2023BPix']
REGIONS = [('SRL', 'mass_SRL'), ('CR', 'mass_CR'), ('SRH', 'mass_SRH')]

def dump(path, data):
    path = Path(path)
    def clean(x):
        if isinstance(x, dict): return {str(k):clean(v) for k,v in x.items()}
        if isinstance(x, (list,tuple)): return [clean(v) for v in x]
        if isinstance(x, (float,np.floating)): return float(x) if math.isfinite(x) else None
        if isinstance(x,np.integer): return int(x)
        return x
    temp=path.with_suffix(path.suffix+'.tmp.'+str(os.getpid()));temp.write_text(json.dumps(clean(data),indent=2,allow_nan=False));temp.replace(path)

def apply_fit(w, fit):
    for a in list(fit.floatParsFinal())+list(fit.constPars()):
        if a.InheritsFrom('RooAbsCategory'):
            if w.cat(a.GetName()): w.cat(a.GetName()).setIndex(a.getIndex())
        elif w.var(a.GetName()): w.var(a.GetName()).setVal(a.getVal())
    indices={}; sources={}
    for era in ERAS:
        name='pdfindex_'+era;cat=w.cat(name)
        if not cat: continue
        stored=fit.constPars().find(name)
        if stored and stored.InheritsFrom('RooAbsCategory'):
            idx=stored.getIndex();sources[era]='fit.constPars'
        else:
            c=fit.floatParsFinal().find('local_c1_'+era)
            b=fit.floatParsFinal().find('local_b1_'+era)
            idx=1 if b and b.getError()>0 and (not c or c.getError()==0) else 0
            sources[era]='active_parameter_errors'
        cat.setIndex(idx);indices[era]=idx
    return indices,sources

def polynomial_minimum(w):
    out={}
    for era in ERAS:
        for reg,degree,suffix in [('SR',7,''),('CR',4,'_CR')]:
            if reg=='SR' and w.cat('pdfindex_'+era).getIndex()==1:
                out[f'{reg}_{era}']=min([1.]+[w.var(f'local_b{k}_{era}').getVal() for k in range(1,degree+1)])
            else:
                coeff=[1.]+[w.var(f'local_c{k}_{era}{suffix}').getVal() for k in range(1,degree+1)]
                p=np.polynomial.Chebyshev(coeff);roots=p.deriv().roots()
                points=[-1.,1.]+[float(q.real) for q in roots if abs(q.imag)<1e-9 and -1<q.real<1]
                out[f'{reg}_{era}']=float(min(p(q) for q in points))
    return out

def inspect_model(w, stage):
    mc=w.obj('ModelConfig' if stage=='s' else 'ModelConfig_bonly')
    data=w.data('data_obs');bins={};deviance=0.;minimum=math.inf
    for region,xname in REGIONS:
        x=w.var(xname);obs=ROOT.RooArgSet(x);width=(x.getMax()-x.getMin())/x.getBins()
        for i,era in enumerate(ERAS,1):
            channel=f'{region}_ch{i}';cp=mc.GetPdf().getPdf(channel)
            adds=[p for p in cp.getComponents() if p.InheritsFrom('RooAddPdf')]
            pdf=cp if cp.InheritsFrom('RooAddPdf') else next(p for p in adds if p.canBeExtended())
            expected=pdf.expectedEvents(obs)
            sub=data.reduce(ROOT.RooFit.Cut(f'CMS_channel==CMS_channel::{channel}'))
            curves={p.GetName():[] for p in pdf.pdfList()};yields=[a.getVal() for a in pdf.coefList()]
            centers=[];counts=[];total=[];ns=0.;term=0.;valid=True
            for j in range(sub.numEntries()):
                row=sub.get(j);xx=row.getRealValue(xname);n=sub.weight();x.setVal(xx)
                mu=pdf.getVal(obs)*expected*width;minimum=min(minimum,mu)
                centers.append(xx);counts.append(n);total.append(mu);ns+=n
                if mu<=0 or not math.isfinite(mu): valid=False
                elif n:term+=n*math.log(n/mu)
                for p,yield_ in zip(pdf.pdfList(),yields):curves[p.GetName()].append(p.getVal(obs)*yield_*width)
            d=2*(expected-ns+term) if valid else None
            if d is None:deviance=math.inf
            else:deviance+=d
            bins[channel]=dict(region=region,era=era,mass=centers,data=counts,total=total,components=curves,
                               expected_events=expected,observed_events=ns,bin_width=width,deviance=d)
    aux=0.
    for glob in mc.GetGlobalObservables():
        name=glob.GetName()
        if not name.endswith('_In') or not w.var(name[:-3]):raise ValueError('Unknown constraint '+name)
        aux+=(w.var(name[:-3]).getVal()-glob.getVal())**2
    return dict(exact_deviance=deviance+aux,data_deviance=deviance,aux_penalty=aux,minimum_bin=minimum),bins

def read_stage(ws_path, fd_path, stage, output=None):
    source=ROOT.TFile.Open(str(ws_path));w=source.Get('w');f=ROOT.TFile.Open(str(fd_path))
    if not f or f.IsZombie() or not f.Get('fit_'+stage):
        source.Close();return dict(readable=False,source=str(fd_path),failure='missing_fit_'+stage),None
    fit=f.Get('fit_'+stage);indices,idx_source=apply_fit(w,fit)
    if stage=='b':w.var('r').setVal(0)
    params=[dict(name=v.GetName(),value=v.getVal(),error=v.getError(),low=v.getMin(),high=v.getMax()) for v in fit.floatParsFinal()]
    metric,bins=inspect_model(w,stage);mins=polynomial_minimum(w)
    r=fit.floatParsFinal().find('r')
    rv=r.getVal() if r else 0.;re_=r.getError() if r else None
    rboundary=bool(r and min(abs(rv-r.getMin()),abs(rv-r.getMax()))<1e-4*max(1,r.getMax()-r.getMin()))
    a=dict(readable=True,source=str(fd_path),status=fit.status(),covQual=fit.covQual(),edm=fit.edm(),
           parameters=params,indices=indices,indices_source=idx_source,r=rv,error=re_,r_boundary=rboundary,
           continuum_minima=mins,minimum=min(mins.values()),**metric)
    a['edm_valid']=math.isfinite(a['edm']) and 0<=a['edm']<.01
    a['point_valid']=a['status']==0 and a['edm_valid'] and a['minimum']>0 and a['minimum_bin']>0 and math.isfinite(a['exact_deviance']) and not rboundary
    a['numerical_ok']=a['point_valid'] and a['covQual']==3
    a['floating_means']={v.GetName():v.getVal() for v in w.allVars() if v.GetName().startswith('Ymean_')}
    assert len(a['floating_means'])==8
    assert all(not w.var(name).isConstant() and fit.floatParsFinal().find(name) for name in a['floating_means'])
    a['mean_boundary']=[name for name in a['floating_means'] if min(abs(w.var(name).getVal()-w.var(name).getMin()),abs(w.var(name).getVal()-w.var(name).getMax())) < 1e-4*(w.var(name).getMax()-w.var(name).getMin())]
    if output:
        # Generation templates are immutable once their physical state is fixed.
        # Repeated extraction must not change ROOT headers/UUIDs and file hashes.
        output=Path(output)
        def state(ws):
            return dict(variables={v.GetName():[float(v.getVal()),bool(v.isConstant())] for v in ws.allVars()},
                        categories={v.GetName():int(v.getIndex()) for v in ws.allCats()})
        desired=state(w);same=False
        if output.exists():
            old=ROOT.TFile.Open(str(output))
            if old and not old.IsZombie() and old.Get('w'):same=state(old.Get('w'))==desired
            if old:old.Close()
        if not same:
            w.saveSnapshot('postfit_'+stage,w.allVars());w.writeToFile(str(output))
        dump(output.with_suffix('.state.json'),desired)
    f.Close();source.Close();return a,bins

def read_tree(path):
    f=ROOT.TFile.Open(str(path))
    if not f or f.IsZombie():return []
    tree=f.Get('limit')
    out=[dict(limit=float(a.limit),quantile=float(a.quantileExpected)) for a in tree] if tree else []
    f.Close();return out
