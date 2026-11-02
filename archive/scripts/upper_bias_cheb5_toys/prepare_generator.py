"""Fit the Chebyshev5 SR / Chebyshev4 CR background truth on observed data."""
from pathlib import Path
import json, math, hashlib, subprocess, time, sys
import ROOT
import numpy as np

P=Path(__file__).resolve().parents[2]
R=P/'runs/20261030_upper_bias_cheb5_toys'
OLD=P/'runs/20261030_upper_envelope_scan'
ERAS=['2022','2022EE','2023','2023BPix']
ROOT.gROOT.SetBatch(True)
ROOT.gSystem.Load('libHiggsAnalysisCombinedLimit')

def sha(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def customize(w, fixed=False, parameters=None):
    """Keep the CR and all signal/Υ nodes, replace only the SR continuum."""
    keep=[];replacements=[];constraints=[]
    for era in ERAS:
        for j in [4,5]:
            v=ROOT.RooRealVar(f'local_c{j}_{era}','',0,-1,1)
            getattr(w,'import')(v);keep.append(v)
        coeff=ROOT.RooArgList(*[w.var(f'local_c{j}_{era}') for j in range(1,6)])
        pdf=ROOT.RooChebychev(f'generator_cheb5_{era}','',w.var('mass'),coeff)
        getattr(w,'import')(pdf,ROOT.RooFit.RecycleConflictNodes());keep.append(pdf)
        replacements.append(f'local_envelope_{era}={pdf.GetName()}')
        checks=[]
        for t in np.linspace(-1,1,513):
            terms=[np.polynomial.chebyshev.chebval(t,[0]*j+[1]) for j in range(1,6)]
            checks.append('(1'+''.join(f'+@{j}*({v:.17g})' for j,v in enumerate(terms))+'>1e-8)')
        constraint=ROOT.RooGenericPdf(f'positive_generator_cheb5_{era}','&&'.join(checks),coeff)
        getattr(w,'import')(constraint);keep.append(constraint);constraints.append(constraint.GetName())
    for name in ['ModelConfig','ModelConfig_bonly']:
        mc=w.obj(name)
        new=w.factory(f'EDIT::{name}_cheb5_generator('+mc.GetPdf().GetName()+','+','.join(replacements)+')')
        assert new
        for old in list(new._extraConstraints):
            if old.GetName().startswith('positive_continuum_') and not old.GetName().endswith('_CR'):
                new._extraConstraints.remove(old)
        for cn in constraints:new.addExtraConstraint(w.pdf(cn))
        mc.SetPdf(new)
        assert len(list(new._extraConstraints))==8
        names=[v.GetName() for v in new._extraConstraints]
        assert sum(n.startswith('positive_generator_cheb5_') for n in names)==4,names
    w.genobj('discreteParams').removeAll()
    for c in w.allCats():
        if c.GetName().startswith('pdfindex'):c.setConstant(True)
    if parameters:
        for a in parameters:
            if w.var(a['name']):w.var(a['name']).setVal(a['value'])
    if fixed:
        for v in w.allVars():
            if v.GetName().startswith(('local_','Ymean_','CMS_EXO25020_bkgScale_')):v.setConstant(True)
    w.var('r').setVal(0);w.var('mass').setVal(11.2);w.var('mass_CR').setVal(9.8)
    return keep

def export(w, fit=None):
    if fit:
        for v in fit.floatParsFinal():w.var(v.GetName()).setVal(v.getVal())
    w.var('r').setVal(0)
    parameters=[dict(name=v.GetName(),value=v.getVal(),error=v.getError(),low=v.getMin(),high=v.getMax()) for v in fit.floatParsFinal()] if fit else []
    checks=[];channels=[]
    for region,suffix,degree,obsname in [('SR','',5,'mass'),('CR','_CR',4,'mass_CR')]:
        for i,era in enumerate(ERAS,1):
            co=[1]+[w.var(f'local_c{j}_{era}{suffix}').getVal() for j in range(1,degree+1)]
            roots=np.polynomial.chebyshev.chebroots(np.polynomial.chebyshev.chebder(co))
            pts=[-1.,1.]+[float(z.real) for z in roots if abs(z.imag)<1e-8 and -1<z.real<1]
            check=dict(region=region,era=era,coefficients=co,minimum=float(min(np.polynomial.chebyshev.chebval(pts,co))))
            checks.append(check)
            ch=f'{region}_ch{i}';cp=w.obj('ModelConfig_bonly').GetPdf().getPdf(ch)
            pdf=next(p for p in cp.getComponents() if p.InheritsFrom('RooAddPdf'))
            x=w.var(obsname);obs=ROOT.RooArgSet(x);bw=(x.getMax()-x.getMin())/x.getBins()
            data=w.data('data_obs').reduce(ROOT.RooFit.Cut(f'CMS_channel==CMS_channel::{ch}'))
            total=pdf.expectedEvents(obs);nsum=0;term=0;points=[]
            for k in range(data.numEntries()):
                row=data.get(k);x.setVal(row.getRealValue(obsname));n=data.weight();mu=pdf.getVal(obs)*total*bw
                assert mu>0,(ch,k,mu)
                nsum+=n;term+=n*math.log(n/mu) if n else 0
                components={}
                for pp,rr in zip(pdf.pdfList(),pdf.coefList()):
                    pn=pp.GetName();key='upsilon1s' if 'upsilon1s' in pn else 'upsilon2s' if 'upsilon2s' in pn else 'dy'
                    components[key]=rr.getVal()*pp.getVal(obs)*bw
                assert abs(sum(components.values())-mu)<1e-7*max(1,mu)
                points.append(dict(mass=x.getVal(),observed=n,expected=mu,components=components))
            channels.append(dict(era=era,region=region,events=nsum,expected_events=total,range=[x.getMin(),x.getMax()],bins=x.getBins(),deviance=2*(total-nsum+term),points=points))
    return dict(parameters=parameters,positivity=checks,channels=channels,deviance=sum(c['deviance'] for c in channels),status=fit.status() if fit else None,covQual=fit.covQual() if fit else None,edm=fit.edm() if fit else None)

def main():
    D=R/'generator';D.mkdir(parents=True,exist_ok=True)
    if (D/'results.json').exists():
        print('GENERATOR ALREADY PREPARED',flush=True);return
    source=OLD/'degree3/m10.8/joint/workspace_seeded.root'
    source_hash=sha(source);attempts=[]
    for variant in [0,1,2]:
        wd=D/f'start{variant}';wd.mkdir(exist_ok=True)
        f=ROOT.TFile.Open(str(source));w=f.Get('w');keep=customize(w)
        if variant:
            for era in ERAS:
                w.var(f'local_c4_{era}').setVal(.025 if variant==1 else -.025)
                w.var(f'local_c5_{era}').setVal(-.015 if variant==1 else .015)
        w.writeToFile(str(wd/'workspace.root'));f.Close()
        cmd=[str(P/'bin/combine'),'-M','FitDiagnostics',str(wd/'workspace.root'),'-m','10.8','-n','_generator','--skipSBFit','--minos','none','--robustFit','0','--keepFailures','--cminDefaultMinimizerStrategy','2','--cminDefaultMinimizerTolerance','0.0001','--cminDefaultMinimizerPrecision','1e-12','--X-rtd','MINIMIZER_freezeDisassociatedParams','-v','0']
        (wd/'command.json').write_text(json.dumps(cmd,indent=2));start=time.time()
        with (wd/'combine.log').open('w') as log:rc=subprocess.run(cmd,cwd=wd,stdout=log,stderr=subprocess.STDOUT,timeout=600).returncode
        (wd/'receipt.json').write_text(json.dumps(dict(exit_code=rc,seconds=time.time()-start),indent=2));assert rc==0
        ff=ROOT.TFile.Open(str(wd/'fitDiagnostics_generator.root'));fit=ff.Get('fit_b');wf=ROOT.TFile.Open(str(wd/'workspace.root'));w=wf.Get('w')
        q=export(w,fit);q['variant']=variant;q['directory']=str(wd)
        q['valid']=q['status']==0 and q['covQual']==3 and min(a['minimum'] for a in q['positivity'])>0
        (wd/'results.json').write_text(json.dumps(q,indent=2));attempts.append(q)
        print('GENERATOR',variant,q['status'],q['covQual'],q['deviance'],q['valid'],flush=True)
        wf.Close();ff.Close()
    valid=[q for q in attempts if q['valid']];assert len(valid)>=2,'Generator requires two successful independent starts'
    best=min(valid,key=lambda q:q['deviance'])
    f=ROOT.TFile.Open(str(Path(best['directory'])/'workspace.root'));w=f.Get('w')
    for a in best['parameters']:w.var(a['name']).setVal(a['value'])
    means=[dict(name=v.GetName(),value=v.getVal(),constant=v.isConstant()) for v in w.allVars() if v.GetName().startswith('Ymean_')]
    assert len(means)==8 and all(a['constant'] for a in means)
    assert sum(c['events'] for c in best['channels'])==52105
    w.var('r').setVal(0);w.var('mass').setVal(11.2);w.var('mass_CR').setVal(9.8);w.writeToFile(str(D/'workspace_postfit.root'));f.Close()
    assert sha(source)==source_hash
    best.update(generator_degree=5,fit_CR_degree=4,source=str(source),source_sha256=source_hash,means=means,attempts=[{k:v for k,v in a.items() if k not in ['channels','parameters']} for a in attempts],max_valid_deviance_difference=max(a['deviance'] for a in valid)-min(a['deviance'] for a in valid),interpretation='Alternative Chebyshev5 continuum fitted B-only to the observed SR+CR, fixed μ from observed CR-only fit; generator truth conditional on this fitted model.')
    (D/'results.json').write_text(json.dumps(best,indent=2))
    print('SELECTED',best['variant'],best['deviance'],'deltaD',best['max_valid_deviance_difference'],flush=True)

if __name__=='__main__':main()
