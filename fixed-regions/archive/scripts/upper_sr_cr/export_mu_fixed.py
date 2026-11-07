from pathlib import Path
import ROOT, json, math, hashlib, csv
import numpy as np
ROOT.gROOT.SetBatch(True);ROOT.gSystem.Load('libHiggsAnalysisCombinedLimit')
P=Path.cwd();R=P/'runs/20261030_upper_sr_cr_mu_fixed';eras=['2022','2022EE','2023','2023BPix'];results=[]
for degree in [1,2,3]:
    D=R/f'degree{degree}';J=D/'joint'
    f=ROOT.TFile.Open(str(J/'fitdiagnostics/fitDiagnostics_upper.root'));fit=f.Get('fit_b')
    fw=ROOT.TFile.Open(str(J/'workspace_seeded.root'));w=fw.Get('w')
    pars=[];bounds=[]
    for v in fit.floatParsFinal():
        old=w.var(v.GetName());assert old;pre=old.getVal();old.setVal(v.getVal())
        a=dict(name=v.GetName(),prefit=pre,value=v.getVal(),error=v.getError(),low=v.getMin(),high=v.getMax())
        a['near_boundary']=min(a['value']-a['low'],a['high']-a['value'])<1e-4*(a['high']-a['low'])
        if a['near_boundary']:bounds.append(a['name'])
        pars.append(a)
    w.var('r').setVal(0)
    for a in json.loads((R/'cr_only/results.json').read_text())['means']:
        v=w.var(a['name']);assert v.isConstant() and abs(v.getVal()-a['value'])<1e-12
        assert not fit.floatParsFinal().find(a['name'])
    checks=[];transfer=[]
    for era in eras:
        for region,suffix,deg in [('SR','',degree),('CR','_CR',4)]:
            coeff=[1]+[w.var(f'local_c{k}_{era}{suffix}').getVal() for k in range(1,deg+1)]
            roots=np.polynomial.chebyshev.chebroots(np.polynomial.chebyshev.chebder(coeff))
            pts=[-1,1]+[z.real for z in roots if abs(z.imag)<1e-8 and -1<z.real<1]
            checks.append(dict(era=era,region=region,minimum=float(min(np.polynomial.chebyshev.chebval(pts,coeff)))))
        for proc in ['upsilon1s','upsilon2s']:
            mu=w.var(f'Ymean_{proc}_{era}');sig=w.function(f'Ysigma_{proc}_{era}');tails=[w.function(f'Ytail_{proc}_{era}_{j}') for j in range(4)]
            xref=ROOT.RooRealVar('xref','',6,12);xref.setRange('SR',10.506666666666668,12);xref.setRange('CR',9.093333333333334,10.506666666666668)
            pdf=ROOT.RooCrystalBall('reference_Y','',xref,mu,sig,sig,*tails);obs=ROOT.RooArgSet(xref)
            for region,suffix in [('SR',''),('CR','_CR')]:
                actual=w.function(f'local_acceptance_{proc}_{era}{suffix}').getVal()
                truth=pdf.createIntegral(obs,ROOT.RooFit.NormSet(obs),ROOT.RooFit.Range(region)).getVal()
                assert abs(actual-truth)<1e-8,(degree,era,proc,region,actual,truth)
                rate=w.function(f'local_rate_{proc}_{era}{suffix}')
                assert rate.dependsOn(mu) and rate.dependsOn(w.var(f'CMS_EXO25020_bkgScale_{proc}_{era}'))
                transfer.append(dict(era=era,process=proc,region=region,fraction=actual,root_integral=truth,events=rate.getVal()))
    channels=[];data=w.data('data_obs')
    for region,varname in [('SR','mass'),('CR','mass_CR')]:
        for i,era in enumerate(eras,1):
            ch=f'{region}_ch{i}';cp=w.obj('ModelConfig_bonly').GetPdf().getPdf(ch)
            pdf=next(p for p in cp.getComponents() if p.InheritsFrom('RooAddPdf'))
            x=w.var(varname);obs=ROOT.RooArgSet(x);bw=(x.getMax()-x.getMin())/x.getBins()
            d=data.reduce(ROOT.RooFit.Cut(f'CMS_channel==CMS_channel::{ch}'));total=pdf.expectedEvents(obs)
            points=[];dev=0;nobs=0
            for k in range(d.numEntries()):
                row=d.get(k);m=row.getRealValue(varname);x.setVal(m);n=d.weight();pred=total*pdf.getVal(obs)*bw
                assert pred>0;dev+=n*math.log(n/pred) if n else 0;nobs+=n
                parts={}
                for pp,rr in zip(pdf.pdfList(),pdf.coefList()):
                    name=pp.GetName();proc='upsilon1s' if 'upsilon1s' in name else 'upsilon2s' if 'upsilon2s' in name else 'dy'
                    parts[proc]=rr.getVal()*pp.getVal(obs)*bw
                assert abs(sum(parts.values())-pred)<1e-7*max(1,pred)
                points.append(dict(mass=m,observed=n,expected=pred,components=parts))
            channels.append(dict(era=era,region=region,events=nobs,expected_events=total,range=[x.getMin(),x.getMax()],bins=x.getBins(),deviance=2*(total-nobs+dev),points=points))
    g=ROOT.TFile.Open(str(next((J/'gof_observed').glob('higgsCombine*.root'))));t=g.Get('limit');t.GetEntry(0);gof=float(t.limit)
    a=dict(degree=degree,status=fit.status(),covQual=fit.covQual(),edm=fit.edm(),nfloat=fit.floatParsFinal().getSize(),minNll=fit.minNll(),parameters=pars,boundaries=bounds,positivity=checks,transfer_checks=transfer,channels=channels,D_FD=sum(c['deviance'] for c in channels),D_GoF=gof)
    a['numerical_ok']=a['status']==0 and a['covQual']==3 and all(c['minimum']>0 for c in checks)
    assert abs(a['D_FD']-gof)<.05,(degree,a['D_FD'],gof)
    for region in ['signal','control']:
        au=json.loads((D/region/'audit.json').read_text());assert hashlib.sha256(Path(au['source']).read_bytes()).hexdigest()==au['source_sha256']
    assert sum(c['events'] for c in channels)==52105
    w.writeToFile(str(J/'workspace_postfit.root'))
    (D/'results.json').write_text(json.dumps(a,indent=2))
    with (D/'parameters.csv').open('w') as stream:
        wr=csv.DictWriter(stream,fieldnames=pars[0].keys());wr.writeheader();wr.writerows(pars)
    results.append(a);g.Close();f.Close();fw.Close()
    print(degree,a['status'],a['covQual'],a['nfloat'],a['D_FD'],a['D_GoF'],a['boundaries'],flush=True)
(R/'results.json').write_text(json.dumps(results,indent=2))
