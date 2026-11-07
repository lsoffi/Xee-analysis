from pathlib import Path
import ROOT, math, json
def readfit(d):
    f=ROOT.TFile.Open(str(d/'fitDiagnostics_upper.root'));fit=f.Get('fit_b')
    params=[]
    for v in fit.floatParsFinal():
        params.append(dict(name=v.GetName(),value=v.getVal(),error=v.getError(),low=v.getMin(),high=v.getMax()))
    source=d/'workspace.root' if (d/'workspace.root').exists() else d.parent/'workspace_seeded.root'
    wf=ROOT.TFile.Open(str(source));w=wf.Get('w')
    for a in params:w.var(a['name']).setVal(a['value'])
    w.var('r').setVal(0);D=0
    regions=[('SR_','mass'),('CR_','mass_CR')] if w.var('mass_CR') else [('', 'mass')]
    for prefix,obsname in regions:
        for i in range(1,5):
            ch=prefix+f'ch{i}';cp=w.obj('ModelConfig_bonly').GetPdf().getPdf(ch)
            pdf=next(p for p in cp.getComponents() if p.InheritsFrom('RooAddPdf'))
            x=w.var(obsname);obs=ROOT.RooArgSet(x);bw=(x.getMax()-x.getMin())/x.getBins();total=pdf.expectedEvents(obs)
            data=(w.data('toy_data') or w.data('data_obs')).reduce(ROOT.RooFit.Cut(f'CMS_channel==CMS_channel::{ch}'));nsum=0;term=0
            for k in range(data.numEntries()):
                row=data.get(k);x.setVal(row.getRealValue(obsname));n=data.weight();mu=pdf.getVal(obs)*total*bw;nsum+=n;term+=n*math.log(n/mu) if n else 0
            D+=2*(total-nsum+term)
    a=dict(directory=str(d),status=fit.status(),covQual=fit.covQual(),minNll=fit.minNll(),deviance=D,edm=fit.edm(),parameters=params);wf.Close();f.Close();return a
