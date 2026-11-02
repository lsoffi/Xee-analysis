"""CR-only mean estimates followed by simultaneous fits with the eight means fixed."""
from pathlib import Path
import ROOT, json, subprocess, time, shutil, math
ROOT.gROOT.SetBatch(True);ROOT.gSystem.Load('libHiggsAnalysisCombinedLimit')
P=Path.cwd();OLD=P/'runs/20261030_upper_sr_cr';R=P/'runs/20261030_upper_sr_cr_mu_fixed';R.mkdir(exist_ok=True)
CR=R/'cr_only';CR.mkdir(exist_ok=True)
base=json.loads((OLD/'degree2/joint/fitdiagnostics/combine.command.json').read_text())
def execute(cmd,d):
    (d/'command.json').write_text(json.dumps(cmd,indent=2));t=time.time()
    with (d/'combine.log').open('w') as h:rc=subprocess.run(cmd,cwd=d,stdout=h,stderr=subprocess.STDOUT,timeout=300).returncode
    (d/'receipt.json').write_text(json.dumps(dict(exit_code=rc,seconds=time.time()-t)))
    assert rc==0
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
            data=w.data('data_obs').reduce(ROOT.RooFit.Cut(f'CMS_channel==CMS_channel::{ch}'));nsum=0;term=0
            for k in range(data.numEntries()):
                row=data.get(k);x.setVal(row.getRealValue(obsname));n=data.weight();mu=pdf.getVal(obs)*total*bw;nsum+=n;term+=n*math.log(n/mu) if n else 0
            D+=2*(total-nsum+term)
    a=dict(directory=str(d),status=fit.status(),covQual=fit.covQual(),minNll=fit.minNll(),deviance=D,edm=fit.edm(),parameters=params);wf.Close();f.Close();return a
# The only data used here are CR bins. All four eras are fitted simultaneously.
seed=json.loads((OLD/'degree2/results.json').read_text())['parameters'];attempts=[]
for variant in [0,1,2]:
    d=CR/f'start{variant}';d.mkdir(exist_ok=True)
    f=ROOT.TFile.Open(str(OLD/'degree2/control/workspace.root'));w=f.Get('w')
    for a in seed:
        name=a['name'];target=name[:-3] if name.endswith('_CR') else name
        if (name.endswith('_CR') or name.startswith('Ymean_') or 'bkgScale_upsilon' in name) and w.var(target):w.var(target).setVal(a['value'])
    if variant:
        for v in w.allVars():
            name=v.GetName()
            if name.startswith('Ymean_'):v.setVal((v.getMin()+v.getMax())/2+(.07 if variant==1 else -.07))
            if 'bkgScale_upsilon2s' in name:v.setVal(.5 if variant==1 else 1.5)
    assert abs(w.data('data_obs').sumEntries()-48038)<1e-8
    assert w.var('mass').getMin()==9.093333333333334 and w.var('mass').getMax()==10.506666666666668
    w.var('r').setVal(0);w.writeToFile(str(d/'workspace.root'));f.Close()
    cmd=base.copy();cmd[cmd.index('-d')+1]=str(d/'workspace.root');execute(cmd,d);a=readfit(d);attempts.append(a)
    print('CR',variant,a['status'],a['covQual'],a['deviance'],flush=True)
best=min((a for a in attempts if a['status']==0 and a['covQual']==3),key=lambda a:a['deviance'])
means=[a for a in best['parameters'] if a['name'].startswith('Ymean_')];assert len(means)==8
for a in means:a['near_boundary']=min(a['value']-a['low'],a['high']-a['value'])<1e-4*(a['high']-a['low'])
(CR/'results.json').write_text(json.dumps(dict(selected=best,attempts=attempts,means=means,events=48038,range=[9.093333333333334,10.506666666666668]),indent=2))
print('CR MEANS',json.dumps(means),flush=True)
for degree in [1,2,3]:
    D=R/f'degree{degree}';J=D/'joint';J.mkdir(parents=True,exist_ok=True)
    for region in ['signal','control']:
        (D/region).mkdir(exist_ok=True);shutil.copy2(OLD/f'degree{degree}'/region/'audit.json',D/region/'audit.json')
    f=ROOT.TFile.Open(str(OLD/f'degree{degree}/joint/workspace_postfit.root'));w=f.Get('w')
    # Reset CR nuisance and yields to the CR-only solution as initial values.
    for a in best['parameters']:
        name=a['name'];target=name+'_CR' if name.startswith('local_') else name
        if w.var(target):w.var(target).setVal(a['value'])
    for a in means:
        v=w.var(a['name']);assert v;v.setVal(a['value']);v.setConstant(True)
    for era in ['2022','2022EE','2023','2023BPix']:
        for proc in ['upsilon1s','upsilon2s']:assert not w.var(f'CMS_EXO25020_bkgScale_{proc}_{era}').isConstant()
    w.var('r').setVal(0);w.writeToFile(str(J/'workspace_seeded.root'));f.Close()
    fd=J/'fitdiagnostics';fd.mkdir(exist_ok=True)
    cmd=base.copy();cmd[cmd.index('-d')+1]=str(J/'workspace_seeded.root')
    cmd+=['--freezeParameters',','.join(a['name'] for a in means)]
    execute(cmd,fd);a=readfit(fd)
    assert all(not p['name'].startswith('Ymean_') for p in a['parameters'])
    print('JOINT FIXED',degree,a['status'],a['covQual'],a['deviance'],flush=True)
    gd=J/'gof_observed';gd.mkdir(exist_ok=True)
    gc=json.loads((OLD/f'degree{degree}/joint/gof_observed/combine.command.json').read_text());gc[gc.index('-d')+1]=str(J/'workspace_seeded.root');gc+=['--freezeParameters',','.join(a['name'] for a in means)];execute(gc,gd)
config=dict(status='completed',source_strategy='CR-only means frozen in simultaneous fit, per user correction',SR_range=[10.506666666666668,12],CR_range=[9.093333333333334,10.506666666666668],SR_degrees=[1,2,3],CR_degree=4,means=means,means_source=str(CR/'results.json'),upsilon_normalizations='free and shared per era between SR and CR',continuum='free independently per era and region',upsilon_sigma_tails='fixed MC',fit='B-only',mean_uncertainty_propagated=False,prior_from_CR_fit=False,gof_toys=0)
(R/'configuration.json').write_text(json.dumps(config,indent=2))
