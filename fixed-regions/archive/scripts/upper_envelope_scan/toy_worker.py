from pathlib import Path
import ROOT,json,sys,subprocess,time,math,fcntl
from fit_helpers import readfit
ROOT.gROOT.SetBatch(True);ROOT.gSystem.Load('libHiggsAnalysisCombinedLimit')
P=Path.cwd();R=P/'runs/20261030_upper_envelope_scan';mode=sys.argv[1];degree=int(sys.argv[2]);key=sys.argv[3];C=str(P/'bin/combine')
if mode=='gof':
    toyid=int(key);D=R/f'gof/degree{degree}/toy{toyid:03d}';G=R/f'gof/degree{degree}/generation';masses=sys.argv[4:] or [f'{k/10:.1f}' for k in range(105,121)];toyname=f'toy_{toyid}'
else:
    mass=key;label=sys.argv[4];D=R/f'bias/degree{degree}/m{mass}/{label}';G=D/'generation';masses=[mass];toyname='toy_asimov'
D.mkdir(parents=True,exist_ok=True);lockfile=(D/'worker.lock').open('a');fcntl.flock(lockfile,fcntl.LOCK_EX);truth=json.loads((G/'truth.json').read_text());tf=ROOT.TFile.Open(str(next(G.glob('higgsCombine*.root'))));toy=tf.Get('toys/'+toyname);snap=tf.Get('toys/'+toyname+'_snapshot');assert toy and snap
def execute(cmd,d):
    d.mkdir(parents=True,exist_ok=True)
    if (d/'receipt.json').exists() and json.loads((d/'receipt.json').read_text())['exit_code']==0:return
    (d/'command.json').write_text(json.dumps(cmd,indent=2));t=time.time()
    with (d/'combine.log').open('w') as h:rc=subprocess.run(cmd,cwd=d,stdout=h,stderr=subprocess.STDOUT,timeout=300).returncode
    (d/'receipt.json').write_text(json.dumps(dict(exit_code=rc,seconds=time.time()-t)));assert rc==0,(d,rc)
crresult=D/'cr_result.json'
if not crresult.exists():
    attempts=[];ref=json.loads((P/'runs/20261030_upper_sr_cr_mu_fixed/cr_only/results.json').read_text())['selected']['parameters']
    for variant in [0,1]:
        wd=D/f'cr_start{variant}';wd.mkdir(exist_ok=True);f=ROOT.TFile.Open(str(R/'degree2/m10.8/control/workspace.root'));w=f.Get('w')
        for a in ref:
            if w.var(a['name']):w.var(a['name']).setVal(a['value'])
        if variant:
            for v in w.allVars():
                if v.GetName().startswith('Ymean_'):v.setVal((v.getMin()+v.getMax())/2-.07)
                if 'bkgScale_upsilon2s' in v.GetName():v.setVal(1.5)
        weight=ROOT.RooRealVar('toy_weight','',0,1e12);obs=ROOT.RooArgSet(w.var('mass'),w.cat('CMS_channel'),weight);data=ROOT.RooDataSet('toy_data','',obs,ROOT.RooFit.WeightVar(weight));nt=0
        for k in range(toy.numEntries()):
            row=toy.get(k);ch=row.find('CMS_channel').getLabel()
            if not ch.startswith('CR_'):continue
            w.cat('CMS_channel').setLabel(ch[3:]);w.var('mass').setVal(row.getRealValue('mass_CR'));weight.setVal(toy.weight());data.add(obs,toy.weight());nt+=toy.weight()
        getattr(w,'import')(data);w.allVars().assignValueOnly(snap);w.var('r').setVal(0);w.writeToFile(str(wd/'workspace.root'));f.Close()
        cmd=[C,'-M','FitDiagnostics',str(wd/'workspace.root'),'-D','toy_data','-n','_upper','--skipSBFit','--minos','none','--robustFit','0','--keepFailures','--cminDefaultMinimizerStrategy','2','--cminDefaultMinimizerTolerance','0.001','--cminDefaultMinimizerPrecision','1e-12','--X-rtd','MINIMIZER_freezeDisassociatedParams','-v','0'];execute(cmd,wd);q=readfit(wd);q['variant']=variant;attempts.append(q)
    converged=[q for q in attempts if q['status']==0];best=min(converged or attempts,key=lambda a:a['deviance']);means=[a for a in best['parameters'] if a['name'].startswith('Ymean_')];assert len(means)==8
    crresult.write_text(json.dumps(dict(selected=best,means=means,attempts=attempts,events=nt,valid=best['status']==0),indent=2))
cr=json.loads(crresult.read_text())
for mass in masses:
    wd=D/f'm{mass}';wd.mkdir(exist_ok=True)
    if (wd/'result.json').exists():continue
    source=R/f'degree{degree}/m{mass}/joint/workspace_seeded.root'
    while not source.exists():time.sleep(5)
    f=ROOT.TFile.Open(str(source));w=f.Get('w');td=toy.Clone('toy_data');getattr(w,'import')(td);w.allVars().assignValueOnly(snap)
    for a in cr['selected']['parameters']:
        target=a['name']+'_CR' if a['name'].startswith('local_') else a['name']
        if w.var(target) and not target.startswith('Ymean_'):w.var(target).setVal(a['value'])
    for a in cr['means']:w.var(a['name']).setVal(a['value']);w.var(a['name']).setConstant(True)
    rtrue=truth['rtrue'];rmax=max(100,5*rtrue);w.var('r').setRange(-rmax,rmax);w.var('r').setVal(rtrue);w.var('mass').setVal(11.2);w.writeToFile(str(wd/'workspace.root'));f.Close()
    common=[str(wd/'workspace.root'),'-D','toy_data','-m',mass,'-n','_upper','--rMin',str(-rmax),'--rMax',str(rmax),'--cminDefaultMinimizerStrategy','0' if mode=='gof' else '2','--cminDefaultMinimizerTolerance','0.001','--cminDefaultMinimizerPrecision','1e-12','--X-rtd','MINIMIZER_freezeDisassociatedParams','--cminRunAllDiscreteCombinations','-v','0']
    if mode=='gof':
        execute([C,'-M','GoodnessOfFit']+common+['--algo','saturated','--trackParameters','r'],wd)
        ff=ROOT.TFile.Open(str(next(wd.glob('higgsCombine*.root'))));tr=ff.Get('limit');assert tr.GetEntries()==1;tr.GetEntry(0);score=float(tr.limit);r=float(tr.trackedParam_r);ff.Close();assert math.isfinite(score)
        q=dict(mode=mode,degree=degree,mass=float(mass),toy=toyid,score=score,r=r,r_boundary=abs(r)>.99*rmax,CR_valid=cr['valid'],CR_status=cr['selected']['status'],CR_covQual=cr['selected']['covQual'])
    else:
        execute([C,'-M','FitDiagnostics']+common+['--skipBOnlyFit','--minos','none','--robustFit','0','--keepFailures','--preFitValue',str(rtrue)],wd)
        ff=ROOT.TFile.Open(str(wd/'fitDiagnostics_upper.root'));fit=ff.Get('fit_s');v=fit.floatParsFinal().find('r');q=dict(mode=mode,degree=degree,mass=float(mass),label=label,rtrue=rtrue,r=v.getVal(),error=v.getError(),status=fit.status(),covQual=fit.covQual(),CR_valid=cr['valid'],CR_status=cr['selected']['status'],CR_covQual=cr['selected']['covQual']);q['bias']=q['r']-rtrue;q['bias_sigma']=q['bias']/q['error'] if q['error']>0 else None;q['valid']=q['status']==0 and q['covQual']==3 and q['error']>0 and abs(q['r'])<.99*rmax and cr['valid'];ff.Close()
    (wd/'result.json').write_text(json.dumps(q,indent=2));print(mode,degree,key,mass,q,flush=True)
