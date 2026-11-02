"""Targeted FD retries for new degree3 only; no expected-band recomputation."""
from pathlib import Path
import ROOT,json,subprocess,shutil,math,time
ROOT.gROOT.SetBatch(True);ROOT.gSystem.Load('libHiggsAnalysisCombinedLimit')
P=Path.cwd();R=P/'runs/20261030_upper_envelope_scan/degree3';S=P/'scripts/upper_envelope_scan'

def quality(path,workspace):
    f=ROOT.TFile.Open(str(path/'fitDiagnostics_upper.root'));fit=f.Get('fit_s');t=f.Get('tree_fit_sb');t.GetEntry(0)
    a=dict(status=fit.status(),covQual=fit.covQual(),r=fit.floatParsFinal().find('r').getVal(),error=fit.floatParsFinal().find('r').getError(),nll_nll0=float(t.nll_nll0))
    b=f.Get('fit_b');wf=ROOT.TFile.Open(str(workspace));w=wf.Get('w')
    for v in b.floatParsFinal():w.var(v.GetName()).setVal(v.getVal())
    for era in ['2022','2022EE','2023','2023BPix']:
        bp=b.floatParsFinal().find('local_b1_'+era);cp=b.floatParsFinal().find('local_c1_'+era)
        idx=1 if bp and bp.getError()>0 and (not cp or cp.getError()==0) else 0;w.cat('pdfindex_'+era).setIndex(idx)
    w.var('r').setVal(0);deviance=0;pred=[]
    for region,obsname in [('SR','mass'),('CR','mass_CR')]:
        for i in range(1,5):
            ch=f'{region}_ch{i}';cp=w.obj('ModelConfig_bonly').GetPdf().getPdf(ch);pdf=next(q for q in cp.getComponents() if q.InheritsFrom('RooAddPdf'));x=w.var(obsname);obs=ROOT.RooArgSet(x);bw=(x.getMax()-x.getMin())/x.getBins();total=pdf.expectedEvents(obs);data=w.data('data_obs').reduce(ROOT.RooFit.Cut(f'CMS_channel==CMS_channel::{ch}'));s=0;term=0
            for j in range(data.numEntries()):
                row=data.get(j);x.setVal(row.getRealValue(obsname));n=data.weight();mu=pdf.getVal(obs)*total*bw;pred.append(mu);s+=n;term+=n*math.log(n/mu) if n else 0
            deviance+=2*(total-s+term)
    constrained={v.GetName():v.getVal() for v in b.floatParsFinal() if v.GetName().startswith(('CMS_eff_','CMS_scale_','CMS_EXO25020_signalModel','lumi_'))}
    a.update(B_status=b.status(),B_covQual=b.covQual(),B_deviance=deviance,B_prediction=pred,B_constrained=constrained);wf.Close();f.Close();return a

for result in sorted(R.glob('m*/results.json')):
    a=json.loads(result.read_text());J=result.parent/'joint'
    if a['numerical_ok'] or (J/'refinement_degree3.json').exists():continue
    original=quality(J/'fitdiagnostics',J/'workspace_seeded.root');trials=[]
    for start in [1,10]:
        D=J/f'fd_start{start}';D.mkdir(exist_ok=True);cmd=json.loads((J/'fitdiagnostics/combine.command.json').read_text())
        if '--preFitValue' in cmd:cmd[cmd.index('--preFitValue')+1]=str(start)
        else:cmd+=['--preFitValue',str(start)]
        if '--initFromBonly' not in cmd:cmd+=['--initFromBonly']
        (D/'command.json').write_text(json.dumps(cmd,indent=2))
        start_time=time.time()
        with (D/'combine.log').open('w') as h:rc=subprocess.run(cmd,cwd=D,stdout=h,stderr=subprocess.STDOUT,timeout=300).returncode
        (D/'combine.receipt.json').write_text(json.dumps(dict(exit_code=rc,seconds=time.time()-start_time)))
        q=quality(D,J/'workspace_seeded.root');q.update(directory=str(D),start=start,exit_code=rc)
        constraints_match=q['B_constrained'].keys()==original['B_constrained'].keys() and all(abs(v-original['B_constrained'][name])<1e-6 for name,v in q['B_constrained'].items())
        q['B_reference_matches']=constraints_match and q['B_status']==0 and abs(q['B_deviance']-original['B_deviance'])<1e-4 and max(abs(x-y)/max(1,abs(x)) for x,y in zip(q['B_prediction'],original['B_prediction']))<1e-4
        trials.append(q)
    comparable=[q for q in trials if q['B_reference_matches']];best_nll=min([original['nll_nll0']]+[q['nll_nll0'] for q in comparable])
    candidates=[q for q in comparable if q['exit_code']==0 and q['status']==0 and q['covQual']==3 and q['nll_nll0']<=best_nll+1e-5]
    promoted=None
    if candidates:
        best=min(candidates,key=lambda q:q['nll_nll0']);D=Path(best['directory']);(J/'fitdiagnostics').rename(J/'fitdiagnostics_initial_degree3');D.rename(J/'fitdiagnostics');shutil.copy2(J/'fitdiagnostics/command.json',J/'fitdiagnostics/combine.command.json');promoted=best
        with (result.parent/'extract_refined_degree3.log').open('w') as h:subprocess.run([str(P/'bin/python'),str(S/'extract.py'),f'{a["mass"]:.1f}','3'],cwd=P,stdout=h,stderr=subprocess.STDOUT,check=True)
    for q in [original]+trials:q.pop('B_prediction',None)
    (J/'refinement_degree3.json').write_text(json.dumps(dict(original=original,trials=trials,promoted=promoted),indent=2));print(a['mass'],'PROMOTED',bool(promoted),[(q['status'],q['covQual'],q['B_reference_matches'],q['nll_nll0']) for q in trials],flush=True)
