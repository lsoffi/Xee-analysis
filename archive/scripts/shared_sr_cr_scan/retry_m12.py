"""Bounded numerical continuation for the failed 12-GeV observed S+B fit.

No change to the model, coefficient domains, positivity requirement or r bounds.
Start from the interior B-only snapshot, using r=0 and the prior finite best r.
Try strategy1 then, only after a failed point fit, strategy2 at tighter tolerance.
All original fits and toy input templates remain untouched.
"""
from pathlib import Path
import json, math, sys
import ROOT
from run_observed import execute, common, C, P, R
from observed_helpers import read_stage, dump

ROOT.gROOT.SetBatch(True)
ROOT.gSystem.Load('libHiggsAnalysisCombinedLimit')
mass='12.0';D=R/'observed'/('m'+mass);J=D/'joint'
prior=json.loads((D/'results.json').read_text())
protocol=json.loads((D/'observed_protocol.json').read_text())
bound=max(abs(v) for v in protocol['initial_r_range'])
assert bound>=100
start_selected=prior['fit']['s']['r']
err=prior['fit']['s']['error']
if not math.isfinite(err) or err<=0:err=max(1,abs(start_selected))
source=J/'workspace_postfit_b.root'
trials=[]
for label,initial in [('zero',0.),('selected',start_selected)]:
    for strategy,tolerance in [(1,.001),(2,.00001)]:
        directory=J/f'fitdiagnostics_anchor_retry_{label}'
        if strategy==2:directory=J/f'fitdiagnostics_anchor_retry_{label}_strategy2'
        directory.mkdir(exist_ok=True)
        path=directory/'workspace.root'
        if not (directory/'combine.receipt.json').exists():
            f=ROOT.TFile.Open(str(source));w=f.Get('w')
            w.var('r').setRange(-bound,bound);w.var('r').setVal(initial);w.var('r').setError(err)
            w.writeToFile(str(path));f.Close()
        command=[C,'-M','FitDiagnostics']+common(path,mass,strategy,tolerance)+[
            '--saveShapes','--saveNormalizations','--keepFailures','--robustFit','0',
            '--minos','none','--skipBOnlyFit','--rMin',str(-bound),'--rMax',str(bound),
            '--preFitValue',str(initial)]
        receipt=execute(command,directory,timeout=1200)
        result,bins=read_stage(J/'workspace_seeded.root',directory/'fitDiagnostics_shared.root','s')
        result.update(directory=str(directory),receipt=receipt,initial_r=initial,initial_r_error=err,
                      strategy=strategy,tolerance=tolerance,initial_snapshot=str(source))
        trials.append(result);dump(directory/'retry_result.json',result)
        print(label,strategy,result.get('r'),result.get('status'),result.get('covQual'),result.get('edm'),result.get('exact_deviance'),flush=True)
        if result.get('point_valid') and receipt['exit_code']==0:break
dump(D/'numerical_retry_summary.json',dict(reason='Initial two S+B fits status3; active continuum positivity boundary at12GeV in2022EE',
    unchanged_r_range=[-bound,bound],max_attempts=4,trials=trials,all_original_results_preserved=True))
