from pathlib import Path
import ROOT,json,subprocess,shutil
ROOT.gSystem.Load('libHiggsAnalysisCombinedLimit');P=Path.cwd();R=P/'runs/20261030_upper_envelope_scan';S=P/'scripts/upper_envelope_scan'
def quality(path):
    f=ROOT.TFile.Open(str(path/'fitDiagnostics_upper.root'));fit=f.Get('fit_s');t=f.Get('tree_fit_sb');t.GetEntry(0)
    a=dict(status=fit.status(),covQual=fit.covQual(),r=fit.floatParsFinal().find('r').getVal(),error=fit.floatParsFinal().find('r').getError(),nll_nll0=float(t.nll_nll0));f.Close();return a
for result in sorted(R.glob('degree*/m*/results.json')):
    a=json.loads(result.read_text());J=result.parent/'joint'
    if a['numerical_ok'] or (J/'refinement.json').exists():continue
    original=quality(J/'fitdiagnostics');trials=[]
    for start in [1,10]:
        D=J/f'fd_start{start}';D.mkdir(exist_ok=True);cmd=json.loads((J/'fitdiagnostics/combine.command.json').read_text());cmd+=['--preFitValue',str(start),'--initFromBonly'];(D/'command.json').write_text(json.dumps(cmd,indent=2))
        with (D/'combine.log').open('w') as h:rc=subprocess.run(cmd,cwd=D,stdout=h,stderr=subprocess.STDOUT,timeout=300).returncode
        q=quality(D);q.update(directory=str(D),start=start,exit_code=rc);trials.append(q)
    candidates=[q for q in trials if q['status']==0 and q['covQual']==3 and q['nll_nll0']<=original['nll_nll0']+1e-5]
    promoted=None
    if candidates:
        best=min(candidates,key=lambda q:q['nll_nll0']);D=Path(best['directory']);(J/'fitdiagnostics').rename(J/'fitdiagnostics_initial');D.rename(J/'fitdiagnostics');shutil.copy2(J/'fitdiagnostics/command.json',J/'fitdiagnostics/combine.command.json');promoted=best
        with (result.parent/'extract_refined.log').open('w') as h:subprocess.run([str(P/'bin/python'),str(S/'extract.py'),f'{a["mass"]:.1f}',str(a['degree'])],cwd=P,stdout=h,stderr=subprocess.STDOUT,check=True)
    (J/'refinement.json').write_text(json.dumps(dict(original=original,trials=trials,promoted=promoted),indent=2));print(a['degree'],a['mass'],original,trials,'PROMOTED',bool(promoted),flush=True)
