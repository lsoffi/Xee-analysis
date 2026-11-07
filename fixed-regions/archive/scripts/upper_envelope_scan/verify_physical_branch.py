from pathlib import Path
import json,subprocess,ROOT,math
ROOT.gSystem.Load('libHiggsAnalysisCombinedLimit')
P=Path.cwd();R=P/'runs/20261030_upper_envelope_scan'
for degree,mass in [(1,11.0),(1,11.7),(2,12.0)]:
    J=R/f'degree{degree}/m{mass:.1f}/joint';D=J/'physical_branch_start10';D.mkdir(exist_ok=True)
    cmd=json.loads((J/'fitdiagnostics/combine.command.json').read_text())
    cmd[cmd.index('--rMin')+1]='0'
    for opt in ['--saveShapes','--saveNormalizations']:
        if opt in cmd:cmd.remove(opt)
    if '--preFitValue' in cmd:cmd[cmd.index('--preFitValue')+1]='10'
    else:cmd+=['--preFitValue','10']
    if '--initFromBonly' in cmd:cmd.remove('--initFromBonly')
    (D/'command.json').write_text(json.dumps(cmd,indent=2))
    with (D/'combine.log').open('w') as h:rc=subprocess.run(cmd,cwd=D,stdout=h,stderr=subprocess.STDOUT,timeout=600).returncode
    assert rc==0,(D,rc)
    f=ROOT.TFile.Open(str(D/'fitDiagnostics_upper.root'));fit=f.Get('fit_s');t=f.Get('tree_fit_sb');t.GetEntry(0)
    r=fit.floatParsFinal().find('r');z=math.sqrt(max(0,-2*float(t.nll_nll0))) if r.getVal()>0 else 0
    obs=json.loads((J.parent/'results.json').read_text())
    q=dict(mass=mass,degree=degree,exit_code=rc,r=r.getVal(),error=r.getError(),status=fit.status(),covQual=fit.covQual(),nll_nll0=float(t.nll_nll0),Z_physical=z,Z_significance=obs['Z'],consistent=abs(z-obs['Z'])<.05,note='Independent FD restricted to r >= 0; unconstrained FD retained separately. Multimodal likelihood flagged.')
    (D/'results.json').write_text(json.dumps(q,indent=2));print(q,flush=True);f.Close()
