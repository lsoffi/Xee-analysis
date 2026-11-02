from pathlib import Path
import json, subprocess, ROOT
ROOT.gSystem.Load('libHiggsAnalysisCombinedLimit')
P=Path.cwd();R=P/'runs/20261030_upper_sr_cr';out=[]
for degree in [1,2,3]:
    J=R/f'degree{degree}/joint'
    primary=json.loads((R/f'degree{degree}/results.json').read_text())
    for variant in [1,2]:
        D=J/f'multistart{variant}';D.mkdir(exist_ok=True)
        f=ROOT.TFile.Open(str(J/'workspace_seeded.root'));w=f.Get('w')
        for v in w.allVars():
            name=v.GetName()
            if name.startswith('Ymean_'):v.setVal((v.getMin()+v.getMax())/2 + (.08 if variant==1 else -.08))
            if 'bkgScale_upsilon2s' in name:v.setVal(.5 if variant==1 else 1.5)
            if name.startswith('local_c') and not name.endswith('_CR'):v.setVal(0 if variant==1 else -.65 if name.startswith('local_c1_') else 0)
        w.writeToFile(str(D/'start.root'));f.Close()
        cmd=json.loads((J/'fitdiagnostics/combine.command.json').read_text());cmd[cmd.index('-d')+1]=str(D/'start.root')
        (D/'command.json').write_text(json.dumps(cmd,indent=2))
        with (D/'combine.log').open('w') as h:res=subprocess.run(cmd,cwd=D,stdout=h,stderr=subprocess.STDOUT,timeout=180)
        ff=ROOT.TFile.Open(str(D/'fitDiagnostics_upper.root'));fit=ff.Get('fit_b')
        a=dict(degree=degree,start=variant,exit_code=res.returncode,status=fit.status(),covQual=fit.covQual(),deltaNLL=fit.minNll()-primary['minNll']);out.append(a);print(a,flush=True);ff.Close()
(R/'multistart.json').write_text(json.dumps(out,indent=2))
