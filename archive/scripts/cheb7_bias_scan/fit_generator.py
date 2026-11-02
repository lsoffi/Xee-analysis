from pathlib import Path
import ROOT,json,subprocess
ROOT.gSystem.Load('libHiggsAnalysisCombinedLimit');P=Path.cwd();R=P/'runs/20261030_cheb7_bias_scan';d=R/'generator_fit';d.mkdir(exist_ok=True);f=ROOT.TFile.Open(str(R/'m7.5/joint/workspace.root'));w=f.Get('w');q=json.loads((P/'runs/20261030_fixed_sr_scan/m7.5/joint/quality.json').read_text())
for p in q['b']['parameters']:
 if w.var(p['name']):w.var(p['name']).setVal(p['value'])
w.writeToFile(str(d/'workspace.root'));cmd=[str(P/'bin/combine'),'-M','FitDiagnostics','-d',str(d/'workspace.root'),'-m','7.5','--skipSBFit','--saveShapes','--saveNormalizations','--keepFailures','--minos','none','--cminDefaultMinimizerStrategy','1','--cminDefaultMinimizerTolerance','0.001','--X-rtd','MINIMIZER_freezeDisassociatedParams','-n','_generator','-v','1'];(d/'command.json').write_text(json.dumps(cmd,indent=2))
with (d/'combine.log').open('w') as out:subprocess.run(cmd,cwd=d,stdout=out,stderr=subprocess.STDOUT)
f.Close()
