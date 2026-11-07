"""Read-only pilot check of fixed-zero versus injected-r initialisation."""
from pathlib import Path
import json,ROOT
from worker_helpers import execute,read_fit,atomic_json
ROOT.gROOT.SetBatch(True);ROOT.gSystem.Load('libHiggsAnalysisCombinedLimit')
P=Path(__file__).resolve().parents[2];R=P/'runs/20261030_upper_bias_cheb5_toys';results=[]
for m,label in [('10.8','r0'),('10.5','rExpected'),('12.0','rExpected')]:
 source=R/f'pilot_profile/m{m}/{label}/toy001/result.json';q=json.loads(source.read_text());d=source.parent/'zero_seed_audit';d.mkdir(exist_ok=True)
 f=ROOT.TFile.Open(str(Path(q['directory'])/'workspace.root'));w=f.Get('w');w.var('r').setVal(0);w.writeToFile(str(d/'workspace.root'));f.Close()
 cmd=[str(P/'bin/combine'),'-M','FitDiagnostics',str(d/'workspace.root'),'-D','toy_data','-m',m,'-n','_upper','--skipBOnlyFit','--preFitValue','0','--minos','none','--robustFit','0','--keepFailures','--cminDefaultMinimizerStrategy','2','--cminDefaultMinimizerTolerance','0.00001','--cminDefaultMinimizerPrecision','1e-12','--cminRunAllDiscreteCombinations','--X-rtd','MINIMIZER_freezeDisassociatedParams','-v','0']
 receipt=execute(cmd,d,180);fit=read_fit(d,True,'hesse68');fit.update(receipt=receipt,source=str(source),original_r=q['r'],original_deviance_with_aux=q['deviance_with_aux']);fit['delta_r']=fit.get('r',0)-q['r'];fit['delta_deviance_with_aux']=fit.get('deviance_with_aux',0)-q['deviance_with_aux'];atomic_json(d/'audit.json',fit);print(m,fit['delta_r'],fit['delta_deviance_with_aux'],fit.get('status'),flush=True);results.append(fit)
atomic_json(R/'zero_seed_pilot_audit.json',results)
