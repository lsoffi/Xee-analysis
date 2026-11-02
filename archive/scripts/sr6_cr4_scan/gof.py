from pathlib import Path
import ROOT,json,subprocess,time,os
ROOT.gSystem.Load('libHiggsAnalysisCombinedLimit');P=Path.cwd();R=P/'runs/20261030_sr6_cr4_scan';d=R/'gof';d.mkdir(exist_ok=True);f=ROOT.TFile.Open(str(R/'m7.5/joint/workspace.root'));w=f.Get('w');ff=ROOT.TFile.Open(str(R/'m7.5/joint/fitdiagnostics/fitDiagnostics_local.root'));fit=ff.Get('fit_b');assert fit.status()==0
for v in fit.floatParsFinal():
 if w.var(v.GetName()):w.var(v.GetName()).setVal(v.getVal())
# Refit the observed data through Combine before generation, so snapshot branch selection is consistent.
w.var('r').setVal(0);w.saveSnapshot('background_fit',w.allVars());w.writeToFile(str(d/'workspace.root'));C=str(P/'bin/combine');opts=['-d',str(d/'workspace.root'),'-m','7.5','--snapshotName','background_fit','--cminDefaultMinimizerStrategy','1','--cminDefaultMinimizerTolerance','0.001','--X-rtd','MINIMIZER_freezeDisassociatedParams','--cminRunAllDiscreteCombinations']
def run(cmd,dd):
 dd.mkdir(exist_ok=True);(dd/'command.json').write_text(json.dumps(cmd));t=time.time()
 with (dd/'combine.log').open('w') as out:r=subprocess.run(cmd,cwd=dd,stdout=out,stderr=subprocess.STDOUT,timeout=7200)
 (dd/'receipt.json').write_text(json.dumps(dict(exit_code=r.returncode,seconds=time.time()-t)));print(dd.name,r.returncode,time.time()-t,flush=True)
run([C,'-M','GoodnessOfFit','--algo','saturated','--fixedSignalStrength','0','-n','_observed']+opts,d/'observed')
run([C,'-M','GenerateOnly','--saveToys','--toysFrequentist','--expectSignal','0','--genBinnedChannels',','.join(f'{region}_ch{i}' for region in ['SR','CR'] for i in range(1,5)),'-t','100','-s','306604','-n','_toys']+opts,d/'generate')
tf=next((d/'generate').glob('higgsCombine*.root'));run([C,'-M','GoodnessOfFit','--algo','saturated','--fixedSignalStrength','0','--toysFile',str(tf),'--toysFrequentist','--bypassFrequentistFit','-t','100','-s','306604','-n','_toys']+opts,d/'toys')
