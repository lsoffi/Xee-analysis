from pathlib import Path
import ROOT,json,subprocess
ROOT.gSystem.Load('libHiggsAnalysisCombinedLimit');P=Path.cwd();R=P/'runs/20261030_cheb7_bias_scan';d=R/'generator_gof';d.mkdir(exist_ok=True);f=ROOT.TFile.Open(str(R/'generator_fit/workspace.root'));w=f.Get('w');q=json.loads((R/'generator_fit/fit.json').read_text())
for k,v in q['parameters'].items():
 if w.var(k):w.var(k).setVal(v)
w.var('r').setVal(0);w.saveSnapshot('background_fit',w.allVars());w.writeToFile(str(d/'workspace.root'));C=str(P/'bin/combine');opts=['-d',str(d/'workspace.root'),'-m','7.5','--snapshotName','background_fit','--cminDefaultMinimizerStrategy','1','--cminDefaultMinimizerTolerance','0.001','--X-rtd','MINIMIZER_freezeDisassociatedParams']
def run(cmd,dd):
 dd.mkdir(exist_ok=True);(dd/'command.json').write_text(json.dumps(cmd))
 with (dd/'combine.log').open('w') as out:subprocess.run(cmd,cwd=dd,stdout=out,stderr=subprocess.STDOUT,timeout=900)
run([C,'-M','GoodnessOfFit','--algo','saturated','--fixedSignalStrength','0','-n','_observed']+opts,d/'observed')
run([C,'-M','GenerateOnly','--saveToys','--toysFrequentist','--bypassFrequentistFit','--expectSignal','0','--genBinnedChannels',','.join(f'{reg}_ch{i}' for reg in ['SR','CR'] for i in range(1,5)),'-t','100','-s','307509','-n','_toys']+opts,d/'generate');tf=next((d/'generate').glob('higgsCombine*.root'))
run([C,'-M','GoodnessOfFit','--algo','saturated','--fixedSignalStrength','0','--toysFile',str(tf),'--toysFrequentist','--bypassFrequentistFit','-t','100','-s','307509','-n','_toys']+opts,d/'toys')
