from pathlib import Path
import ROOT,json,subprocess,time
ROOT.gSystem.Load('libHiggsAnalysisCombinedLimit');P=Path.cwd();R=P/'runs/20261030_fixed_sr_scan';d=R/'gof';d.mkdir(exist_ok=True);f=ROOT.TFile.Open(str(R/'m6.5/joint/workspace.root'));w=f.Get('w');q=json.loads((R/'m6.5/joint/quality.json').read_text())
for p in q['b']['parameters']:
 if w.var(p['name']):w.var(p['name']).setVal(p['value'])
w.var('r').setVal(0);w.saveSnapshot('background_fit',w.allVars());w.writeToFile(str(d/'workspace.root'));C=str(P/'bin/combine');opts=['-d',str(d/'workspace.root'),'-m','6.5','--snapshotName','background_fit','--cminDefaultMinimizerStrategy','1','--cminDefaultMinimizerTolerance','0.001','--X-rtd','MINIMIZER_freezeDisassociatedParams']
def run(cmd,dd):
 dd.mkdir(exist_ok=True);(dd/'command.json').write_text(json.dumps(cmd));t=time.time()
 with (dd/'combine.log').open('w') as out:r=subprocess.run(cmd,cwd=dd,stdout=out,stderr=subprocess.STDOUT,timeout=600)
 (dd/'receipt.json').write_text(json.dumps(dict(exit_code=r.returncode,seconds=time.time()-t)))
run([C,'-M','GoodnessOfFit','--algo','saturated','--fixedSignalStrength','0','-n','_observed']+opts,d/'observed')
run([C,'-M','GenerateOnly','--saveToys','--toysFrequentist','--bypassFrequentistFit','--expectSignal','0','--genBinnedChannels',','.join(f'{region}_ch{i}' for region in ['SR','CR'] for i in range(1,5)),'-t','100','-s','306509','-n','_toys']+opts,d/'generate')
tf=next((d/'generate').glob('higgsCombine*.root'));run([C,'-M','GoodnessOfFit','--algo','saturated','--fixedSignalStrength','0','--toysFile',str(tf),'--toysFrequentist','--bypassFrequentistFit','-t','100','-s','306509','-n','_toys']+opts,d/'toys')
