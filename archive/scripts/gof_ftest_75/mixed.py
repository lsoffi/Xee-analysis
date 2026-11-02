from pathlib import Path
exec((Path('scripts/gof_ftest_75/run.py').read_text()).split('for n in [2,3]:')[0])
d=R/'degree2SR_3CR';w=str(d/'workspace.root');fd=d/'fitdiagnostics'
run([C,'-M','FitDiagnostics','-d',w,'--skipSBFit','--saveShapes','--saveNormalizations','--robustFit','0','--minos','none','-n','_mixed']+opts,fd)
code="""
import ROOT
from pathlib import Path
ROOT.gSystem.Load('libHiggsAnalysisCombinedLimit');d=Path('runs/20261029_gof_ftest_75/degree2SR_3CR');f=ROOT.TFile.Open(str(d/'workspace.root'));w=f.Get('w');ff=ROOT.TFile.Open(str(d/'fitdiagnostics/fitDiagnostics_mixed.root'));fit=ff.Get('fit_b');assert fit.status()==0 and fit.covQual()==3
for v in fit.floatParsFinal():
 if w.var(v.GetName()):w.var(v.GetName()).setVal(v.getVal())
w.saveSnapshot('background_fit',w.allVars());w.writeToFile(str(d/'postfit.root'))
"""
subprocess.run([str(P/'bin/python'),'-c',code],check=True,stdout=subprocess.DEVNULL);w=str(d/'postfit.root')
run([C,'-M','GoodnessOfFit','-d',w,'--algo','saturated','--fixedSignalStrength','0','-n','_observed']+opts,d/'observed')
seed=75023;g=d/f'toys_{N}';run([C,'-M','GenerateOnly','-d',w,'--saveToys','--toysFrequentist','--bypassFrequentistFit','--expectSignal','0','--genBinnedChannels',','.join(f'{reg}_ch{i}' for reg in ['SR','CR'] for i in range(1,5)),'-t',str(N),'-s',str(seed),'-n','_toys']+opts,g)
tf=str(g/f'higgsCombine_toys.GenerateOnly.mH7.5.{seed}.root')
for fitted,ws in [('mixed',w),('3',str(R/'degree3/workspace.root'))]:
 run([C,'-M','GoodnessOfFit','-d',ws,'--algo','saturated','--fixedSignalStrength','0','--toysFile',tf,'--toysFrequentist','--bypassFrequentistFit','-t',str(N),'-s',str(seed),'-n','_gof']+opts,g/f'fit_degree{fitted}')
