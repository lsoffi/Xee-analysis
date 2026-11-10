import ROOT as R,sys,pathlib,shutil,json
R.gSystem.Load('libHiggsAnalysisCombinedLimit')
m,model,dest=sys.argv[1:];p=pathlib.Path(dest)
f=R.TFile.Open('workspace.root');w=f.Get('w');ff=R.TFile.Open('fitDiagnosticsTest.root');fit=ff.Get('fit_b')
assert fit and fit.status()==0 and fit.covQual()==3,'Generator prefit must have status 0 and covQual 3'
for v in fit.floatParsFinal():w.var(v.GetName()).setVal(v.getVal())
w.var('r').setRange(-80,80);w.var('r').setVal(0);w.saveSnapshot('background',w.allVars());w.writeToFile(str(p/f'{model}_{m}.root'))
w.saveSnapshot('injection',w.allVars());w.writeToFile(str(p/f'gen_{model}.root'))
shutil.copy2('model.json',p/f'{model}.json');shutil.copy2('fitDiagnosticsTest.root',p/f'prefit_{model}.root')
