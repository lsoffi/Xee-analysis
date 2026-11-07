from pathlib import Path
import json,ROOT
ROOT.gSystem.Load('libHiggsAnalysisCombinedLimit')
p=Path.cwd()
folders={k:str(p.parent/k/'results'/p.name) for k in ('C4','L4')}
jobs=[[],[]]
for j,(name,folder) in enumerate(folders.items()):
 d=Path(folder)
 wf=ROOT.TFile.Open(str(d/'workspace.root'));w=wf.Get('w')
 ff=ROOT.TFile.Open(str(d/'fitDiagnosticsTest.root'));fit=ff.Get('fit_b')
 assert fit.status()==0 and fit.covQual()==3
 for par in fit.floatParsFinal():w.var(par.GetName()).setVal(par.getVal())
 w.var('r').setVal(0);w.saveSnapshot('background',w.allVars());w.writeToFile(str(p/(name+'_2.3.root')))
 (p/(name+'.json')).write_bytes((d/'model.json').read_bytes())
 w.var('r').setRange(-80,80);w.var('r').setVal(1.341796875);w.saveSnapshot('injection',w.allVars());w.writeToFile(str(p/('gen_'+name+'.root')))
 ff.Close();wf.Close()
