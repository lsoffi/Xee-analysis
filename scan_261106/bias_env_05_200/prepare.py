from pathlib import Path
import json,ROOT
ROOT.gSystem.Load('libHiggsAnalysisCombinedLimit')
p=Path.cwd()
folders={'C2':'low_cheb2_totalC0','L2':'low_logpoly2','C3':'low_cheb3_totalC0','L3':'low_logpoly3'}
jobs=[[],[]]
for j,(name,folder) in enumerate(folders.items()):
 d=p.parent/folder/'results/m0.5'
 wf=ROOT.TFile.Open(str(d/'workspace.root'));w=wf.Get('w')
 ff=ROOT.TFile.Open(str(d/'fitDiagnosticsTest.root'));fit=ff.Get('fit_b')
 assert fit.status()==0 and fit.covQual()==3
 for par in fit.floatParsFinal():w.var(par.GetName()).setVal(par.getVal())
 w.var('r').setVal(0);w.saveSnapshot('background',w.allVars());w.writeToFile(str(p/(name+'_0.5.root')))
 (p/(name+'.json')).write_bytes((d/'model.json').read_bytes())
 ff.Close();wf.Close()
 for t in range(200):jobs[int(t>0)].append(f'0.5 {name} 0 {730000+j*1000+t}')
for n,rows in zip(['pilot','remaining'],jobs):(p/(n+'.txt')).write_text('\n'.join(rows)+'\n')
