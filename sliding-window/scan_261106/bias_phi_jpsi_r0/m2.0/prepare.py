from pathlib import Path
import json,ROOT
ROOT.gSystem.Load('libHiggsAnalysisCombinedLimit')
p=Path.cwd()
folders={'C2': '/afs/cern.ch/work/s/soffi/DIELE-WORK/Analisi-statistica-Xee/scan_261106/middle_m2_deg2/low_cheb2_totalC0/results/m2.0', 'C3': '/afs/cern.ch/work/s/soffi/DIELE-WORK/Analisi-statistica-Xee/scan_261106/middle_m2_deg3/low_cheb3_totalC0/results/m2.0', 'L2': '/afs/cern.ch/work/s/soffi/DIELE-WORK/Analisi-statistica-Xee/scan_261106/middle_m2_deg2/low_logpoly2/results/m2.0', 'L3': '/afs/cern.ch/work/s/soffi/DIELE-WORK/Analisi-statistica-Xee/scan_261106/middle_m2_deg3/low_logpoly3_totalC0/results/m2.0'}
jobs=[[],[]]
for j,(name,folder) in enumerate(folders.items()):
 d=Path(folder)
 wf=ROOT.TFile.Open(str(d/'workspace.root'));w=wf.Get('w')
 ff=ROOT.TFile.Open(str(d/'fitDiagnosticsTest.root'));fit=ff.Get('fit_b')
 assert fit.status()==0 and fit.covQual()==3
 for par in fit.floatParsFinal():w.var(par.GetName()).setVal(par.getVal())
 w.var('r').setVal(0);w.saveSnapshot('background',w.allVars());w.writeToFile(str(p/(name+'_2.0.root')))
 (p/(name+'.json')).write_bytes((d/'model.json').read_bytes())
 ff.Close();wf.Close()
 for t in range(200):jobs[int(t>0)].append(f'2.0 {name} 0 {1450000+j*1000+t}')
for n,rows in zip(['pilot','remaining'],jobs):(p/(n+'.txt')).write_text('\n'.join(rows)+'\n')
