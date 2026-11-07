from pathlib import Path
import json,ROOT
ROOT.gSystem.Load('libHiggsAnalysisCombinedLimit')
p=Path.cwd()
folders={'C2': '/afs/cern.ch/work/s/soffi/DIELE-WORK/Analisi-statistica-Xee/scan_261106/bias_phi_jpsi_r0/C2/results/m1.9', 'C3': '/afs/cern.ch/work/s/soffi/DIELE-WORK/Analisi-statistica-Xee/scan_261106/bias_phi_jpsi_r0/C3/results/m1.9', 'L2': '/afs/cern.ch/work/s/soffi/DIELE-WORK/Analisi-statistica-Xee/scan_261106/bias_phi_jpsi_r0/L2/results/m1.9', 'L3': '/afs/cern.ch/work/s/soffi/DIELE-WORK/Analisi-statistica-Xee/scan_261106/bias_phi_jpsi_r0/L3/results/m1.9'}
jobs=[[],[]]
for j,(name,folder) in enumerate(folders.items()):
 d=Path(folder)
 wf=ROOT.TFile.Open(str(d/'workspace.root'));w=wf.Get('w')
 ff=ROOT.TFile.Open(str(d/'fitDiagnosticsTest.root'));fit=ff.Get('fit_b')
 assert fit.status()==0 and fit.covQual()==3
 for par in fit.floatParsFinal():w.var(par.GetName()).setVal(par.getVal())
 w.var('r').setVal(0);w.saveSnapshot('background',w.allVars());w.writeToFile(str(p/(name+'_1.9.root')))
 (p/(name+'.json')).write_bytes((d/'model.json').read_bytes())
 ff.Close();wf.Close()
 for t in range(200):jobs[int(t>0)].append(f'1.9 {name} 0 {1440000+j*1000+t}')
for n,rows in zip(['pilot','remaining'],jobs):(p/(n+'.txt')).write_text('\n'.join(rows)+'\n')
