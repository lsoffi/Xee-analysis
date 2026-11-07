from pathlib import Path
import json,ROOT
ROOT.gSystem.Load('libHiggsAnalysisCombinedLimit')
p=Path.cwd(); web=Path('/eos/user/s/soffi/www/Xee/261106/scan_SR_CR')
limits=json.loads((web/'results.json').read_text())
jobs=[]
for i,m in enumerate(['0.4','0.5','0.6']):
 exp=limits['limits'][limits['masses'].index(float(m))][3]
 for j,(model,folder) in enumerate([('cheb','low_cheb3_totalC0'),('log','low_logpoly3_totalC0' if m=='0.6' else 'low_logpoly3')]):
  orig=p.parent/folder/'results'/('m'+m)
  fin=ROOT.TFile.Open(str(orig/'workspace.root'));w=fin.Get('w')
  ff=ROOT.TFile.Open(str(orig/'fitDiagnosticsTest.root'));fit=ff.Get('fit_b')
  assert fit.status()==0 and fit.covQual()==3
  for par in fit.floatParsFinal():w.var(par.GetName()).setVal(par.getVal())
  w.var('r').setVal(0);w.saveSnapshot('background',w.allVars())
  w.writeToFile(str(p/f'{model}_{m}.root'));ff.Close();fin.Close()
  for k,inj in enumerate([0,exp]):
   for t in range(20):jobs.append(f'{m} {model} {inj} {610000+i*1000+j*100+k*20+t}')
(p/'points.txt').write_text('\n'.join(jobs)+'\n')
