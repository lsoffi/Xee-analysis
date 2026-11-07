from pathlib import Path
import ROOT,json
ROOT.gSystem.Load('libHiggsAnalysisCombinedLimit');P=Path.cwd();R=P/'runs/20261029_gof_ftest_75'
for n,folder in [(2,'20261029_upsilon_control'),(3,'20261029_upsilon_control_cheb3')]:
 d=R/f'degree{n}';d.mkdir(exist_ok=True);f=ROOT.TFile.Open(str(P/'runs'/folder/'joint/workspace.root'));w=f.Get('w');q=json.loads((P/'runs'/folder/'parameters.json').read_text())['joint']['b']
 for p in q['parameters']:
  if w.var(p['name']):w.var(p['name']).setVal(p['value'])
 w.var('r').setVal(0);w.saveSnapshot('background_fit',w.allVars());w.writeToFile(str(d/'workspace.root'));f.Close()
