from pathlib import Path
import ROOT,json,sys,hashlib
ROOT.gSystem.Load('libHiggsAnalysisCombinedLimit');P=Path.cwd();R=P/'runs/20261030_cheb7_bias_scan';m=sys.argv[1];d=R/('m'+m);q=json.loads((R/'generator_fit/fit.json').read_text());f=ROOT.TFile.Open(str(d/'joint/workspace.root'));w=f.Get('w')
for name,value in q['parameters'].items():
 if w.var(name):w.var(name).setVal(value)
w.var('r').setVal(0);w.saveSnapshot('generator_truth',w.allVars());w.writeToFile(str(d/'generator.root'));f.Close();sources=[]
for method,folder in [('fixed','20261030_fixed_sr_scan'),('sliding','20261029_cr91_dcb_scan')]:
 src=P/'runs'/folder/('m'+m)/'joint/workspace.root';f=ROOT.TFile.Open(str(src));w=f.Get('w');qp=json.loads((P/'runs'/folder/('m'+m)/'joint/quality.json').read_text())
 for p in qp['b']['parameters']:
  if w.var(p['name']):w.var(p['name']).setVal(p['value'])
 w.var('r').setVal(0);w.writeToFile(str(d/f'fit_{method}.root'));sources.append(dict(method=method,path=str(src),sha256=hashlib.sha256(src.read_bytes()).hexdigest(),range=[w.var('mass').getMin(),w.var('mass').getMax()]));f.Close()
base=json.loads((P/'runs/20261029_cr91_dcb_scan'/('m'+m)/'results.json').read_text());(d/'study.json').write_text(json.dumps(dict(mass=float(m),toys=30,injection=base['limits']['expected_50'],injection_definition='median expected95CLs from previous sliding95 DCB CR baseline; same r for both fits',sources=sources),indent=2))
