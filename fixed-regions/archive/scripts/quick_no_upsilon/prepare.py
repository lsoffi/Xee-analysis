from pathlib import Path
import ROOT,json,hashlib
ROOT.gSystem.Load('libHiggsAnalysisCombinedLimit');P=Path.cwd();R=P/'runs/20261030_quick_no_upsilon';cases=[]
for mass in ['7.0','7.5','8.0']:
 old=P/'runs/20261029_cr91_dcb_scan'/('m'+mass);d=R/('m'+mass);d.mkdir(exist_ok=True);src=old/'signal/workspace.root';f=ROOT.TFile.Open(str(src));w=f.Get('w');q=json.loads((old/'joint/quality.json').read_text())
 for p in q['b']['parameters']:
  v=w.var(p['name'])
  if v:v.setVal(p['value'])
 w.var('r').setVal(0);w.saveSnapshot('truth',w.allVars());w.writeToFile(str(d/'generator.root'));ys=[]
 for v in w.allVars():
  if v.GetName().startswith('CMS_EXO25020_bkgScale_upsilon'):
   ys.append(dict(name=v.GetName(),truth=v.getVal()));v.setVal(0);v.setConstant(True)
 w.writeToFile(str(d/'fit_noY.root'));a=json.loads((old/'results.json').read_text());case=dict(mass=mass,range=a['range'],injection=a['limits']['expected_50'],injection_definition='Median expected 95% CLs limit of previous SR+CR DCB Cheb3 baseline',Y_scales=ys,source=str(src),source_sha256=hashlib.sha256(src.read_bytes()).hexdigest(),toys_per_injection=30);(d/'audit.json').write_text(json.dumps(case,indent=2));cases.append(case);f.Close()
(R/'cases.json').write_text(json.dumps(cases,indent=2))
