from pathlib import Path
import ROOT,json,sys
ROOT.gSystem.Load('libHiggsAnalysisCombinedLimit');P=Path.cwd();R=P/'runs/20261030_sr7_cr3_scan';d=R/('m'+sys.argv[1])/'joint';f=ROOT.TFile.Open(str(d/'workspace.root'));w=f.Get('w');a=json.loads((R/'anchor/fit.json').read_text());assert a['status']==0 and a['covQual']==3 and a['minimum']>0
for k,v in a['parameters'].items():
 if w.var(k):w.var(k).setVal(v)
w.var('r').setVal(.01);w.writeToFile(str(d/'workspace_seeded.root'));f.Close();(d/'workspace.root').rename(d/'workspace_unseeded.root');(d/'workspace_seeded.root').rename(d/'workspace.root')
