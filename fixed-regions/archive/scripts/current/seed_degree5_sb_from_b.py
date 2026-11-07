import ROOT,json
from pathlib import Path
ROOT.gROOT.SetBatch(True);ROOT.gSystem.Load('libHiggsAnalysisCombinedLimit')
eras=['2022','2022EE','2023','2023BPix'];fits={};files=[];audit=[]
for era in eras:
 f=ROOT.TFile.Open(f'outputs/fitDiagnostics_degree5_stageB_{era}.root');files.append(f);fits[era]=f.Get('fit_b');assert fits[era].status()==0
for case in eras+['combined']:
 f=ROOT.TFile.Open(f'outputs/Xee_degree5_eraStart_{case}.root');w=f.Get('w')
 for era in (eras if case=='combined' else [case]):
  for a in fits[era].floatParsFinal():
   name=a.GetName();target=w.var(name)
   # Use era-specific fitted background coefficients and yields. Shared nuisance starts stay nominal.
   if target and (name.startswith('t') and name.endswith('_'+era+'_envelope') or name.startswith('CMS_EXO25020_bkgScale_') and name.endswith('_'+era)):
    target.setVal(a.getVal());target.setError(a.getError());audit.append(dict(case=case,source_B_era=era,parameter=name,B_value=a.getVal(),SB_initial=target.getVal()))
 w.var('r').setVal(0);w.var('r').setError(.1);w.var('r').setConstant(False)
 for name in [a.GetName() for a in fits[eras[0]].floatParsFinal()]:
  v=w.var(name)
  if v and (name.startswith('CMS_eff_') or name in ['lumi_1','lumi_2','CMS_scale_e']):v.setVal(0)
 w.saveSnapshot('background_initialized',w.allVars());w.writeToFile(f'outputs/Xee_degree5_Binitialized_{case}.root')
Path('outputs/degree5_B_to_SB_initialization_audit.json').write_text(json.dumps(audit,indent=2));print('Seeded all era and combined workspaces from separate B-only fits')
