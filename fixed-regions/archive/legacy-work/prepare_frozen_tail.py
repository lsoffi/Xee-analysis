import ROOT,json
from pathlib import Path
ROOT.gSystem.Load('libHiggsAnalysisCombinedLimit');ROOT.gROOT.SetBatch(True)
chosen={'2022':'t4_2022_envelope','2022EE':'t3_2022EE_envelope','2023':'t3_2023_envelope','2023BPix':'t3_2023BPix_envelope'};info=[]
for era in list(chosen)+['combined']:
 f=ROOT.TFile.Open(f'outputs/Xee_degree5_native_{era}.root');w=f.Get('w')
 for e in (list(chosen) if era=='combined' else [era]):
  v=w.var(chosen[e]);v.setConstant(True);info.append(dict(fit=era,era=e,parameter=v.GetName(),frozen_value=v.getVal(),reference='initial degree5 workspace value; numerical seed, not independent background measurement'))
 w.writeToFile(f'outputs/Xee_degree5_native_frozenTail_{era}.root')
Path('outputs/degree5_frozen_tail_choices.json').write_text(json.dumps(info,indent=2));print(json.dumps(info,indent=2))
