import ROOT,pathlib,json,hashlib,numpy as np
ROOT.gSystem.Load('libHiggsAnalysisCombinedLimit');P=pathlib.Path('/Users/imb/Documents/Codex/Analisi-statistica-Xee');R=P/'runs/20261029_background_validation';rows=[]
for mass in ['7.5','9.5']:
 src=P/f'data/mass_scan/261028/cards_region2_data_envelope_allCorrections_binned/ee/{mass}/Xee_ee_4_allYears.root';fo=ROOT.TFile.Open(str(src));wo=fo.Get('w');wo.var('r').setVal(.1)
 for i in range(1,5):wo.cat(f'CMS_EXO25020_bkgEnvelopeIdx_'+['2022','2022EE','2023','2023BPix'][i-1]).setIndex(1)
 for model in ['original','degree5','cheb7','bernstein']:
  path=R/f'workspaces/m{mass}_{model}{".clean" if model!="original" else ""}.root';f=ROOT.TFile.Open(str(path));w=f.Get('w');w.var('r').setVal(.1)
  obs=w.obj('ModelConfig').GetObservables();params=w.obj('ModelConfig').GetPdf().getParameters(obs);active=[v.GetName() for v in params];checks=[]
  for i in range(1,5):
   for proc in ['Zd','dy','upsilon1s','upsilon2s']:
    name=f'n_exp_binch{i}_proc_{proc}';checks.append(abs(w.function(name).getVal()-wo.function(name).getVal())<1e-9)
   for name in [f'shapeSig_Zd_ch{i}',f'shapeBkg_upsilon1s_ch{i}',f'shapeBkg_upsilon2s_ch{i}']:
    for x in [6.5,7.5,9.5,10.5,11.5]:
     w.var('mass').setVal(x);wo.var('mass').setVal(x);checks.append(abs(w.pdf(name).getVal(ROOT.RooArgSet(w.var('mass')))-wo.pdf(name).getVal(ROOT.RooArgSet(wo.var('mass'))))<1e-9)
  assert all(checks)
  assert w.var('mass').getMin()==6 and w.var('mass').getMax()==12 and w.var('mass').getBins()==450
  for era in ['2022','2022EE','2023','2023BPix']:assert f'CMS_EXO25020_bkgScale_upsilon1s_{era}' in active
  if model=='degree5':assert sum(n.startswith('t') and n.endswith('_envelope') for n in active)==20
  rows.append(dict(mass=mass,model=model,input_sha256=hashlib.sha256(src.read_bytes()).hexdigest(),data_entries=w.data('data_obs').numEntries(),data_sum=w.data('data_obs').sumEntries(),signal_and_yields_and_resonant_shapes_identical=True,active_chebyshev_coefficients=sum(n.startswith('t') and n.endswith('_envelope') for n in active),active_bernstein_coefficients=sum(n.startswith('a') and n.endswith('_envelope') for n in active),discrete_parameters=w.genobj('discreteParams').getSize() if w.genobj('discreteParams') else 0))
(R/'model_integrity_audit.json').write_text(json.dumps(rows,indent=2));print('Verified',len(rows),'model workspaces')
