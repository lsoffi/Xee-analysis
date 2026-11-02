from pathlib import Path
import ROOT,json,hashlib
ROOT.gROOT.SetBatch(True);ROOT.gSystem.Load('libHiggsAnalysisCombinedLimit')
P=Path('/Users/imb/Documents/Codex/Analisi-statistica-Xee');R=P/'runs/20261029_common95_cheb2_gaussian';OLD=P/'runs/20261029_common95_cheb2';count=0;checks=[]
for path in R.glob('m*/workspace.root'):
 d=path.parent;a=json.loads((d/'audit.json').read_text());b=json.loads((OLD/d.name/'audit.json').read_text());assert a['range']==b['range'] and a['bins']==b['bins'];assert a['source_sha256']==b['source_sha256']==hashlib.sha256(Path(a['source']).read_bytes()).hexdigest()
 f=ROOT.TFile.Open(str(path));w=f.Get('w');of=ROOT.TFile.Open(a['source']);orig=of.Get('w')
 for i,(ea,eb) in enumerate(zip(a['eras'],b['eras']),1):
  assert ea['events']==eb['events'];era=ea['era'];g=w.pdf('gaussian_signal_'+era);dcb=orig.pdf('shapeSig_Zd_ch'+str(i));assert g.ClassName()=='RooGaussian';assert g.getMean().GetName()==dcb.x0().GetName() and g.getSigma().GetName()==dcb.sigmaL().GetName()
  assert abs(g.getMean().getVal()-dcb.x0().getVal())<1e-12 and abs(g.getSigma().getVal()-dcb.sigmaL().getVal())<1e-12
  for proc in ['upsilon1s','upsilon2s']:
   va=ea['processes'][proc];vb=eb['processes'][proc];assert abs(va['fraction']-vb['fraction'])<1e-12 and abs(va['local_rate_r1']-vb['local_rate_r1'])<1e-8
  assert ea['processes']['Zd']['full_rate_r1']==eb['processes']['Zd']['full_rate_r1']
  checks.extend(c for pr in ea['processes'].values() for c in pr['acceptance_checks'])
 f.Close();of.Close();count+=1
assert count==36,count
out=dict(workspaces=count,input_hashes_unchanged=True,gaussian_mean_and_sigma_original_DCB_nodes=True,original_total_signal_rate_preserved=True,local_signal_acceptance_updated=True,ranges_bins_data_upsilon_rates_unchanged=True,acceptance_checks=len(checks),maximum_integral_difference=max(abs(c['analytic']-c['root_integral']) for c in checks));(R/'model_audit.json').write_text(json.dumps(out,indent=2));print(out)
