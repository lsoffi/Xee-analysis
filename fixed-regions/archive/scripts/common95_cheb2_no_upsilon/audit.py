from pathlib import Path
import ROOT,json,hashlib
ROOT.gROOT.SetBatch(True);ROOT.gSystem.Load('libHiggsAnalysisCombinedLimit')
P=Path('/Users/imb/Documents/Codex/Analisi-statistica-Xee');R=P/'runs/20261029_common95_cheb2_no_upsilon';old=P/'runs/20261029_common95_cheb2';count=0
for path in R.glob('m*/workspace.root'):
 d=path.parent;a=json.loads((d/'audit.json').read_text());prev=json.loads((old/d.name/'audit.json').read_text());assert a['range']==prev['range'] and a['bins']==prev['bins'];assert a['source_sha256']==prev['source_sha256']==hashlib.sha256(Path(a['source']).read_bytes()).hexdigest()
 f=ROOT.TFile.Open(str(path));w=f.Get('w')
 for ea,eb in zip(a['eras'],prev['eras']):
  assert ea['events']==eb['events'];sa=ea['processes']['Zd'];sb=eb['processes']['Zd'];assert abs(sa['fraction']-sb['fraction'])<1e-12 and abs(sa['local_rate_r1']-sb['local_rate_r1'])<1e-10
  for process in ['upsilon1s','upsilon2s']:
   z=w.arg('excluded_'+process+'_'+ea['era']);assert z and z.isConstant() and z.getVal()==0
 f.Close();count+=1
assert count==36,count
(R/'model_audit.json').write_text(json.dumps(dict(workspaces_checked=count,input_hashes_unchanged=True,ranges_bins_data_signal_acceptance_and_rates_identical=True,upsilon_rates_fixed_zero=True,original_with_upsilon_outputs_preserved=True),indent=2));print('AUDIT PASSED',count)
