from pathlib import Path
import ROOT,json,hashlib
ROOT.gROOT.SetBatch(True);ROOT.gSystem.Load('libHiggsAnalysisCombinedLimit');P=Path.cwd();R=P/'runs/20261029_upsilon_peak_windows';fractions=json.loads((R/'signal_fractions.json').read_text());common=json.loads((R/'common_windows.json').read_text());cores=json.loads((R/'core_windows.json').read_text());cand={'current':[[9.093333333333334,10.506666666666668]],**{str(round(100*r['target_fraction'])):r['union'] for r in common},**{r['label']:r['union'] for r in cores}};expected={(r['mass'],r['era'],r['candidate']):r['signal_fraction_in_CR'] for r in fractions};eras=['2022','2022EE','2023','2023BPix'];maxdiff=0;checks=0
base=P/'data/mass_scan/261028/cards_region2_data_envelope_allCorrections_binned/ee'
for mass in sorted({r['mass'] for r in fractions}):
 f=ROOT.TFile.Open(str(base/f'{mass:.1f}'/'Xee_ee_4_allYears.root'));w=f.Get('w');x=w.var('mass');obs=ROOT.RooArgSet(x)
 for i,era in enumerate(eras,1):
  pdf=w.pdf(f'shapeSig_Zd_ch{i}')
  for name,regions in cand.items():
   total=0
   for lo,hi in regions:x.setRange('candidate_check',lo,hi);integ=pdf.createIntegral(obs,ROOT.RooFit.NormSet(obs),ROOT.RooFit.Range('candidate_check'));total+=integ.getVal()
   diff=abs(total-expected[(mass,era,name)]);maxdiff=max(maxdiff,diff);checks+=1;assert diff<1e-8,(mass,era,name,diff)
 f.Close()
prov=json.loads((R/'provenance.json').read_text());assert all(hashlib.sha256(Path(p).read_bytes()).hexdigest()==h for p,h in prov['hashes'].items())
v=json.loads((R/'validation.json').read_text());v.update(root_integral_checks=checks,max_root_integral_difference=maxdiff,source_hashes_unchanged=True);(R/'validation.json').write_text(json.dumps(v,indent=2));print(v)
