from pathlib import Path
import ROOT,json,hashlib
ROOT.gSystem.Load('libHiggsAnalysisCombinedLimit');P=Path.cwd();R=P/'runs/20261030_sr_degree_selection';out=[]
for k in range(65,86):
 m=f'{k/10:.1f}';p=P/'data/mass_scan/261028/cards_region2_data_envelope_allCorrections_binned/ee'/m/'Xee_ee_4_allYears.root';f=ROOT.TFile.Open(str(p));w=f.Get('w');w.var('r').setVal(1);out.append(dict(mass=k/10,rates=[w.function(f'n_exp_binch{i}_proc_Zd').getVal() for i in range(1,5)],source=str(p),sha256=hashlib.sha256(p.read_bytes()).hexdigest()));f.Close()
(R/'signal_rates.json').write_text(json.dumps(out,indent=2));print('21 rates exported')
