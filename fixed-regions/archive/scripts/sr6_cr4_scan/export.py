from pathlib import Path
import ROOT,json,hashlib
ROOT.gSystem.Load('libHiggsAnalysisCombinedLimit');P=Path.cwd();R=P/'runs/20261030_sr6_cr4_scan';out=[]
for k in range(60,91):
 mass=k/10;src=P/'data/mass_scan/261028/cards_region2_data_envelope_allCorrections_binned/ee'/f'{mass:.1f}'/'Xee_ee_4_allYears.root';f=ROOT.TFile.Open(str(src));w=f.Get('w');w.var('r').setVal(1);x=w.var('mass');obs=ROOT.RooArgSet(x);x.setRange('SR',6,6+232/75);x.setRange('CR',6+232/75,6+338/75);rows=[]
 for i,era in enumerate(['2022','2022EE','2023','2023BPix'],1):
  pdf=w.pdf(f'shapeSig_Zd_ch{i}');rows.append(dict(era=era,full_rate_r1=w.function(f'n_exp_binch{i}_proc_Zd').getVal(),mean=pdf.x0().getVal(),sigmaL=pdf.sigmaL().getVal(),sigmaR=pdf.sigmaR().getVal(),alphaL=pdf.alphaL().getVal(),alphaR=pdf.alphaR().getVal(),nL=pdf.nL().getVal(),nR=pdf.nR().getVal(),fraction_SR=pdf.createIntegral(obs,ROOT.RooFit.NormSet(obs),ROOT.RooFit.Range('SR')).getVal(),fraction_CR=pdf.createIntegral(obs,ROOT.RooFit.NormSet(obs),ROOT.RooFit.Range('CR')).getVal()))
 out.append(dict(mass=mass,source=str(src),sha256=hashlib.sha256(src.read_bytes()).hexdigest(),eras=rows));f.Close()
(R/'signals.json').write_text(json.dumps(out,indent=2));print('31 mass points exported')
