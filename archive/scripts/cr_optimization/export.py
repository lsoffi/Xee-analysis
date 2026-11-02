from pathlib import Path
import ROOT,json,hashlib
ROOT.gSystem.Load('libHiggsAnalysisCombinedLimit'); ROOT.gROOT.SetBatch(True)
P=Path.cwd(); src=P/'data/mass_scan/261028/cards_region2_data_envelope_allCorrections_binned/ee/7.5/Xee_ee_4_allYears.root'; f=ROOT.TFile.Open(str(src));w=f.Get('w');x=w.var('mass');obs=ROOT.RooArgSet(x); edges=[6+i/75 for i in range(451)];out={'edges':edges,'source':str(src),'sha256':hashlib.sha256(src.read_bytes()).hexdigest(),'eras':[]}
for i,era in enumerate(['2022','2022EE','2023','2023BPix'],1):
 ch=f'ch{i}'; d=w.data('data_obs').reduce(ROOT.RooFit.Cut(f'CMS_channel==CMS_channel::{ch}')); counts=[0.]*450
 for j in range(d.numEntries()):
  row=d.get(j); k=int((row.getRealValue('mass')-6)*75);counts[k]+=d.weight()
 rec={'era':era,'counts':counts,'templates':[],'fractions':[]}
 for proc in ['upsilon1s','upsilon2s']:
  pdf=w.pdf('shapeBkg_'+proc+'_'+ch);rate=w.function(f'n_exp_bin{ch}_proc_{proc}').getVal(); vals=[]
  for a,b in zip(edges[:-1],edges[1:]):
   x.setRange('exportbin',a,b);vals.append(pdf.createIntegral(obs,ROOT.RooFit.NormSet(obs),ROOT.RooFit.Range('exportbin')).getVal())
  assert abs(sum(vals)-1)<1e-8
  rec['fractions'].append(vals);rec['templates'].append([v*rate for v in vals])
 out['eras'].append(rec)
(P/'runs/20261030_cr_optimization/input.json').write_text(json.dumps(out));print('exported')
