import ROOT,json,shutil
from pathlib import Path
ROOT.gSystem.Load('libHiggsAnalysisCombinedLimit');ROOT.gROOT.SetBatch(True)
f=ROOT.TFile.Open('outputs/Xee_ee_4_2022_chebyshev_rebin8_degree5.root');w=f.Get('w');bounds=[]
for nsigma in [1,2]:
 for i in range(5):
  v=w.var('t%d_2022_envelope'%i);center=v.getVal();sigma=v.getError();lo=center-nsigma*sigma;hi=center+nsigma*sigma
  assert sigma>0
  v.setRange(lo,hi);bounds.append(dict(nsigma=nsigma,parameter=v.GetName(),prefit=center,saved_sigma=sigma,low=lo,high=hi))
 tag='rebin8_degree5_limit'+str(nsigma)
 w.writeToFile('outputs/Xee_ee_4_2022_chebyshev_'+tag+'.root')
 shutil.copyfile('outputs/rebin8_rebin8_degree5_binning.json','outputs/rebin8_'+tag+'_binning.json')
Path('outputs/chebyshev_degree5_bounds.json').write_text(json.dumps(bounds,indent=2))
print(json.dumps(bounds,indent=2))
