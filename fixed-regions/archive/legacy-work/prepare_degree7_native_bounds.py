import ROOT,json
from pathlib import Path
ROOT.gSystem.Load('libHiggsAnalysisCombinedLimit');ROOT.gROOT.SetBatch(True)
f=ROOT.TFile.Open('outputs/Xee_ee_4_2022_chebyshev.root');w=f.Get('w');bounds=[]
for n in [0,1,2]:
 for i in range(7):
  v=w.var('t%d_2022_envelope'%i);c=v.getVal();s=v.getError();lo,hi=(-5,5) if n==0 else (c-n*s,c+n*s);v.setRange(lo,hi);bounds.append(dict(nsigma=n,parameter=v.GetName(),prefit=c,sigma=s,low=lo,high=hi))
 tag='degree7_native_'+('free' if n==0 else 'limit'+str(n));w.writeToFile('outputs/Xee_ee_4_2022_chebyshev_'+tag+'.root')
Path('outputs/degree7_native_bounds.json').write_text(json.dumps(bounds,indent=2))
print('450 entries:',w.data('data_obs').numEntries(),'events:',w.data('data_obs').sumEntries())
