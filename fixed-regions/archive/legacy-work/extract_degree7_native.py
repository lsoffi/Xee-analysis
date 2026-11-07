import ROOT,csv,json
from pathlib import Path
ROOT.gROOT.SetBatch(True);ROOT.gSystem.Load('libHiggsAnalysisCombinedLimit')
for case in ['free','limit1','limit2']:
 tag='degree7_native_'+case;prefix='combine_2022_'+tag
 code=Path('work/extract_combine_results.py').read_text().replace('inclusive_2022_chebyshev.root','inclusive_2022_chebyshev_'+tag+'.root').replace('Xee_ee_4_2022_chebyshev.root','Xee_ee_4_2022_chebyshev_'+tag+'.root').replace('combine_2022_pre_post',prefix+'_pre_post')
 exec(compile(code,'extract_parameters','exec'))
 f=ROOT.TFile.Open('outputs/fitDiagnostics_inclusive_2022_chebyshev_'+tag+'.root')
 for stage in ['prefit','fit_b','fit_s']:
  dr=f.Get('shapes_'+stage+'/ch1');h=dr.Get('total');data=dr.Get('data');countsum=0
  with open('outputs/'+prefix+'_'+stage+'_shapes.csv','w') as o:
   wr=csv.writer(o);wr.writerow(['mass','low','high','width','data','total','signal','dy','upsilon1s','upsilon2s'])
   for i in range(1,h.GetNbinsX()+1):
    width=h.GetBinWidth(i);count=data.GetPointY(i-1)*width;countsum+=count
    wr.writerow([h.GetBinCenter(i),h.GetBinLowEdge(i),h.GetBinLowEdge(i)+width,width,count]+[dr.Get(k).GetBinContent(i)*width for k in ['total','Zd','dy','upsilon1s','upsilon2s']])
  assert h.GetNbinsX()==450
  assert abs(countsum-33341)<1e-6
