import ROOT,csv,sys,json
from pathlib import Path
ROOT.gROOT.SetBatch(True);ROOT.gSystem.Load('libHiggsAnalysisCombinedLimit')
tag=sys.argv[1];prefix='combine_2022_'+tag
# Reuse the same parameter extraction, changing only output paths and range labels.
s=Path('work/extract_combine_results.py').read_text().replace('inclusive_2022_chebyshev.root','inclusive_2022_chebyshev_'+tag+'.root').replace('Xee_ee_4_2022_chebyshev.root','Xee_ee_4_2022_chebyshev_'+tag+'.root').replace('combine_2022_pre_post',prefix+'_pre_post')
s=s.replace('450 conteggi originali, somma 33341','57 bin dopo rebin ×8, somma 33341' if not tag.startswith('to11') else '47 bin dopo rebin ×8, somma 33056')
if tag.startswith('to11'):s=s.replace('6–12 GeV','6–11 GeV')
s=s.replace("v=w.var(name);row=","v=w.var(name)\n if not v:continue\n row=")
exec(compile(s,'extract_rebin8_parameters','exec'))
f=ROOT.TFile.Open('outputs/fitDiagnostics_inclusive_2022_chebyshev_'+tag+'.root')
for stage in ['prefit','fit_b','fit_s']:
 dr=f.Get('shapes_'+stage+'/ch1');h=dr.Get('total');data=dr.Get('data');countsum=0
 with open('outputs/'+prefix+'_'+stage+'_shapes.csv','w') as o:
  wr=csv.writer(o);wr.writerow(['mass','low','high','width','data','total','signal','dy','upsilon1s','upsilon2s'])
  for i in range(1,h.GetNbinsX()+1):
   width=h.GetBinWidth(i);count=data.GetPointY(i-1)*width;countsum+=count
   wr.writerow([h.GetBinCenter(i),h.GetBinLowEdge(i),h.GetBinLowEdge(i)+width,width,count]+[dr.Get(k).GetBinContent(i)*width for k in ['total','Zd','dy','upsilon1s','upsilon2s']])
 print('EXPORT',tag,stage,'bins',h.GetNbinsX(),'data sum',countsum)
