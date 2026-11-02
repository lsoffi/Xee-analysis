import ROOT,csv
ROOT.gSystem.Load('libHiggsAnalysisCombinedLimit')
f=ROOT.TFile.Open('outputs/fitDiagnostics_inclusive_2022_chebyshev_to11.root')
for tag in ['prefit','fit_b','fit_s']:
 dr=f.Get('shapes_'+tag+'/ch1');h=dr.Get('total');data=dr.Get('data')
 with open('outputs/combine_2022_to11_'+tag+'_shapes.csv','w') as o:
  wr=csv.writer(o);wr.writerow(['mass','data','total','signal','dy','upsilon1s','upsilon2s'])
  for i in range(1,h.GetNbinsX()+1):
   width=h.GetBinWidth(i);wr.writerow([h.GetBinCenter(i),data.GetPointY(i-1)*width]+[dr.Get(k).GetBinContent(i)*width for k in ['total','Zd','dy','upsilon1s','upsilon2s']])
 print(tag,'sum counts',sum(h.GetBinContent(i)*h.GetBinWidth(i) for i in range(1,h.GetNbinsX()+1)))
