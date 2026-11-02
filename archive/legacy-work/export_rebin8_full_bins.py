import ROOT,json,csv,sys
ROOT.gSystem.Load('libHiggsAnalysisCombinedLimit');ROOT.gROOT.SetBatch(True)
tag=sys.argv[1];prefix='combine_2022_'+tag
f=ROOT.TFile.Open('outputs/fitDiagnostics_inclusive_2022_chebyshev_'+tag+'.root');wf=ROOT.TFile.Open('outputs/Xee_ee_4_2022_chebyshev_'+tag+'.root');w=wf.Get('w');x=w.var('mass');obs=ROOT.RooArgSet(x);model=w.pdf('pdf_2022_nuis');ps=model.pdfList();ns=model.coefList();initial=w.allVars().snapshot();binning=json.load(open('outputs/rebin8_'+tag+'_binning.json'));edges=binning['edges'];counts=binning['counts']
for stage in ['prefit','fit_b','fit_s']:
 w.allVars().assignValueOnly(initial)
 if stage!='prefit':w.allVars().assignValueOnly(f.Get(stage).floatParsFinal())
 if stage=='fit_b':w.var('r').setVal(0)
 total=0
 with open('outputs/'+prefix+'_'+stage+'_shapes.csv','w') as o:
  wr=csv.writer(o);wr.writerow(['mass','low','high','width','data','total','signal','dy','upsilon1s','upsilon2s'])
  for i,count in enumerate(counts):
   lo,hi=edges[i:i+2];width=hi-lo;x.setVal((lo+hi)/2)
   vals=[ps[j].getVal(obs)*ns[j].getVal()*width for j in range(4)]
   wr.writerow([(lo+hi)/2,lo,hi,width,count,sum(vals)]+vals);total+=sum(vals)
 print(tag,stage,'bins',len(counts),'observed',sum(counts),'expected',total)
assert len(counts)==w.data('data_obs').numEntries()
assert sum(counts)==w.data('data_obs').sumEntries()
