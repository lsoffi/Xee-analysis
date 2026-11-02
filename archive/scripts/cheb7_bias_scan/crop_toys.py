from pathlib import Path
import ROOT,json,sys
ROOT.gSystem.Load('libHiggsAnalysisCombinedLimit');src,dest,workspace=sys.argv[1:4];f=ROOT.TFile.Open(src);wf=ROOT.TFile.Open(workspace);w=wf.Get('w');x=w.var('mass');y=w.var('mass_CR');cat=w.cat('CMS_channel');lo=x.getMin();hi=x.getMax();weight=ROOT.RooRealVar('toy_weight','',0,1e12);obs=ROOT.RooArgSet(x,y,cat,weight);out=ROOT.TFile.Open(dest,'RECREATE');out.mkdir('toys');out.cd('toys');audit=[]
for key in f.Get('toys').GetListOfKeys():
 name=key.GetName();obj=key.ReadObj()
 if name.endswith('_snapshot'):obj.Write(name);continue
 if not obj.InheritsFrom('RooAbsData'):continue
 data=ROOT.RooDataSet(name,'',obs,ROOT.RooFit.WeightVar(weight));count_before={};count_after={}
 for i in range(obj.numEntries()):
  row=obj.get(i);label=row.find('CMS_channel').getLabel();wt=obj.weight();count_before[label]=count_before.get(label,0)+wt
  if label.startswith('SR_') and not lo-1e-9<=row.getRealValue('mass')<hi-1e-9:continue
  cat.setLabel(label);x.setVal(row.getRealValue('mass') if label.startswith('SR_') else (lo+hi)/2);y.setVal(row.getRealValue('mass_CR') if label.startswith('CR_') else (y.getMin()+y.getMax())/2);weight.setVal(wt);data.add(obs,wt);count_after[label]=count_after.get(label,0)+wt
 for label,n in count_before.items():
  if label.startswith('CR_'):assert abs(count_after[label]-n)<1e-7
 data.Write(name);audit.append(dict(toy=name,before=count_before,after=count_after))
out.Close();Path(dest+'.audit.json').write_text(json.dumps(dict(range=[lo,hi],toys=audit),indent=2));print('cropped',len(audit),lo,hi)
