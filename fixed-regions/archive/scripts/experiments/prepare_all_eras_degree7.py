import ROOT,json
from pathlib import Path
ROOT.gROOT.SetBatch(True);ROOT.gSystem.Load('libHiggsAnalysisCombinedLimit')
ERAS=[('2022','ch1'),('2022EE','ch2'),('2023','ch3'),('2023BPix','ch4')]
f=ROOT.TFile.Open('/Users/imb/Downloads/Xee_ee_4_allYears.root')
for selection in [[e] for e in ERAS]+[ERAS]:
 w=f.Get('w').Clone();x=w.var('mass');cat=ROOT.RooCategory('CMS_channel','CMS_channel');ws=ROOT.RooWorkspace('w');imp=getattr(ws,'import');model=ROOT.RooSimultaneous('model_s','',cat);modelb=ROOT.RooSimultaneous('model_b','',cat);keep=[];allkeys=set()
 seeds=[w.var(f't{i}_2022_envelope').getVal() for i in range(7)]
 for era,ch in selection:
  cat.defineType(ch,int(ch[2:])-1)
  coeff=ROOT.RooArgList()
  for i in range(7):
   v=w.var(f't{i}_{era}_envelope');v.setVal(seeds[i]);coeff.add(v)
  cheb=ROOT.RooChebychev(f'dy_degree7_{era}','',x,coeff)
  ps=ROOT.RooArgList(cheb,w.pdf('shapeBkg_upsilon1s_'+ch),w.pdf('shapeBkg_upsilon2s_'+ch));ns=ROOT.RooArgList(*[w.function('n_exp_bin'+ch+'_proc_'+p) for p in ['dy','upsilon1s','upsilon2s']]);bp=ROOT.RooAddPdf('pdf_'+era+'_bonly_nuis','',ps,ns)
  sp=ROOT.RooArgList(w.pdf('shapeSig_Zd_'+ch));sp.add(ps);sn=ROOT.RooArgList(w.function('n_exp_bin'+ch+'_proc_Zd'));sn.add(ns);s=ROOT.RooAddPdf('pdf_'+era+'_nuis','',sp,sn)
  keys=['CMS_eff_e_id','CMS_eff_e_reco','CMS_eff_trigger_e','lumi_1','CMS_scale_e','CMS_EXO25020_signalModelSigmaNuisance_'+era]
  if era.startswith('2023'):keys+=['lumi_2']
  allkeys.update(keys);cl=ROOT.RooArgList(*[w.pdf(k+'_Pdf') for k in keys]);cl.add(s);prod=ROOT.RooProdPdf('pdf_'+era,'',cl);clb=ROOT.RooArgList(*[w.pdf(k+'_Pdf') for k in keys]);clb.add(bp);prodb=ROOT.RooProdPdf('pdf_'+era+'_bonly','',clb);model.addPdf(prod,ch);modelb.addPdf(prodb,ch);keep += [cheb,ps,ns,bp,sp,sn,s,cl,prod,clb,prodb]
 imp(model);imp(modelb,ROOT.RooFit.RecycleConflictNodes())
 weight=ROOT.RooRealVar('_weight_','_weight_',1,0,1e9);data=ROOT.RooDataSet('data_obs','',ROOT.RooArgSet(x,cat,weight),ROOT.RooFit.WeightVar(weight));orig=w.data('data_obs');chs={ch for _,ch in selection}
 for i in range(orig.numEntries()):
  row=orig.get(i);label=row.getCatLabel('CMS_channel')
  if label not in chs:continue
  x.setVal(row.getRealValue('mass'));cat.setLabel(label);data.add(ROOT.RooArgSet(x,cat,weight),orig.weight())
 imp(data);ws.var('r').setRange(-30,30);obs=ROOT.RooArgSet(ws.var('mass'),ws.cat('CMS_channel'));nuis=ROOT.RooArgSet();params=ws.pdf('model_s').getParameters(obs)
 for a in params:
  if not a.isConstant() and a.GetName()!='r':nuis.add(a)
 glob=ROOT.RooArgSet(*[ws.var(k+'_In') for k in sorted(allkeys)])
 for name,pdf in [('ModelConfig','model_s'),('ModelConfig_bonly','model_b')]:
  mc=ROOT.RooStats.ModelConfig(name,ws);mc.SetPdf(ws.pdf(pdf));mc.SetObservables(obs);mc.SetParametersOfInterest(ROOT.RooArgSet(ws.var('r')));mc.SetNuisanceParameters(nuis);mc.SetGlobalObservables(glob);imp(mc)
 label=selection[0][0] if len(selection)==1 else 'combined';ws.writeToFile(f'outputs/Xee_degree7_native_{label}.root');print('PREPARED',label,data.numEntries(),data.sumEntries())
