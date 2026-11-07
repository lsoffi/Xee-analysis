import ROOT
from pathlib import Path
ROOT.gROOT.SetBatch(True)
ROOT.gSystem.Load('libHiggsAnalysisCombinedLimit')
f=ROOT.TFile.Open('/Users/imb/Downloads/Xee_ee_4_allYears.root');w=f.Get('w');x=w.var('mass');cat=ROOT.RooCategory('CMS_channel','CMS_channel');cat.defineType('ch1',0)
cheb=w.pdf('shapeBkg_dy_ch1').getPdf(1)
# Preserve the original Chebyshev basis and full-range yield definitions.
x.setRange('originalFullRange',6,12)
x.setRange('fitRange',6,11)
cheb.selectNormalizationRange('originalFullRange',True)
ranged_norms={};keep=[]
for proc,pdf in [('Zd',w.pdf('shapeSig_Zd_ch1')),('dy',cheb),('upsilon1s',w.pdf('shapeBkg_upsilon1s_ch1')),('upsilon2s',w.pdf('shapeBkg_upsilon2s_ch1'))]:
 num=pdf.createIntegral(ROOT.RooArgSet(x),ROOT.RooFit.Range('fitRange'))
 den=pdf.createIntegral(ROOT.RooArgSet(x),ROOT.RooFit.Range('originalFullRange'))
 orig=w.function('n_exp_binch1_proc_'+proc)
 norm=ROOT.RooFormulaVar('n_exp_binch1_proc_'+proc+'_ranged','@0*@1/@2',ROOT.RooArgList(orig,num,den))
 for attr in ['combine.channel','combine.process']:
  val=orig.getStringAttribute(attr)
  if val:norm.setStringAttribute(attr,val)
 norm.setAttribute('combine.signal',proc=='Zd')
 ranged_norms[proc]=norm;keep.extend([num,den,norm])
 print('Range yield',proc,norm.getVal(),'full yield',orig.getVal())
x.setRange(6,11);x.setBins(375)
ps=ROOT.RooArgList(cheb,w.pdf('shapeBkg_upsilon1s_ch1'),w.pdf('shapeBkg_upsilon2s_ch1'))
ns=ROOT.RooArgList(*[ranged_norms[p] for p in ['dy','upsilon1s','upsilon2s']])
b=ROOT.RooAddPdf('pdf_2022_bonly_nuis','',ps,ns)
sps=ROOT.RooArgList(w.pdf('shapeSig_Zd_ch1'));sps.add(ps)
sns=ROOT.RooArgList(ranged_norms['Zd']);sns.add(ns)
s=ROOT.RooAddPdf('pdf_2022_nuis','',sps,sns)
keys=['CMS_eff_e_id','CMS_eff_e_reco','CMS_eff_trigger_e','lumi_1','CMS_scale_e','CMS_EXO25020_signalModelSigmaNuisance_2022']
cl=ROOT.RooArgList(*[w.pdf(k+'_Pdf') for k in keys]);cl.add(s)
prod=ROOT.RooProdPdf('pdf_2022','',cl)
model=ROOT.RooSimultaneous('model_s','model_s',cat);model.addPdf(prod,'ch1')
clb=ROOT.RooArgList(*[w.pdf(k+'_Pdf') for k in keys]);clb.add(b)
prodb=ROOT.RooProdPdf('pdf_2022_bonly','',clb)
modelb=ROOT.RooSimultaneous('model_b','model_b',cat);modelb.addPdf(prodb,'ch1')
ws=ROOT.RooWorkspace('w');imp=getattr(ws,'import');imp(model);imp(modelb,ROOT.RooFit.RecycleConflictNodes())
weight=ROOT.RooRealVar('_weight_','_weight_',1,0,1e9)
data=ROOT.RooDataSet('data_obs','data_obs',ROOT.RooArgSet(x,cat,weight),ROOT.RooFit.WeightVar(weight))
original=w.data('data_obs')
for i in range(original.numEntries()):
 row=original.get(i)
 if row.getCatLabel('CMS_channel')!='ch1' or not 6 <= row.getRealValue('mass') < 11:continue
 x.setVal(row.getRealValue('mass'));cat.setLabel('ch1');data.add(ROOT.RooArgSet(x,cat,weight),original.weight())
imp(data)
ws.var('r').setRange(-30,30)
ob=ROOT.RooArgSet(ws.var('mass'),ws.cat('CMS_channel'));params=ws.pdf('model_s').getParameters(ob);nu=ROOT.RooArgSet()
for p in params:
 if not p.isConstant() and p.GetName()!='r':nu.add(p)
gl=ROOT.RooArgSet(*[ws.var(k+'_In') for k in keys])
mc=ROOT.RooStats.ModelConfig('ModelConfig',ws);mc.SetPdf(ws.pdf('model_s'));mc.SetObservables(ob);mc.SetParametersOfInterest(ROOT.RooArgSet(ws.var('r')));mc.SetNuisanceParameters(nu);mc.SetGlobalObservables(gl);imp(mc)
mcb=ROOT.RooStats.ModelConfig('ModelConfig_bonly',ws);mcb.SetPdf(ws.pdf('model_b'));mcb.SetObservables(ob);mcb.SetParametersOfInterest(ROOT.RooArgSet(ws.var('r')));mcb.SetNuisanceParameters(nu);mcb.SetGlobalObservables(gl);imp(mcb)
ws.writeToFile('outputs/Xee_ee_4_2022_chebyshev_to11.root')
print('WROTE 2022-only Chebyshev workspace, entries',data.numEntries(),'sum',data.sumEntries())
