from pathlib import Path
import ROOT,json
ROOT.gSystem.Load('libHiggsAnalysisCombinedLimit');P=Path.cwd();R=P/'runs/20261029_cr_windows';summ=json.loads((R/'control_comparison.json').read_text());sr=json.loads((P/'runs/20261029_upsilon_control_cheb3/signal/audit.json').read_text());eras=['2022','2022EE','2023','2023BPix']
for tag,vals in summ.items():
 d=R/tag/'degree3';f=ROOT.TFile.Open(str(d/'gof/workspace.root'));w=f.Get('w');w.loadSnapshot('control_fit')
 for arg in list(w.allFunctions())+list(w.allVars())+list(w.allPdfs()):
  for attr in list(arg.attributes()):
   if str(attr).startswith('ORIGNAME:'):arg.setAttribute(str(attr),False)
 names=[];gs=[]
 for i,era in enumerate(eras):
  for proc in ['upsilon1s','upsilon2s']:names.append(f'CMS_EXO25020_bkgScale_{proc}_{era}');gs.append(sr['eras'][i]['processes'][proc]['local_rate_r1'])
 target=2;T=ROOT.RooRealVar('SR_upsilon_events','Expected Upsilon events in SR',vals['3']['SR_sum'],0,5000);args=ROOT.RooArgList(T);expr='@0'
 for i,(name,g) in enumerate(zip(names,gs)):
  if i==target:continue
  args.add(w.var(name));expr+=f'-({g!r})*@{args.getSize()-1}'
 scale=ROOT.RooFormulaVar('profile_scale1s_2022EE',f'({expr})/({gs[target]!r})',args);c=ROOT.RooCustomizer(w.obj('ModelConfig').GetPdf(),'profile');c.replaceArg(w.var(names[target]),scale);pdf=c.build();posit=ROOT.RooGenericPdf('profile_positive_scale','(@0>=0)&&(@0<=10)',ROOT.RooArgList(scale));pdf.addExtraConstraint(posit)
 nw=ROOT.RooWorkspace('w');getattr(nw,'import')(pdf,ROOT.RooFit.RecycleConflictNodes());getattr(nw,'import')(w.data('data_obs'));mc=ROOT.RooStats.ModelConfig('ModelConfig',nw);mc.SetPdf(nw.pdf(pdf.GetName()));mc.SetParametersOfInterest(ROOT.RooArgSet(nw.var('SR_upsilon_events')))
 for getter,setter in [('GetObservables','SetObservables'),('GetNuisanceParameters','SetNuisanceParameters'),('GetGlobalObservables','SetGlobalObservables')]:
  a=ROOT.RooArgSet()
  for v in getattr(w.obj('ModelConfig'),getter)():
   if nw.arg(v.GetName()):a.add(nw.arg(v.GetName()))
  getattr(mc,setter)(a)
 getattr(nw,'import')(mc);nw.saveSnapshot('profile_start',nw.allVars());nw.writeToFile(str(d/'tail_profile.root'));f.Close()
