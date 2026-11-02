from pathlib import Path
import ROOT,json,math
ROOT.gSystem.Load('libHiggsAnalysisCombinedLimit');R=Path('runs/20261029_gof_ftest_75');out={}
for key,filename in [('degree2','workspace.root'),('degree3','workspace.root'),('degree2SR_3CR','postfit.root')]:
 f=ROOT.TFile.Open(str(R/key/filename));w=f.Get('w');w.loadSnapshot('background_fit');data=w.data('data_obs');out[key]={}
 for region,varname in [('SR','mass'),('CR','mass_CR')]:
  vals=[]
  for i in range(1,5):
   ch=f'{region}_ch{i}';p=w.obj('ModelConfig').GetPdf().getPdf(ch);pdf=next(x for x in p.getComponents() if x.InheritsFrom('RooAddPdf'));x=w.var(varname);obs=ROOT.RooArgSet(x);bw=(x.getMax()-x.getMin())/x.getBins();total=pdf.expectedEvents(obs);nobs=0;term=0
   d=data.reduce(ROOT.RooFit.Cut(f'CMS_channel==CMS_channel::{ch}'))
   for j in range(d.numEntries()):
    row=d.get(j);x.setVal(row.getRealValue(varname));n=d.weight();mu=total*pdf.getVal(obs)*bw;nobs+=n;term+=n*math.log(n/mu) if n else 0
   vals.append(2*(total-nobs+term))
  out[key][region]=dict(total=sum(vals),eras=vals)
 f.Close()
(R/'regional_exact.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2))
