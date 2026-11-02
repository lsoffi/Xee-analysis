from pathlib import Path
import ROOT,json,math
ROOT.gSystem.Load('libHiggsAnalysisCombinedLimit')
P=Path.cwd();R=P/'runs/20261030_upper_sr_cr';a=json.loads((R/'multistart.json').read_text());primary=json.loads((R/'results.json').read_text())
for q in a:
 d=R/f'degree{q["degree"]}/joint/multistart{q["start"]}'
 f=ROOT.TFile.Open(str(d/'fitDiagnostics_upper.root'));fit=f.Get('fit_b');fw=ROOT.TFile.Open(str(d/'start.root'));w=fw.Get('w')
 for v in fit.floatParsFinal():w.var(v.GetName()).setVal(v.getVal())
 w.var('r').setVal(0);D=0
 for region,obsname in [('SR','mass'),('CR','mass_CR')]:
  for i in range(1,5):
   ch=f'{region}_ch{i}';cp=w.obj('ModelConfig_bonly').GetPdf().getPdf(ch);pdf=next(p for p in cp.getComponents() if p.InheritsFrom('RooAddPdf'));x=w.var(obsname);obs=ROOT.RooArgSet(x);bw=(x.getMax()-x.getMin())/x.getBins();total=pdf.expectedEvents(obs);data=w.data('data_obs').reduce(ROOT.RooFit.Cut(f'CMS_channel==CMS_channel::{ch}'));ns=0;term=0
   for k in range(data.numEntries()):
    row=data.get(k);x.setVal(row.getRealValue(obsname));n=data.weight();mu=total*pdf.getVal(obs)*bw;ns+=n;term+=n*math.log(n/mu) if n else 0
   D+=2*(total-ns+term)
 q['deviance']=D;q['delta_deviance']=D-next(p['D_FD'] for p in primary if p['degree']==q['degree']);q['raw_minNll_difference_not_comparable']=True
 f.Close();fw.Close()
(R/'multistart.json').write_text(json.dumps(a,indent=2));print([(q['degree'],q['start'],q['delta_deviance']) for q in a])
