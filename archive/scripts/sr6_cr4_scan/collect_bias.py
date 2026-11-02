from pathlib import Path
import ROOT,json,csv
ROOT.gSystem.Load('libHiggsAnalysisCombinedLimit');P=Path.cwd();R=P/'runs/20261030_sr6_cr4_scan';out=[]
for receipt in sorted((R/'bias').glob('*/m*/r*/fit.log.receipt.json')):
 d=receipt.parent;src=d/'fitDiagnostics_bias.root';f=ROOT.TFile.Open(str(src));fit=f.Get('fit_s') if f else None
 if not fit or not hasattr(fit,'status'):out.append(dict(family=d.parent.parent.name,mass=float(d.parent.name[1:]),rtrue=int(d.name[1:]),failed=True));continue
 r=fit.floatParsFinal().find('r');rt=int(d.name[1:]);out.append(dict(family=d.parent.parent.name,mass=float(d.parent.name[1:]),rtrue=rt,rhat=r.getVal(),error=r.getError(),bias=r.getVal()-rt,bias_sigma=(r.getVal()-rt)/r.getError() if r.getError()>0 else None,status=fit.status(),covQual=fit.covQual(),edm=fit.edm(),valid=fit.status()==0 and fit.covQual()==3 and r.getError()>0 and abs(r.getVal())<29.99));f.Close()
(R/'bias_summary.json').write_text(json.dumps(out,indent=2));print(len(out),'bias fits,',sum(r.get('valid',False) for r in out),'status0/cov3',flush=True)

with (R/'bias_summary.csv').open('w') as outcsv:
 writer=csv.DictWriter(outcsv,fieldnames=['family','mass','rtrue','rhat','error','bias','bias_sigma','status','covQual','edm','valid'],extrasaction='ignore');writer.writeheader();writer.writerows(out)
