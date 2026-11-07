from pathlib import Path
import json,subprocess,time,shutil,ROOT
ROOT.gSystem.Load('libHiggsAnalysisCombinedLimit');P=Path.cwd();R=P/'runs/20261030_sr6_cr4_scan'
for it in range(240):
 for receipt in sorted((R/'bias').glob('*/m*/r*/fit.log.receipt.json')):
  d=receipt.parent
  if (d/'refinement.json').exists():continue
  wd=d/'precision';wd.mkdir(exist_ok=True);cmd=json.loads((d/'fit.log.command.json').read_text());cmd[cmd.index('--cminDefaultMinimizerStrategy')+1]='2';cmd+=['--cminDefaultMinimizerPrecision','1e-12'];(wd/'command.json').write_text(json.dumps(cmd));start=time.time()
  with (wd/'combine.log').open('w') as out:rr=subprocess.run(cmd,cwd=wd,stdout=out,stderr=subprocess.STDOUT,timeout=600)
  f=ROOT.TFile.Open(str(wd/'fitDiagnostics_bias.root'));old=ROOT.TFile.Open(str(d/'fitDiagnostics_bias.root'));fit=f.Get('fit_s');o=old.Get('fit_s');r=fit.floatParsFinal().find('r');ro=o.floatParsFinal().find('r');q=dict(status=fit.status(),covQual=fit.covQual(),edm=fit.edm(),old_status=o.status(),old_covQual=o.covQual(),r=r.getVal(),old_r=ro.getVal(),delta_nll=fit.minNll()-o.minNll(),seconds=time.time()-start);good=fit.status()==0 and fit.covQual()==3;q['promoted']=good;f.Close();old.Close()
  if good:
   (d/'fitDiagnostics_bias.root').rename(d/'fitDiagnostics_initial.root');shutil.copy2(wd/'fitDiagnostics_bias.root',d/'fitDiagnostics_bias.root')
  (d/'refinement.json').write_text(json.dumps(q,indent=2));print(str(d.relative_to(R)),q,flush=True)
 if len(list((R/'bias').glob('*/m*/r*/refinement.json')))==124:break
 time.sleep(15)
