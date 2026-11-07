from pathlib import Path
import ROOT,json,subprocess,time,shutil
ROOT.gSystem.Load('libHiggsAnalysisCombinedLimit');P=Path.cwd();R=P/'runs/20261030_sr6_cr4_scan'
for it in range(180):
 for ref in sorted((R/'bias').glob('*/m*/r*/refinement.json')):
  d=ref.parent
  if (d/'seeded_refinement.json').exists():continue
  f=ROOT.TFile.Open(str(d/'fitDiagnostics_bias.root'));fit=f.Get('fit_s');rh=fit.floatParsFinal().find('r');oldr=rh.getVal();olde=rh.getError()
  if fit.status()==0 and fit.covQual()==3:f.Close();(d/'seeded_refinement.json').write_text(json.dumps(dict(needed=False)));continue
  m=d.parent.name[1:];wd=d/'seeded_precision';wd.mkdir(exist_ok=True);src=ROOT.TFile.Open(str(R/('m'+m)/'joint/workspace.root'));w=src.Get('w')
  for v in fit.floatParsFinal():
   z=w.var(v.GetName())
   if z:
    z.setVal(v.getVal())
    if v.getError()>0:z.setError(v.getError())
  for era in ['2022','2022EE','2023','2023BPix']:
   bv=fit.floatParsFinal().find('local_b1_'+era);cv=fit.floatParsFinal().find('local_c1_'+era);idx=1 if bv and bv.getError()>0 and (not cv or cv.getError()==0) else 0;w.cat('pdfindex_'+era).setIndex(idx)
  w.writeToFile(str(wd/'workspace.root'));src.Close();f.Close();cmd=json.loads((d/'fit.log.command.json').read_text());cmd[cmd.index('-d')+1]=str(wd/'workspace.root');cmd[cmd.index('--preFitValue')+1]=str(oldr);cmd[cmd.index('--cminDefaultMinimizerStrategy')+1]='2';cmd+=['--cminDefaultMinimizerPrecision','1e-12'];data_cmd=json.loads((R/('m'+m)/'joint/fitdiagnostics/command.json').read_text());rmax=max(3,float(data_cmd[data_cmd.index('--rMax')+1]));cmd[cmd.index('--rMin')+1]=str(-rmax);cmd[cmd.index('--rMax')+1]=str(rmax);(wd/'command.json').write_text(json.dumps(cmd));t=time.time()
  with (wd/'combine.log').open('w') as out:res=subprocess.run(cmd,cwd=wd,stdout=out,stderr=subprocess.STDOUT,timeout=600)
  f=ROOT.TFile.Open(str(wd/'fitDiagnostics_bias.root'));fit=f.Get('fit_s');r=fit.floatParsFinal().find('r');good=fit.status()==0 and fit.covQual()==3 and abs(r.getVal()-oldr)<.02*max(olde,r.getError(),.001) and abs(r.getVal())<.95*rmax;q=dict(needed=True,status=fit.status(),covQual=fit.covQual(),edm=fit.edm(),r=r.getVal(),old_r=oldr,error=r.getError(),old_error=olde,promoted=good,rmax=rmax,seconds=time.time()-t);f.Close()
  if good:
   (d/'fitDiagnostics_bias.root').rename(d/'fitDiagnostics_before_seed.root');shutil.copy2(wd/'fitDiagnostics_bias.root',d/'fitDiagnostics_bias.root')
  (d/'seeded_refinement.json').write_text(json.dumps(q,indent=2));print(str(d.relative_to(R)),q,flush=True)
 if len(list((R/'bias').glob('*/m*/r*/seeded_refinement.json')))==124:break
 time.sleep(20)
