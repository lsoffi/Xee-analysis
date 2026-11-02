from pathlib import Path
import json,subprocess,time,os,shutil,ROOT
ROOT.gSystem.Load('libHiggsAnalysisCombinedLimit');P=Path.cwd();R=P/'runs/20261030_sr6_cr4_scan';S=P/'scripts/sr6_cr4_scan';PY=str(P/'bin/python');PLOT='/Users/imb/Documents/Codex/2026-10-27/referenced-chatgpt-conversation-this-is-an/work/fit-env/bin/python'
for it in range(180):
 for result in sorted(R.glob('m*/results.json')):
  m=result.parent.name[1:];d=result.parent/'joint'
  if (d/'fd_refinement.json').exists():continue
  wd=d/'fd_precision';wd.mkdir(exist_ok=True);cmd=json.loads((d/'fitdiagnostics/command.json').read_text());cmd[cmd.index('--cminDefaultMinimizerStrategy')+1]='2';cmd[cmd.index('--cminDefaultMinimizerTolerance')+1]='0.001';cmd+=['--cminDefaultMinimizerPrecision','1e-12'];(wd/'command.json').write_text(json.dumps(cmd));start=time.time()
  with (wd/'combine.log').open('w') as out:rr=subprocess.run(cmd,cwd=wd,stdout=out,stderr=subprocess.STDOUT,timeout=600)
  f=ROOT.TFile.Open(str(wd/'fitDiagnostics_local.root'));old=ROOT.TFile.Open(str(d/'fitdiagnostics/fitDiagnostics_local.root'));q={};good=True
  for stage in ['b','s']:
   fit=f.Get('fit_'+stage);o=old.Get('fit_'+stage);q[stage]=dict(status=fit.status(),covQual=fit.covQual(),edm=fit.edm(),old_status=o.status(),old_covQual=o.covQual(),r=fit.floatParsFinal().find('r').getVal() if stage=='s' else 0,delta_nll=fit.minNll()-o.minNll());good &= fit.status()==0 and fit.covQual()==3
  f.Close();old.Close();q['promoted']=bool(good);q['seconds']=time.time()-start
  if good:
   (d/'fitdiagnostics').rename(d/'fitdiagnostics_initial');wd.rename(d/'fitdiagnostics')
   for script,interpreter in [('quality.py',PY),('report.py',PLOT)]:
    with (result.parent/('refine_'+script+'.log')).open('w') as out:subprocess.run([interpreter,str(S/script),m],cwd=P,stdout=out,stderr=subprocess.STDOUT,env=dict(os.environ,MPLCONFIGDIR=str(P/'work/mpl-cache')))
  (d/'fd_refinement.json').write_text(json.dumps(q,indent=2));print(m,q,flush=True)
 if len(list(R.glob('m*/joint/fd_refinement.json')))==31:break
 time.sleep(15)
