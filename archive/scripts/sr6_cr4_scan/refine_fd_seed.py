from pathlib import Path
import ROOT,json,subprocess,time,shutil,os
ROOT.gSystem.Load('libHiggsAnalysisCombinedLimit');P=Path.cwd();R=P/'runs/20261030_sr6_cr4_scan';S=P/'scripts/sr6_cr4_scan';PLOT='/Users/imb/Documents/Codex/2026-10-27/referenced-chatgpt-conversation-this-is-an/work/fit-env/bin/python'
for it in range(120):
 for ref in sorted(R.glob('m*/joint/fd_refinement.json')):
  d=ref.parent
  if (d/'fd_seed_refinement.json').exists():continue
  a=json.loads((d.parent/'results.json').read_text());m=d.parent.name[1:]
  if a['fit_ok']:(d/'fd_seed_refinement.json').write_text(json.dumps(dict(needed=False)));continue
  f=ROOT.TFile.Open(str(d/'fitdiagnostics/fitDiagnostics_local.root'));fit=f.Get('fit_s');rh=fit.floatParsFinal().find('r');oldr=rh.getVal();olde=rh.getError();src=ROOT.TFile.Open(str(d/'workspace.root'));w=src.Get('w');wd=d/'fd_seed';wd.mkdir(exist_ok=True)
  for v in fit.floatParsFinal():
   z=w.var(v.GetName())
   if z:
    z.setVal(v.getVal())
    if v.getError()>0:z.setError(v.getError())
  for e in ['2022','2022EE','2023','2023BPix']:
   b=fit.floatParsFinal().find('local_b1_'+e);c=fit.floatParsFinal().find('local_c1_'+e);w.cat('pdfindex_'+e).setIndex(1 if b and b.getError()>0 and (not c or c.getError()==0) else 0)
  w.writeToFile(str(wd/'workspace.root'));src.Close();f.Close();cmd=json.loads((d/'fitdiagnostics/command.json').read_text());cmd[cmd.index('-d')+1]=str(wd/'workspace.root');cmd[cmd.index('--preFitValue')+1]=str(oldr);cmd[cmd.index('--cminDefaultMinimizerStrategy')+1]='2'
  if '--cminDefaultMinimizerPrecision' not in cmd:cmd+=['--cminDefaultMinimizerPrecision','1e-12']
  (wd/'command.json').write_text(json.dumps(cmd));start=time.time()
  with (wd/'combine.log').open('w') as out:res=subprocess.run(cmd,cwd=wd,stdout=out,stderr=subprocess.STDOUT,timeout=600)
  f=ROOT.TFile.Open(str(wd/'fitDiagnostics_local.root'));q={};good=True
  for st in ['b','s']:
   ft=f.Get('fit_'+st);q[st]=dict(status=ft.status(),covQual=ft.covQual(),edm=ft.edm());good &= ft.status()==0 and ft.covQual()==3
  rr=f.Get('fit_s').floatParsFinal().find('r');q.update(needed=True,old_r=oldr,r=rr.getVal(),error=rr.getError(),seconds=time.time()-start);good &= abs(rr.getVal()-oldr)<.02*max(olde,rr.getError(),.001);q['promoted']=bool(good);f.Close()
  if good:
   (d/'fitdiagnostics').rename(d/'fitdiagnostics_before_seed');shutil.copytree(wd,d/'fitdiagnostics')
   for script,py in [('quality.py',str(P/'bin/python')),('report.py',PLOT)]:
    with (d.parent/('seed_'+script+'.log')).open('w') as out:subprocess.run([py,str(S/script),m],cwd=P,stdout=out,stderr=subprocess.STDOUT,env=dict(os.environ,MPLCONFIGDIR=str(P/'work/mpl-cache')))
  (d/'fd_seed_refinement.json').write_text(json.dumps(q,indent=2));print(m,q,flush=True)
 if len(list(R.glob('m*/joint/fd_seed_refinement.json')))==31:break
 time.sleep(20)
