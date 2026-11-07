from pathlib import Path
import ROOT,json,subprocess,shutil,os
ROOT.gSystem.Load('libHiggsAnalysisCombinedLimit');P=Path.cwd();R=P/'runs/20261030_sr7_cr3_scan';S=P/'scripts/sr7_cr3_scan';PLOT='/Users/imb/Documents/Codex/2026-10-27/referenced-chatgpt-conversation-this-is-an/work/fit-env/bin/python'
for result in sorted(R.glob('m*/results.json')):
 a=json.loads(result.read_text());d=result.parent/'joint'
 if a['significance_consistent'] or not a['fit_ok'] or (d/'fd_minima.json').exists():continue
 wf=ROOT.TFile.Open(str(d/'workspace.root'));w=wf.Get('w');mc=w.obj('ModelConfig');nll=mc.GetPdf().createNLL(w.data('data_obs'),ROOT.RooFit.Extended(True),ROOT.RooFit.Constrain(mc.GetNuisanceParameters()),ROOT.RooFit.GlobalObservables(mc.GetGlobalObservables()),ROOT.RooFit.Offset(False));baseline={v.GetName():v.getVal() for v in w.allVars()};rows=[]
 for label,prefit in [('fitdiagnostics',None),('fd_negative_small',-.1),('fd_negative_large',-1.)]:
  dd=d/label
  if prefit is not None:
   dd.mkdir(exist_ok=True);cmd=json.loads((d/'fitdiagnostics/command.json').read_text());cmd[cmd.index('--preFitValue')+1]=str(prefit);(dd/'command.json').write_text(json.dumps(cmd,indent=2))
   with (dd/'combine.log').open('w') as out:subprocess.run(cmd,cwd=dd,stdout=out,stderr=subprocess.STDOUT,timeout=120)
  f=ROOT.TFile.Open(str(dd/'fitDiagnostics_local.root'));fit=f.Get('fit_s')
  if not fit:continue
  for name,v in baseline.items():w.var(name).setVal(v)
  for p in fit.floatParsFinal():w.var(p.GetName()).setVal(p.getVal())
  rows.append(dict(folder=label,status=fit.status(),covQual=fit.covQual(),nll=nll.getVal(),r=w.var('r').getVal()));f.Close()
 valid=[r for r in rows if r['status']==0 and r['covQual']==3];best=min(valid,key=lambda r:r['nll']);(d/'fd_minima.json').write_text(json.dumps(dict(attempts=rows,selected=best),indent=2));print(a['mass'],rows,'BEST',best,flush=True)
 if best['folder']!='fitdiagnostics':shutil.move(d/'fitdiagnostics',d/'fitdiagnostics_initial');shutil.copytree(d/best['folder'],d/'fitdiagnostics')
 for py,script,log,args in [(str(P/'bin/python'),'quality.py','quality_minima.log',[]),(PLOT,'report.py','report_minima.log',[])]:
  with (result.parent/log).open('w') as out:subprocess.run([py,str(S/script),f"{a['mass']:.1f}"]+args,cwd=P,stdout=out,stderr=subprocess.STDOUT,env=dict(os.environ,MPLCONFIGDIR=str(P/'work/mpl-cache')))
 wf.Close()
