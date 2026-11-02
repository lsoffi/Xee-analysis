from pathlib import Path
import ROOT,json,subprocess
ROOT.gSystem.Load('libHiggsAnalysisCombinedLimit');P=Path.cwd();R=P/'runs/20261030_fixed_sr_scan';m='7.2';d=R/('m'+m)/'joint';q=json.loads((d/'quality.json').read_text());rows=[]
f=ROOT.TFile.Open(str(d/'workspace.root'));w=f.Get('w');mc=w.obj('ModelConfig');nll=mc.GetPdf().createNLL(w.data('data_obs'),ROOT.RooFit.Extended(True),ROOT.RooFit.Constrain(mc.GetNuisanceParameters()),ROOT.RooFit.GlobalObservables(mc.GetGlobalObservables()),ROOT.RooFit.Offset(False))
for i,r in enumerate([-.1,0,.1,1]):
 dd=d/f'multistart{i}';dd.mkdir(exist_ok=True);params={p['name']:p['value'] for p in q['b']['parameters']};params['r']=r;cmd=[str(P/'bin/combine'),'-M','MultiDimFit','-d',str(d/'workspace.root'),'-m',m,'--algo','none','--saveFitResult','--saveWorkspace','--rMin','-100','--rMax','100','--cminDefaultMinimizerStrategy','2','--cminDefaultMinimizerTolerance','0.00001','--X-rtd','MINIMIZER_freezeDisassociatedParams','--setParameters',','.join(f'{k}={v}' for k,v in params.items()),'-n','_check'];(dd/'command.json').write_text(json.dumps(cmd))
 with (dd/'combine.log').open('w') as out:subprocess.run(cmd,cwd=dd,stdout=out,stderr=subprocess.STDOUT)
 if not (dd/'multidimfit_check.root').exists():rows.append(dict(start=r,failed=True));continue
 ff=ROOT.TFile.Open(str(dd/'multidimfit_check.root'));res=ff.Get('fit_mdf');values={p.GetName():p.getVal() for p in res.floatParsFinal()}
 for k,v in values.items():w.var(k).setVal(v)
 rows.append(dict(start=r,status=res.status(),covQual=res.covQual(),nll=nll.getVal(),params=values));ff.Close()
(d/'multistart.json').write_text(json.dumps(rows,indent=2));print([(r['start'],r.get('params',{}).get('r'),r.get('nll'),r.get('status')) for r in rows])
