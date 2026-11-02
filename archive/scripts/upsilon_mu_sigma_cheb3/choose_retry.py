from pathlib import Path
import ROOT,json,subprocess,shutil
ROOT.gSystem.Load('libHiggsAnalysisCombinedLimit');P=Path.cwd();R=P/'runs/20261030_upsilon_mu_sigma_cheb3'
for a in [a for a in json.loads((R/'configs.json').read_text()) if a['tag'] in ['k2_2022EE_b3','k3_2022_b3']]:
 d=R/a['tag'];wf=ROOT.TFile.Open(str(d/'workspace.root'));w=wf.Get('w');nll=w.pdf('model').createNLL(w.data('data_obs'),ROOT.RooFit.Extended(True),ROOT.RooFit.Offset(False));rows=[]
 for path in [d/'multidimfit_shape.root']+list(d.glob('start*/multidimfit_seed.root')):
  if not path.exists():continue
  f=ROOT.TFile.Open(str(path));r=f.Get('fit_mdf')
  if not r:continue
  vals={p.GetName():p.getVal() for p in r.floatParsFinal()}
  for k,v in vals.items():w.var(k).setVal(v)
  rows.append(dict(path=str(path),status=r.status(),covQual=r.covQual(),nll=nll.getVal(),params=vals));f.Close()
 valid=[r for r in rows if r['status'] in [0,1] and r['nll']==r['nll']]
 if not valid:print('NO FIT',a['tag'],flush=True);(d/'selection.json').write_text(json.dumps(rows,indent=2));continue
 best=min(valid,key=lambda r:r['nll']);(d/'selection.json').write_text(json.dumps(dict(attempts=rows,selected=best),indent=2));archive=d/'initial';archive.mkdir(exist_ok=True)
 for name in ['multidimfit_shape.root','higgsCombine_shape.MultiDimFit.mH120.root','combine.log','command.json','combine_logger.out']:
  if (d/name).exists():shutil.move(str(d/name),str(archive/name))
 cmd=[str(P/'bin/combine'),'-M','MultiDimFit','-d',str(d/'workspace.root'),'--algo','singles','--redefineSignalPOIs','upsilon_mu,upsilon_sigma','--floatOtherPOIs','1','--noMCbonly','1','--saveFitResult','--saveWorkspace','--cl','0.68268949','--cminDefaultMinimizerStrategy','0','--cminDefaultMinimizerTolerance','0.001','--setParameters',','.join(f'{k}={v}' for k,v in best['params'].items()),'-n','_shape','-v','1'];(d/'command.json').write_text(json.dumps(cmd))
 with (d/'combine.log').open('w') as out:subprocess.run(cmd,cwd=d,stdout=out,stderr=subprocess.STDOUT)
 if not (d/'multidimfit_shape.root').exists():shutil.copy(best['path'] if Path(best['path']).exists() else str(archive/'multidimfit_shape.root'),d/'multidimfit_shape.root');(d/'profile_missing').touch()
 print(a['tag'],best['nll'],flush=True);wf.Close()
