import ROOT,pathlib,json,subprocess,time,re,concurrent.futures
ROOT.gSystem.Load('libHiggsAnalysisCombinedLimit');ROOT.gROOT.SetBatch(True);P=pathlib.Path('/Users/imb/Documents/Codex/Analisi-statistica-Xee');R=P/'runs/20261029_gaussian_signal_scan'
for mass in ['10.9','11.0']:
 f=ROOT.TFile.Open(str(R/f'workspaces/dcb/m{mass}.root'));w=f.Get('w');reps=[];keep=[];discretes=w.genobj('discreteParams');discretes.removeAll()
 for i,era in enumerate(['2022','2022EE','2023','2023BPix'],1):
  pdf=w.pdf('shapeBkg_dy_ch'+str(i)).getPdf(1);cat=ROOT.RooCategory('fixed_cheb7_index_'+era,'');one=ROOT.RooMultiPdf('fixed_cheb7_'+era,'',cat,ROOT.RooArgList(pdf));getattr(w,'import')(one);keep.extend([cat,one]);discretes.add(w.cat(cat.GetName()));reps.append(pdf.GetName()+'='+one.GetName())
  x=w.var('mass');obs=ROOT.RooArgSet(x)
  for value in [6.1,7.5,9.5,11.9]:
   x.setVal(value);assert abs(pdf.getVal(obs)-w.pdf(one.GetName()).getVal(obs))<1e-12
 for name in ['ModelConfig','ModelConfig_bonly']:
  mc=w.obj(name);old=mc.GetPdf();new=w.factory('EDIT::'+old.GetName()+'_fixedwrapper('+old.GetName()+','+','.join(reps)+')');assert new;mc.SetPdf(new)
 w.writeToFile(str(R/f'workspaces/dcb/m{mass}_fixedwrapper.root'))
def run(mass):
 base=json.loads((R/f'limit_repairs/dcb/m{mass}/command.json').read_text());cmd=base['command'];cmd[cmd.index('-d')+1]=str(R/f'workspaces/dcb/m{mass}_fixedwrapper.root');cmd[cmd.index('-n')+1]='_fixedwrapper';wd=R/f'limit_wrapper_checks/dcb/m{mass}';wd.mkdir(parents=True,exist_ok=True);base.update(command=cmd,cwd=str(wd));(wd/'command.json').write_text(json.dumps(base,indent=2));start=time.time()
 with (wd/'combine.log').open('w') as f:
  try:rc=subprocess.run(cmd,cwd=wd,stdout=f,stderr=subprocess.STDOUT,timeout=180).returncode
  except subprocess.TimeoutExpired:rc=124
 matches=re.findall(r'At r =\s*([\d.eE+-]+):[^\n]*CLs\s*=\s*([\d.eE+-]+)',(wd/'combine.log').read_text());base.update(exit_code=rc,elapsed_seconds=time.time()-start,last_cls=matches[-1] if matches else None);(wd/'receipt.json').write_text(json.dumps(base,indent=2));print(mass,base['last_cls'],flush=True)
with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:list(pool.map(run,['10.9','11.0']))
