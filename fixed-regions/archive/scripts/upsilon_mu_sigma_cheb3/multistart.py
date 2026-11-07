from pathlib import Path
import json,subprocess,concurrent.futures
P=Path.cwd();R=P/'runs/20261030_upsilon_mu_sigma_cheb3';A=json.loads((R/'configs.json').read_text())
def job(item):
 a,i=item;d=R/a['tag']/f'start{i}';d.mkdir(exist_ok=True);N=a['events'];seeds=[(a['nominal_mu'],a['nominal_sigma'],.7,.01,-.15,.0,.01),(9.43,.20,.5,1,-.3,.1,-.05),(9.44,.12,.3,5,-.5,.05,.05)]
 mu,si,frac,y2,c1,c2,c3=seeds[i];params=f'upsilon_mu={mu},upsilon_sigma={si},N_upsilon1s={N*frac},N_continuum={N*(1-frac)},scale_upsilon2s={y2},cheb_c1={c1},cheb_c2={c2},cheb_c3={c3}'
 cmd=[str(P/'bin/combine'),'-M','MultiDimFit','-d',str(d.parent/'workspace.root'),'--algo','none','--redefineSignalPOIs','upsilon_mu,upsilon_sigma','--noMCbonly','1','--saveFitResult','--saveWorkspace','--cminDefaultMinimizerStrategy','1','--cminDefaultMinimizerTolerance','0.001','--setParameters',params,'-n','_seed','-v','1'];(d/'command.json').write_text(json.dumps(cmd))
 with (d/'combine.log').open('w') as f:subprocess.run(cmd,cwd=d,stdout=f,stderr=subprocess.STDOUT)
with concurrent.futures.ThreadPoolExecutor(max_workers=2) as ex:list(ex.map(job,[(a,i) for a in A for i in range(3)]))
