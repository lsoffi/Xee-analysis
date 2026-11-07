from pathlib import Path
import subprocess,json,time,concurrent.futures,os,sys
P=Path.cwd();S=P/'scripts/sr6_cr4_scan';R=P/'runs/20261030_sr6_cr4_scan';C=str(P/'bin/combine');PY=str(P/'bin/python')
def run(cmd,d,log,project=False):
 d.mkdir(parents=True,exist_ok=True);(d/(log+'.command.json')).write_text(json.dumps(cmd));t=time.time()
 with (d/log).open('w') as out:r=subprocess.run(cmd,cwd=P if project else d,stdout=out,stderr=subprocess.STDOUT,timeout=1800)
 receipt=dict(exit_code=r.returncode,seconds=time.time()-t);(d/(log+'.receipt.json')).write_text(json.dumps(receipt));return r.returncode

def job(arg):
 mass,family=arg;d=R/'bias'/family/('m'+mass)
 if all((d/f'r{r}'/'fitDiagnostics_bias.root').exists() for r in [0,1]):return
 # Wait for the mass workspace to be constructed by the data scan.
 while not (R/('m'+mass)/'joint/workspace_unseeded.root').exists():time.sleep(5)
 if run([PY,str(S/'bias_prepare.py'),mass,family],d,'prepare.log',True):return print('PREP FAIL',mass,family,flush=True)
 for rtrue in [0,1]:
  wd=d/f'r{rtrue}';opts=['-d',str(d/'generator.root'),'-m',mass,'--snapshotName','generation'];cmd=[C,'-M','GenerateOnly','--saveToys','--toysFrequentist','--bypassFrequentistFit','--expectSignal',str(rtrue),'--genBinnedChannels',','.join(f'{region}_ch{i}' for region in ['SR','CR'] for i in range(1,5)),'-t','-1','-n','_asimov']+opts
  if run(cmd,wd,'generate.log'):continue
  tf=next(wd.glob('higgsCombine_asimov*.root'));cmd=[C,'-M','FitDiagnostics','-d',str(R/('m'+mass)/'joint/workspace.root'),'-m',mass,'-n','_bias','-t','-1','--toysFile',str(tf),'--expectSignal',str(rtrue),'--skipBOnlyFit','--minos','none','--robustFit','0','--keepFailures','--preFitValue',str(rtrue),'--rMin','-30','--rMax','30','--cminDefaultMinimizerStrategy','1','--cminDefaultMinimizerTolerance','0.001','--X-rtd','MINIMIZER_freezeDisassociatedParams','--cminRunAllDiscreteCombinations']
  run(cmd,wd,'fit.log');print('BIAS',mass,family,rtrue,flush=True)
if __name__=='__main__':
 jobs=[(f'{k/10:.1f}',f) for k in range(60,91) for f in ['bern7','cheb7']]
 if len(sys.argv)>1:jobs=[(sys.argv[1],f) for f in ['bern7','cheb7']]
 with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:list(pool.map(job,jobs))
 print('BIAS END',flush=True)
