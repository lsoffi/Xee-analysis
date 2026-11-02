from pathlib import Path
import json,subprocess,concurrent.futures,time,os
P=Path.cwd();R=P/'runs/20261030_cheb7_bias_scan';S=P/'scripts/cheb7_bias_scan';PY=str(P/'bin/python');C=str(P/'bin/combine')
def run(cmd,d,log='combine.log',cwd=None):
 d.mkdir(parents=True,exist_ok=True);(d/(log+'.command.json')).write_text(json.dumps(cmd))
 with (d/log).open('w') as out:r=subprocess.run(cmd,cwd=cwd or d,stdout=out,stderr=subprocess.STDOUT,timeout=300,env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1'))
 assert r.returncode==0,(cmd,r.returncode)
def job(m):
 d=R/('m'+m)
 while not (d/'study.json').exists():time.sleep(2)
 a=json.loads((d/'study.json').read_text())
 for label,inj in [('zero',0),('injected',a['injection'])]:
  dd=d/label/'asimov';opts=['-m',m,'--X-rtd','MINIMIZER_freezeDisassociatedParams'];run([C,'-M','GenerateOnly','-d',str(d/'generator.root'),'--snapshotName','generator_truth','--saveToys','--toysFrequentist','--bypassFrequentistFit','--expectSignal',str(inj),'--genBinnedChannels',','.join(f'{reg}_ch{i}' for reg in ['SR','CR'] for i in range(1,5)),'-t','-1','-s','307701','-n','_asimov']+opts,dd/'generate');tf=next((dd/'generate').glob('higgsCombine*.root'));crop=dd/'sliding.root';run([PY,str(S/'crop_toys.py'),str(tf),str(crop),str(d/'fit_sliding.root')],dd,'crop.log',P)
  for model,toys in [('fixed',tf),('sliding',crop)]:run([C,'-M','FitDiagnostics','-d',str(d/f'fit_{model}.root'),'--toysFile',str(toys),'--toysFrequentist','-t','-1','--rMin','-100','--rMax','100','--keepFailures','--saveShapes','--saveNormalizations','--minos','none','--cminDefaultMinimizerStrategy','1','--cminDefaultMinimizerTolerance','0.001','-n','_asimov','-v','1']+opts,dd/model)
 print(m,'Asimov done',flush=True)
with concurrent.futures.ThreadPoolExecutor(max_workers=2) as ex:list(ex.map(job,[f'{i/10:.1f}' for i in range(65,86)]))
