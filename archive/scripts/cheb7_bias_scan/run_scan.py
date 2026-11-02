from pathlib import Path
import subprocess,json,time,os,concurrent.futures
P=Path.cwd();R=P/'runs/20261030_cheb7_bias_scan';S=P/'scripts/cheb7_bias_scan';C=str(P/'bin/combine');PY=str(P/'bin/python');env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1')
def execute(cmd,d,log='combine.log',cwd=None):
 d.mkdir(parents=True,exist_ok=True);t=time.time();(d/(log+'.command.json')).write_text(json.dumps(cmd,indent=2))
 with (d/log).open('w') as out:res=subprocess.run(cmd,cwd=cwd or d,stdout=out,stderr=subprocess.STDOUT,timeout=900,env=env)
 receipt=dict(exit_code=res.returncode,seconds=time.time()-t);(d/(log+'.receipt.json')).write_text(json.dumps(receipt));assert res.returncode==0,(cmd,res.returncode);return receipt

def job(m):
 d=R/('m'+m);d.mkdir(exist_ok=True)
 for region in ['signal','control']:execute([PY,str(S/'prepare.py'),m,region],d,region+'_prepare.log',P)
 execute([PY,str(S/'build_joint.py'),m],d,'build.log',P);execute([PY,str(S/'make_workspaces.py'),m],d,'workspaces.log',P);a=json.loads((d/'study.json').read_text())
 for label,inj in [('zero',0),('injected',a['injection'])]:
  dd=d/label;seed=307000+int(float(m)*100)+(label=='injected');opts=['-m',m,'--X-rtd','MINIMIZER_freezeDisassociatedParams'];g=dd/'generate';execute([C,'-M','GenerateOnly','-d',str(d/'generator.root'),'--snapshotName','generator_truth','--saveToys','--toysFrequentist','--bypassFrequentistFit','--expectSignal',str(inj),'--genBinnedChannels',','.join(f'{reg}_ch{i}' for reg in ['SR','CR'] for i in range(1,5)),'-t','30','-s',str(seed),'-n','_toys']+opts,g);tf=next(g.glob('higgsCombine*.root'));cropped=dd/'sliding_toys.root';execute([PY,str(S/'crop_toys.py'),str(tf),str(cropped),str(d/'fit_sliding.root')],dd,'crop.log',P)
  for method,toysfile in [('fixed',tf),('sliding',cropped)]:
   cmd=[C,'-M','FitDiagnostics','-d',str(d/f'fit_{method}.root'),'--toysFile',str(toysfile),'--toysFrequentist','-t','30','-s',str(seed),'--rMin','-100','--rMax','100','--keepFailures','--minos','none','--cminDefaultMinimizerStrategy','1','--cminDefaultMinimizerTolerance','0.001','-n','_toys','-v','1']+opts;execute(cmd,dd/method)
  print(m,label,'complete',flush=True)
 print(m,'DONE',flush=True)
with concurrent.futures.ThreadPoolExecutor(max_workers=2) as ex:list(ex.map(job,[f'{i/10:.1f}' for i in range(65,86)]))
