from pathlib import Path
import subprocess,json,concurrent.futures,time
P=Path.cwd();R=P/'runs/20261030_quick_no_upsilon';C=str(P/'bin/combine')
def run(cmd,d):
 d.mkdir(exist_ok=True,parents=True);(d/'command.json').write_text(json.dumps(cmd,indent=2));t=time.time()
 with (d/'combine.log').open('w') as out:r=subprocess.run(cmd,cwd=d,stdout=out,stderr=subprocess.STDOUT,timeout=300)
 (d/'receipt.json').write_text(json.dumps(dict(exit_code=r.returncode,seconds=time.time()-t)))
def job(pair):
 a,label=pair;m=a['mass'];mu=0 if label=='zero' else a['injection'];d=R/('m'+m)/label;seed=30000+int(float(m)*100)+(label=='injected');opts=['-m',m,'--X-rtd','MINIMIZER_freezeDisassociatedParams'];N=a['toys_per_injection']
 run([C,'-M','GenerateOnly','-d',str(d.parent/'generator.root'),'--snapshotName','truth','--saveToys','--toysFrequentist','--bypassFrequentistFit','--expectSignal',str(mu),'--genBinnedChannels','ch1,ch2,ch3,ch4','-t',str(N),'-s',str(seed),'-n','_toys']+opts,d/'generate')
 tf=d/'generate'/f'higgsCombine_toys.GenerateOnly.mH{m}.{seed}.root'
 if not tf.exists():tf=next((d/'generate').glob('higgsCombine*.root'))
 run([C,'-M','FitDiagnostics','-d',str(d.parent/'fit_noY.root'),'--toysFile',str(tf),'--toysFrequentist','-t',str(N),'-s',str(seed),'--rMin','-5','--rMax','5','--keepFailures','--minos','none','--skipBOnlyFit','--cminDefaultMinimizerStrategy','1','--cminDefaultMinimizerTolerance','0.001','-n','_toys','-v','1']+opts,d/'fit');print(m,label,'done',flush=True)
with concurrent.futures.ThreadPoolExecutor(max_workers=2) as ex:list(ex.map(job,[(a,l) for a in json.loads((R/'cases.json').read_text()) for l in ['zero','injected']]))
