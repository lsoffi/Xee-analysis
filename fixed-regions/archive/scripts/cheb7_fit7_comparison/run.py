from pathlib import Path
import subprocess,json,time,concurrent.futures,os,hashlib
P=Path.cwd();OLD=P/'runs/20261030_cheb7_bias_scan';R=P/'runs/20261030_cheb7_fit7_comparison';env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1')
def execute(cmd,d):
 d.mkdir(parents=True,exist_ok=True);(d/'command.json').write_text(json.dumps(cmd,indent=2));t=time.time()
 with (d/'combine.log').open('w') as out:res=subprocess.run(cmd,cwd=d,stdout=out,stderr=subprocess.STDOUT,timeout=1800,env=env)
 (d/'receipt.json').write_text(json.dumps(dict(exit_code=res.returncode,seconds=time.time()-t)));assert res.returncode==0,(d,res.returncode)
def job(m):
 src=OLD/('m'+m);d=R/('m'+m);d.mkdir(exist_ok=True);a=json.loads((src/'study.json').read_text());a['workspace']=str(src/'generator.root');a['workspace_sha256']=hashlib.sha256((src/'generator.root').read_bytes()).hexdigest();a['toys_sources']=[]
 for label in ['zero','injected']:
  for kind in ['toys','asimov']:
   path=src/label/('fixed' if kind=='toys' else 'asimov/fixed');cmd=json.loads((path/'combine.log.command.json').read_text());cmd[cmd.index('-d')+1]=a['workspace'];cmd+=['--trackParameters','rgx{local_c[1-7]_.*}'];tf=Path(cmd[cmd.index('--toysFile')+1]);a['toys_sources'].append(dict(scenario=label,kind=kind,path=str(tf),sha256=hashlib.sha256(tf.read_bytes()).hexdigest()));execute(cmd,d/label/kind)
  print(m,label,'complete',flush=True)
 (d/'study.json').write_text(json.dumps(a,indent=2));print(m,'DONE',flush=True)
with concurrent.futures.ThreadPoolExecutor(max_workers=2) as ex:list(ex.map(job,[f'{i/10:.1f}' for i in range(65,86)]))
