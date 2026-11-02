from pathlib import Path
import subprocess,json,time
P=Path.cwd();R=P/'runs/20261030_sr6_cr4_scan';C=str(P/'bin/combine')
for it in range(180):
 for rr in sorted(R.glob('m*/results.json')):
  d=rr.parent/'joint/gof_sb'
  if (d/'receipt.json').exists():continue
  a=json.loads(rr.read_text());m=f"{a['mass']:.1f}";d.mkdir(exist_ok=True);cmd=[C,'-M','GoodnessOfFit','--algo','saturated','-d',str(rr.parent/'joint/workspace.root'),'-m',m,'-n','_sb','--rMin','-30','--rMax','30','--cminDefaultMinimizerStrategy','1','--cminDefaultMinimizerTolerance','0.001','--X-rtd','MINIMIZER_freezeDisassociatedParams','--cminRunAllDiscreteCombinations'];(d/'command.json').write_text(json.dumps(cmd));t=time.time()
  with (d/'combine.log').open('w') as out:res=subprocess.run(cmd,cwd=d,stdout=out,stderr=subprocess.STDOUT,timeout=600)
  (d/'receipt.json').write_text(json.dumps(dict(exit_code=res.returncode,seconds=time.time()-t)));print(m,res.returncode,flush=True)
 if len(list(R.glob('m*/joint/gof_sb/receipt.json')))==31:break
 time.sleep(20)
