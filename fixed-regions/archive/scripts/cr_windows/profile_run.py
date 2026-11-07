from pathlib import Path
import json,subprocess
P=Path.cwd();R=P/'runs/20261029_cr_windows';eras=['2022','2022EE','2023','2023BPix']
for tag in json.loads((R/'control_comparison.json').read_text()):
 d=R/tag/'degree3'
 for role,ws,poi,snapshot in [('Y2_profile',d/'gof/workspace.root',','.join('CMS_EXO25020_bkgScale_upsilon2s_'+e for e in eras),'control_fit'),('tail_profile',d/'tail_profile.root','SR_upsilon_events','profile_start')]:
  dest=d/role;dest.mkdir(exist_ok=True);c=[str(P/'bin/combine'),'-M','MultiDimFit','-d',str(ws),'-m','7.5','--algo','singles','--redefineSignalPOIs',poi,'--floatOtherPOIs','1','--snapshotName',snapshot,'--cl','0.68268949','--cminDefaultMinimizerStrategy','1','--cminDefaultMinimizerTolerance','0.00001','--X-rtd','MINIMIZER_freezeDisassociatedParams','--saveFitResult','-n','_profile','-v','1']
  if role=='Y2_profile':c+=['--setParameters','r=0','--freezeParameters','r']
  (dest/'command.json').write_text(json.dumps(c,indent=2))
  with (dest/'combine.log').open('w') as f:res=subprocess.run(c,cwd=dest,stdout=f,stderr=subprocess.STDOUT)
  print(tag,role,res.returncode,flush=True);assert res.returncode==0
