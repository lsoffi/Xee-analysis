from pathlib import Path
import subprocess, json, time, os, sys
P=Path.cwd();S=P/'scripts/upper_sr_cr';R=P/'runs/20261030_upper_sr_cr'
PY=str(P/'bin/python');C=str(P/'bin/combine')
def run(cmd,wd,name):
    wd.mkdir(parents=True,exist_ok=True)
    (wd/(name+'.command.json')).write_text(json.dumps(cmd,indent=2));t=time.time()
    with (wd/(name+'.log')).open('w') as h:
        res=subprocess.run(cmd,cwd=P if name.startswith('prepare') or name in ['build','seed'] else wd,stdout=h,stderr=subprocess.STDOUT,timeout=600,env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1'))
    receipt=dict(exit_code=res.returncode,seconds=time.time()-t)
    (wd/(name+'.receipt.json')).write_text(json.dumps(receipt))
    print(wd.name,name,receipt,flush=True)
    if res.returncode:raise RuntimeError(str(wd)+' '+name)
for degree in map(int,sys.argv[1:] or ['1','2','3']):
    D=R/f'degree{degree}'
    for region in ['signal','control']:run([PY,str(S/'prepare.py'),str(degree),region],D,'prepare_'+region)
    run([PY,str(S/'build.py'),str(degree)],D,'build')
    run([PY,str(S/'seed.py'),str(degree)],D,'seed')
    common=['-d',str(D/'joint/workspace_seeded.root'),'-m','10.8','-n','_upper','--cminDefaultMinimizerStrategy','2','--cminDefaultMinimizerTolerance','0.001','--cminDefaultMinimizerPrecision','1e-12','--X-rtd','MINIMIZER_freezeDisassociatedParams','-v','1']
    run([C,'-M','FitDiagnostics']+common+['--skipSBFit','--saveShapes','--saveNormalizations','--keepFailures','--robustFit','0','--minos','none'],D/'joint/fitdiagnostics','combine')
    run([C,'-M','GoodnessOfFit']+common+['--algo','saturated','--fixedSignalStrength','0'],D/'joint/gof_observed','combine')
