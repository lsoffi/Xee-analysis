from pathlib import Path
import sys,subprocess,json,time,concurrent.futures,os
P=Path.cwd();S=P/'scripts/upper_envelope_scan';R=P/'runs/20261030_upper_envelope_scan';PY=str(P/'bin/python');C=str(P/'bin/combine')
def execute(cmd,d,name,project=False):
    d.mkdir(parents=True,exist_ok=True);receipt=d/(name+'.receipt.json')
    if receipt.exists() and json.loads(receipt.read_text())['exit_code']==0:return
    (d/(name+'.command.json')).write_text(json.dumps(cmd,indent=2));t=time.time()
    with (d/(name+'.log')).open('w') as h:rc=subprocess.run(cmd,cwd=P if project else d,stdout=h,stderr=subprocess.STDOUT,timeout=3600,env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1')).returncode
    receipt.write_text(json.dumps(dict(exit_code=rc,seconds=time.time()-t)));print(d.name,name,rc,round(time.time()-t,1),flush=True)
    if rc:raise RuntimeError(str(d)+' '+name)
def job(arg):
    deg,m=arg;D=R/f'degree{deg}/m{m}';J=D/'joint'
    for region in ['signal','control']:execute([PY,str(S/'prepare.py'),m,str(deg),region],D,'prepare_'+region,True)
    execute([PY,str(S/'build.py'),m,str(deg)],D,'build',True);execute([PY,str(S/'seed.py'),m,str(deg)],D,'seed',True)
    common=['-d',str(J/'workspace_seeded.root'),'-m',m,'-n','_upper','--cminDefaultMinimizerStrategy','2','--cminDefaultMinimizerTolerance','0.001','--cminDefaultMinimizerPrecision','1e-12','--X-rtd','MINIMIZER_freezeDisassociatedParams','--cminRunAllDiscreteCombinations','-v','1']
    execute([C,'-M','FitDiagnostics']+common+['--saveShapes','--saveNormalizations','--keepFailures','--robustFit','0','--minos','none','--rMin','-100','--rMax','100'],J/'fitdiagnostics','combine')
    execute([C,'-M','Significance']+common+['--redefineSignalPOIs','r','--rMin','-100','--rMax','100'],J/'significance','combine')
    lim=[C,'-M','AsymptoticLimits']+common+['--rMax','100','--rRelAcc','0.001','--rAbsAcc','0.00001'];lim[lim.index('--cminDefaultMinimizerTolerance')+1]='0.00001'
    execute(lim,J/'limits','combine')
    gof=common.copy();gof[gof.index('--cminDefaultMinimizerStrategy')+1]='0'
    execute([C,'-M','GoodnessOfFit']+gof+['--algo','saturated','--rMin','-100','--rMax','100','--trackParameters','r'],J/'gof_observed','combine')
    execute([PY,str(S/'extract.py'),m,str(deg)],D,'extract',True)
if __name__=='__main__':
    jobs=[(d,f'{m/10:.1f}') for d in [1,2] for m in range(105,111)]
    if len(sys.argv)>1:
        jobs=[(d,f'{m/10:.1f}') for d in [1,2] for m in range(111,121)] if sys.argv[1]=='--extended' else [(int(sys.argv[1]),sys.argv[2])]
    with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:list(pool.map(job,jobs))
