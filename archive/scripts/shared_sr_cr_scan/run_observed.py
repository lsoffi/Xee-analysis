"""Resumable native Combine scan; preserve all commands, trials and failures."""
from pathlib import Path
import argparse,concurrent.futures,datetime,json,math,os,signal,subprocess,time,traceback

P=Path(__file__).resolve().parents[2];S=P/'scripts/shared_sr_cr_scan';R=P/'runs/20261031_shared_sr_cr_scan'
PY=str(P/'bin/python');C=str(P/'bin/combine')

def atomic(path,value):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    q=path.with_suffix(path.suffix+'.tmp');q.write_text(json.dumps(value,indent=2));q.replace(path)

def execute(cmd,directory,name='combine',timeout=3600,project=False):
    directory=Path(directory);directory.mkdir(parents=True,exist_ok=True)
    receipt=directory/(name+'.receipt.json')
    if receipt.exists():return json.loads(receipt.read_text())
    atomic(directory/(name+'.command.json'),cmd)
    start=time.time()
    with (directory/(name+'.log')).open('w') as log:
        proc=subprocess.Popen(cmd,cwd=P if project else directory,stdout=log,stderr=subprocess.STDOUT,
                              env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1'),start_new_session=True)
        try:rc=proc.wait(timeout=timeout)
        except subprocess.TimeoutExpired:
            os.killpg(proc.pid,signal.SIGTERM)
            try:proc.wait(timeout=5)
            except subprocess.TimeoutExpired:os.killpg(proc.pid,signal.SIGKILL);proc.wait()
            rc=124
    out=dict(exit_code=rc,seconds=time.time()-start,timeout_seconds=timeout,finished_at=datetime.datetime.now().astimezone().isoformat())
    atomic(receipt,out);print(directory.relative_to(R),name,rc,round(out['seconds'],1),flush=True)
    return out

def historical_expected(mass):
    candidates=[P/f'runs/20261030_upper_envelope_scan/degree3/m{mass}/results.json',
                P/f'runs/20261030_sr6_cr4_scan/m{mass}/results.json']
    for path in candidates:
        if path.exists():
            a=json.loads(path.read_text());value=a.get('limits',{}).get('expected_50')
            if value and math.isfinite(value) and value>0:return float(value)
    return 1.

def extract(mass):
    # Extraction is intentionally refreshed after every method.
    result=subprocess.run([PY,str(S/'extract_observed.py'),mass],cwd=P,capture_output=True,text=True,
                          env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1'),timeout=180)
    D=R/'observed'/('m'+mass)
    with (D/'extraction.log').open('a') as h:h.write(result.stdout+result.stderr)
    if result.returncode:raise RuntimeError('Extraction failed at '+mass+'; see '+str(D/'extraction.log'))
    return json.loads((D/'results.json').read_text())

def common(workspace,mass,strategy=2,tolerance=.001):
    return ['-d',str(workspace),'-m',mass,'-n','_shared','--cminDefaultMinimizerStrategy',str(strategy),
            '--cminDefaultMinimizerTolerance',str(tolerance),'--cminDefaultMinimizerPrecision','1e-12',
            '--X-rtd','MINIMIZER_freezeDisassociatedParams','--cminRunAllDiscreteCombinations','-v','1']

def job(mass,prepare=False,fits_only=False):
    mass=f'{float(mass):.1f}';D=R/'observed'/('m'+mass);J=D/'joint'
    try:
        if prepare:
            prep=execute([PY,str(S/'model.py'),mass,'--degree','6','--family','envelope'],D,'prepare',project=True)
            if prep['exit_code']:raise RuntimeError('Builder failed')
        workspace=J/'workspace_seeded.root'
        if not workspace.exists():raise RuntimeError('Missing builder workspace '+str(workspace))
        E=historical_expected(mass);rbound=max(100.,5*E)
        atomic(D/'observed_protocol.json',dict(mass=float(mass),historical_expected_for_initialization=E,
          initial_r_range=[-rbound,rbound],anchors=[0,E],native_methods=True,expected_bands='native provisional',
          minimizer=dict(strategy=2,tolerance=.001,precision=1e-12),
          selection='smallest absolute extended deviance plus auxiliary penalties among numerically valid point fits'))
        for label,start in [('zero',0.),('expected',E)]:
            cmd=[C,'-M','FitDiagnostics']+common(workspace,mass)+['--saveShapes','--saveNormalizations',
                 '--keepFailures','--robustFit','0','--minos','none','--rMin',str(-rbound),'--rMax',str(rbound),
                 '--preFitValue',str(start)]
            if label=='expected':cmd+=['--skipBOnlyFit']
            execute(cmd,J/('fitdiagnostics_anchor_'+label),timeout=1200)
            a=extract(mass)
        # Use independently selected best-fit snapshots as seeds only.
        sb=J/'workspace_postfit_s.root';b=J/'workspace_postfit_b.root'
        sig=execute([C,'-M','Significance']+common(sb,mass)+['--redefineSignalPOIs','r','--rMin',str(-rbound),'--rMax',str(rbound)],J/'significance',timeout=1200)
        a=extract(mass)
        if fits_only:
            atomic(D/'fits_completion.json',dict(status='fits_and_significance_complete',flags=a['flags']))
            return mass,'fits_and_significance_complete'
        sigma=a.get('r_error') or 0.
        if not math.isfinite(sigma) or sigma<0:sigma=0.
        # The native algorithm floors its initial r step at 2% of the range.
        # Reusing the broad signed diagnostic range is unnecessarily expensive
        # and gives a poor initial bracket for small signals.
        limitbound=max(1.,8*sigma,5*abs(a.get('r',0.)),5*E)
        limitbound=min(1e6,limitbound)
        atomic(D/'limits_protocol.json',dict(strategy=0,tolerance=.001,rMax=limitbound,
          reason='Range-scaled initial bracket; observed CLs crossing independently checked; native expected bands remain provisional'))
        limits=execute([C,'-M','AsymptoticLimits']+common(b,mass,strategy=0,tolerance=.001)+
            ['--rMax',str(limitbound),'--rRelAcc','.001','--rAbsAcc','.00001','--setParameters','r=0'],J/'limits',timeout=7200)
        a=extract(mass)
        expected=a.get('limits',{}).get('expected_50')
        gb=max(100.,5*expected) if expected and math.isfinite(expected) else rbound
        atomic(D/'gof_protocol.json',dict(r_range=[-gb,gb],strategy=0,tolerance=.001,precision=1e-12,
                                        r_floating=True,algorithm='saturated',seed='observed postfit B'))
        execute([C,'-M','GoodnessOfFit']+common(b,mass,strategy=0)+['--algo','saturated','--rMin',str(-gb),
                 '--rMax',str(gb),'--setParameters','r=0','--trackParameters','r'],J/'gof_observed',timeout=1200)
        a=extract(mass)
        atomic(D/'completion.json',dict(status=a['status'],numerical_ok=a['numerical_ok'],flags=a['flags'],
              completed_at=datetime.datetime.now().astimezone().isoformat()))
        return mass,a['status']
    except Exception as e:
        atomic(D/'worker_failure.json',dict(error=str(e),traceback=traceback.format_exc(),mass=float(mass)))
        print('FAILED',mass,str(e),flush=True);return mass,'failed'

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--masses',nargs='+');parser.add_argument('--workers',type=int,default=2)
    parser.add_argument('--prepare',action='store_true');parser.add_argument('--fits-only',action='store_true');args=parser.parse_args()
    preferred=['7.5','11.5','9.5','6.0','12.0']
    masses=args.masses or preferred+[f'{n/10:.1f}' for n in range(60,121) if f'{n/10:.1f}' not in preferred]
    assert 1<=args.workers<=4
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as pool:
        finished=list(pool.map(lambda m:job(m,args.prepare,args.fits_only),masses))
    atomic(R/('observed_fit_runner.json' if args.fits_only else 'observed_runner.json'),dict(status='completed',results=finished,workers=args.workers))
