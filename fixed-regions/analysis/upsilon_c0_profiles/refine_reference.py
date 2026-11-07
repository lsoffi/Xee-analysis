"""One unconditional refinement seeded at a deeper valid conditional minimum.

Default is read-only planning. --run requires a closed main profile runner.
Old workspaces/fits/results are never modified. No automatic numerical retries.
"""
from pathlib import Path
import argparse, concurrent.futures, fcntl, json, math, os, signal, subprocess, sys, time
import ROOT

P=Path(__file__).resolve().parents[2]
S=P/'analysis/upsilon_c0_profiles'
R=P/'runs/20261102_upsilon_c0_profiles'
sys.path.insert(0,str(S));sys.path.insert(0,str(P/'analysis/separate_sr_c0_limits'))
from fit_helpers import read_stage, dump
from run import read, sha, stamp, MASSES
ROOT.gROOT.SetBatch(True);ROOT.gSystem.Load('libHiggsAnalysisCombinedLimit')
ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.WARNING)


def plan(side,mass):
    D=R/side/f'm{mass:.1f}';config=read(D/'configuration.json');profile=read(D/'profile.json')
    if profile['status']!='complete':raise RuntimeError('Conditional profile not complete')
    old=[q for q in profile['references'] if q['kind'].startswith('previous_')]
    points=[q for q in profile['points'] if q.get('point_valid')]
    if not old or not points:raise RuntimeError('Missing valid reference or conditional points')
    bestold=min(old,key=lambda q:q['exact_deviance']);best=min(points,key=lambda q:q['exact_deviance'])
    improvement=bestold['exact_deviance']-best['exact_deviance']
    return dict(side=side,mass=mass,needed=improvement>.001,trigger_threshold=.001,
        previous_physical_reference=bestold,conditional_start={k:best[k] for k in ['r','exact_deviance','source']},
        improvement=improvement,workspace=config['workspace'],workspace_sha256=config['workspace_sha256'],
        r_range=[0.,max(100.,1.1*config['rmax'])],profile=str(D/'profile.json'))


def run_one(side,mass):
    p=plan(side,mass);D=R/side/f'm{mass:.1f}';J=D/'refined_reference';J.mkdir(exist_ok=True)
    with (J/'job.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX)
        if (J/'completion.json').exists():return read(J/'completion.json')
        dump(J/'plan.json',p)
        if not p['needed']:
            out=dict(side=side,mass=mass,status='not_needed',improvement=p['improvement'],finished_at=stamp())
            dump(J/'completion.json',out);return out
        assert sha(p['workspace'])==p['workspace_sha256']
        source=ROOT.TFile.Open(p['conditional_start']['source']);fit=source.Get('fit_b')
        assert fit,'Conditional source has no fit_b'
        seed={v.GetName():float(v.getVal()) for v in fit.floatParsFinal()};seed['r']=p['conditional_start']['r'];source.Close()
        parameters=','.join(f'{name}={value:.17g}' for name,value in sorted(seed.items()))
        cmd=[str(P/'bin/combine'),'-M','FitDiagnostics','-d',p['workspace'],'-m',f'{mass:.1f}',
            '-n','_profile_reference','--skipBOnlyFit','--keepFailures','--robustFit','0','--minos','none',
            '--rMin','0','--rMax',str(p['r_range'][1]),'--setParameters',parameters,
            '--preFitValue',str(seed['r']),
            '--cminDefaultMinimizerStrategy','2','--cminDefaultMinimizerTolerance','0.00001',
            '--cminDefaultMinimizerPrecision','1e-12','--X-rtd','MINIMIZER_freezeDisassociatedParams','-v','1']
        receipt=J/'receipt.json'
        if not receipt.exists():
            if (J/'process.json').exists():raise RuntimeError('Unfinished previous reference attempt; no automatic retry')
            dump(J/'command.json',cmd);start=time.monotonic()
            with (J/'combine.log').open('w') as log:
                proc=subprocess.Popen(cmd,cwd=J,stdout=log,stderr=subprocess.STDOUT,start_new_session=True,
                    env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1'))
                dump(J/'process.json',dict(pid=proc.pid,started_at=stamp()))
                try:rc=proc.wait(timeout=180)
                except subprocess.TimeoutExpired:
                    os.killpg(proc.pid,signal.SIGTERM)
                    try:proc.wait(timeout=5)
                    except subprocess.TimeoutExpired:os.killpg(proc.pid,signal.SIGKILL);proc.wait()
                    rc=124
            dump(receipt,dict(exit_code=rc,seconds=time.monotonic()-start,finished_at=stamp(),timeout_seconds=180))
        target=J/'fitDiagnostics_profile_reference.root'
        q,_=read_stage(p['workspace'],target,'s')
        if q.get('readable'):
            f=ROOT.TFile.Open(str(target));fit=f.Get('fit_s');rv=fit.floatParsFinal().find('r')
            assert rv and not rv.isConstant(),'Unconditional refinement must float r'
            lower=abs(rv.getVal()-rv.getMin())<=1e-4*max(1.,rv.getMax()-rv.getMin())
            upper=abs(rv.getVal()-rv.getMax())<=1e-4*max(1.,rv.getMax()-rv.getMin())
            q['status_history']=[dict(label=fit.statusLabelHistory(i),code=int(fit.statusCodeHistory(i))) for i in range(fit.numStatusHistory())]
            q['num_invalid_NLL']=int(fit.numInvalidNLL());f.Close()
            q['physical_lower_r_boundary']=bool(lower);q['upper_r_boundary']=bool(upper)
            q['point_valid']=bool(q['status']==0 and q['edm_valid'] and q['minimum']>0 and q['minimum_bin']>0 and math.isfinite(q['exact_deviance']) and q['r']>=0 and not upper and q['continuity_passed'])
            q['numerical_ok']=q['point_valid'] and q['covQual']==3
        q['receipt']=read(receipt)
        if q['receipt']['exit_code']!=0:q['point_valid']=False;q['numerical_ok']=False
        dump(J/'result.json',q)
        assert sha(p['workspace'])==p['workspace_sha256']
        profile=read(D/'profile.json')
        original=D/'profile_before_reference_refinement.json'
        if not original.exists():dump(original,profile)
        references=[x for x in profile['references'] if x.get('kind')!='refined_unconditional']
        if q.get('point_valid'):
            references.append(dict(exact_deviance=q['exact_deviance'],r=q['r'],source=q['source'],kind='refined_unconditional'))
        best=min(references,key=lambda x:x['exact_deviance'])
        profile.update(references=references,reference_min_deviance=best['exact_deviance'],reference_source=best['source'],reference_r=best['r'],
            reference_description='Lowest valid exact deviance among previous physical fits, fixed-r points, and one new unconditional fit from the deepest conditional point; global minimum not guaranteed.',
            original_reference_improvement=p['previous_physical_reference']['exact_deviance']-best['exact_deviance'],
            reference_refinement=dict(source=str(J/'result.json'),attempted=True,point_valid=q.get('point_valid',False),native_status=q.get('status'),covQual=q.get('covQual'),edm=q.get('edm')),
            updated_at=stamp())
        for point in profile['points']:
            point['q']=point['exact_deviance']-best['exact_deviance'] if point.get('exact_deviance') is not None else None
        dump(D/'profile.json',profile)
        out=dict(side=side,mass=mass,status='completed',point_valid=q.get('point_valid',False),
            new_r=q.get('r'),new_deviance=q.get('exact_deviance'),chosen_reference=best,finished_at=stamp())
        dump(J/'completion.json',out);print(side,mass,out,flush=True);return out


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--run',action='store_true');parser.add_argument('--workers',type=int,default=2)
    parser.add_argument('--side',choices=['low','high']);parser.add_argument('--mass',type=float);args=parser.parse_args()
    assert 1<=args.workers<=4
    cases=[(s,m) for s in ['low','high'] for m in MASSES if (not args.side or s==args.side) and (args.mass is None or m==args.mass)]
    if not args.run:
        for s,m in cases:
            try:print(json.dumps(plan(s,m)))
            except Exception as error:print(s,m,str(error))
    else:
        state=read(R/'runner.json');assert state['status']=='completed','Do not run while main profiles are active'
        reverse=R/'reverse_runner.json'
        if reverse.exists():assert read(reverse)['status']=='completed','Do not run while reverse continuation is active'
        dump(R/'reference_refinement_runner.json',dict(status='running',started_at=stamp(),pid=os.getpid(),cases=cases,workers=args.workers))
        with concurrent.futures.ProcessPoolExecutor(max_workers=args.workers) as pool:
            results=list(pool.map(run_one,[x[0] for x in cases],[x[1] for x in cases]))
        dump(R/'reference_refinement_runner.json',dict(status='completed',finished_at=stamp(),results=results))
