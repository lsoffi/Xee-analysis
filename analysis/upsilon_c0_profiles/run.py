"""Native conditional likelihood scans, preserving the previous C0 analysis."""
from pathlib import Path
import argparse, concurrent.futures, datetime, hashlib, json, math, os, signal, subprocess, sys, time

P=Path(__file__).resolve().parents[2]
S=P/'analysis/upsilon_c0_profiles'
R=P/'runs/20261102_upsilon_c0_profiles'
OLD=P/'runs/20261031_separate_sr_c0_limits'
MASSES=[9.4,9.5,10.0,10.1]

def stamp():return datetime.datetime.now(datetime.timezone.utc).isoformat()
def read(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def dump(p,q):
    p=Path(p);p.parent.mkdir(parents=True,exist_ok=True)
    t=p.with_suffix(p.suffix+'.tmp.'+str(os.getpid()));t.write_text(json.dumps(q,indent=2,allow_nan=False)+'\n');t.replace(p)

def plan(side,mass):
    old=OLD/side/f'm{mass:.1f}';fd=read(old/'fit_results.json');lim=read(old/'results.json')
    observed=lim['limits']['observed'] if lim['observed_limit_check']['verified'] else 0
    scale=max(50.,1.5*max(observed,lim['limits']['expected_50']),max(0.,fd['r'])+6*fd['r_error'])
    rmax=math.ceil(scale/25)*25
    grid=sorted(set([round(rmax*(i/24)**2,9) for i in range(25)]+[max(0.,fd['r'])]))
    return dict(side=side,mass=mass,workspace=str(old/'workspace.root'),workspace_sha256=sha(old/'workspace.root'),
        source_fit=str(old/'fit_results.json'),source_fit_sha256=sha(old/'fit_results.json'),
        grid=grid,rmax=rmax,source_limits=str(old/'results.json'),
        purpose='Observed conditional likelihood profile; no confidence-level calibration or new limits.',
        fixed_parameter='r',floating='All other original floating parameters, including eight Upsilon means and normalizations.',
        seeds=['original_B','original_SB','previous_valid_conditional'],
        conditional_method='FitDiagnostics --customStartingPoint --skipSBFit: stored fit_b minimizes full ModelConfig at supplied fixed r.')

def native(old,side,mass,r,params,directory,rmax,tolerance='0.00001',strategy='2'):
    directory.mkdir(parents=True,exist_ok=True)
    output=directory/'result.json'
    if output.exists():return read(output)
    initial=dict(params,r=r)
    arguments=','.join(f'{n}={v:.17g}' for n,v in sorted(initial.items()))
    cmd=[str(P/'bin/combine'),'-M','FitDiagnostics','-d',str(old/'workspace.root'),'-m',f'{mass:.1f}',
        '-n','_conditional','--customStartingPoint','--skipSBFit','--keepFailures','--robustFit','0','--minos','none',
        '--rMin','0','--rMax',str(max(100.,rmax*1.1)),'--setParameters',arguments,
        '--cminDefaultMinimizerStrategy',strategy,'--cminDefaultMinimizerTolerance',tolerance,
        '--cminDefaultMinimizerPrecision','1e-12','--X-rtd','MINIMIZER_freezeDisassociatedParams','-v','1']
    receipt=directory/'receipt.json'
    if not receipt.exists():
        if (directory/'process.json').exists():raise RuntimeError('Previous attempt has no closed receipt: '+str(directory))
        dump(directory/'command.json',cmd);started=time.monotonic()
        with (directory/'combine.log').open('w') as log:
            p=subprocess.Popen(cmd,cwd=directory,stdout=log,stderr=subprocess.STDOUT,start_new_session=True,
                env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1'))
            dump(directory/'process.json',dict(pid=p.pid,started_at=stamp()))
            try:rc=p.wait(timeout=180)
            except subprocess.TimeoutExpired:
                os.killpg(p.pid,signal.SIGTERM)
                try:p.wait(timeout=5)
                except subprocess.TimeoutExpired:os.killpg(p.pid,signal.SIGKILL);p.wait()
                rc=124
        dump(receipt,dict(exit_code=rc,seconds=time.monotonic()-started,finished_at=stamp()))
    from extract import extract
    q=extract(old/'workspace.root',directory/'fitDiagnostics_conditional.root',r)
    q.update(receipt=read(receipt),anchor=directory.name)
    if q['receipt']['exit_code']!=0:
        q['point_valid']=False;q['flags'].append('native_nonzero_exit')
    # The extractor's writer handles native infinite parameter bounds.
    from fit_helpers import dump as result_dump
    result_dump(output,q)
    return read(output)

def summarize(config,selected,complete=False,pilot=False):
    old=OLD/config['side']/f'm{config["mass"]:.1f}';fd=read(old/'fit_results.json')
    references=[dict(exact_deviance=fd['fit']['b']['exact_deviance'],r=0.,source=fd['fit']['b']['source'],kind='previous_B')]
    if fd['point_valid'] and fd['r']>=0:
        references.append(dict(exact_deviance=fd['fit']['s']['exact_deviance'],r=fd['r'],source=fd['fit']['s']['source'],kind='previous_SB'))
    references += [dict(exact_deviance=q['exact_deviance'],r=q['r'],source=q['source'],kind='conditional') for q in selected if q.get('point_valid')]
    best=min(references,key=lambda q:q['exact_deviance'])
    points=[]
    for q in selected:
        x={k:v for k,v in q.items() if k not in ['parameters','initialization','channels']}
        x['yields']=q.get('selected_yields');x['q']=q['exact_deviance']-best['exact_deviance'] if q.get('exact_deviance') is not None else None
        points.append(x)
    result=dict(side=config['side'],mass=config['mass'],status='complete' if complete else 'running',pilot=pilot,
        rmax=config['rmax'],planned_points=len(config['grid']),points=points,
        reference_min_deviance=best['exact_deviance'],reference_source=best['source'],reference_r=best['r'],references=references,
        reference_description='Lowest valid absolute extended deviance among previous physical B/S fits and conditional points; no statistical calibration.',
        source_workspace=config['workspace'],source_workspace_sha256=config['workspace_sha256'],updated_at=stamp(),
        original_reference_improvement=min(x['exact_deviance'] for x in references if x['kind'].startswith('previous_'))-best['exact_deviance'])
    from fit_helpers import dump as result_dump
    result_dump(R/config['side']/f'm{config["mass"]:.1f}'/('pilot_profile.json' if pilot else 'profile.json'),result)
    return result

def run_case(side,mass,pilot=False):
    import sys
    sys.path.insert(0,str(S))
    # Load ROOT once per isolated worker, not concurrently in Python threads.
    import extract
    config=plan(side,mass);D=R/side/f'm{mass:.1f}';D.mkdir(parents=True,exist_ok=True)
    cpath=D/'configuration.json'
    if cpath.exists():assert read(cpath)==config,'Immutable profile configuration changed'
    else:dump(cpath,config)
    old=OLD/side/f'm{mass:.1f}';fd=read(old/'fit_results.json')
    seeds={tag:{p['name']:p['value'] for p in fd['fit'][stage]['parameters']} for tag,stage in [('original_B','b'),('original_SB','s')]}
    grid=sorted(set([0.,max(0.,fd['r'])])) if pilot else config['grid']
    selected=[];previous=None
    for r in grid:
        pd=D/('r_'+format(r,'.12g').replace('.','p'));trial_seeds=dict(seeds)
        if previous and not pilot:trial_seeds['previous_valid_conditional']=previous['initialization']
        trials=[native(old,side,mass,r,pars,pd/tag,config['rmax']) for tag,pars in trial_seeds.items()]
        if r==0:
            # The previous audited B-only fit is exactly this conditional point.
            # Reuse its native RooFitResult, retaining all new failed attempts.
            existing=extract.extract(old/'workspace.root',fd['fit']['b']['source'],0.)
            existing.update(anchor='audited_existing_B',receipt=fd['fit']['b']['receipt'])
            from fit_helpers import dump as result_dump
            result_dump(pd/'audited_existing_B/result.json',existing)
            trials.append(read(pd/'audited_existing_B/result.json'))
        valid=[q for q in trials if q.get('point_valid')]
        finite=[q for q in trials if q.get('exact_deviance') is not None and math.isfinite(q['exact_deviance'])]
        if not finite:
            q=dict(r=r,exact_deviance=None,point_valid=False,status=None,edm=None,covQual=None,flags=['all_anchors_unreadable'],source=None,anchor=None)
        else:
            q=dict(min(valid or finite,key=lambda q:q['exact_deviance']))
            q['trials']=[dict(source=t.get('source'),point_valid=t.get('point_valid'),exact_deviance=t.get('exact_deviance'),flags=t['flags'],anchor=t['anchor']) for t in trials]
        if valid:previous=q
        selected.append(q);summarize(config,selected,pilot=pilot)
        print(side,f'{mass:.1f}',f'r={r:.6g}',q.get('exact_deviance'),q['point_valid'],q['flags'],flush=True)
    assert sha(old/'workspace.root')==config['workspace_sha256']
    assert sha(old/'fit_results.json')==config['source_fit_sha256']
    q=summarize(config,selected,complete=True,pilot=pilot)
    dump(D/('pilot_completion.json' if pilot else 'completion.json'),dict(completed=True,points=len(selected),valid=sum(p.get('point_valid',False) for p in selected),finished_at=stamp()))
    return dict(side=side,mass=mass,points=len(selected),valid=sum(p.get('point_valid',False) for p in selected))

if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('--workers',type=int,default=4);a.add_argument('--pilot',action='store_true');a.add_argument('--central',action='store_true');args=a.parse_args()
    assert 1<=args.workers<=6
    cases=[('low',9.5)] if args.pilot else [(side,m) for m in ([9.5,10.] if args.central else [9.5,10.,9.4,10.1]) for side in ['low','high']]
    R.mkdir(parents=True,exist_ok=True)
    name='pilot_runner' if args.pilot else 'central_runner' if args.central else 'runner'
    dump(R/(name+'.json'),dict(status='running',pid=os.getpid(),started_at=stamp(),cases=cases,workers=args.workers))
    with concurrent.futures.ProcessPoolExecutor(max_workers=args.workers) as pool:
        futures=[pool.submit(run_case,s,m,args.pilot) for s,m in cases]
        out=[f.result() for f in futures]
    dump(R/(name+'.json'),dict(status='completed',finished_at=stamp(),results=out))
