"""One bounded retry of conditional points with demonstrated lower feasible candidates."""
from pathlib import Path
import concurrent.futures, math
from run import R,OLD,native,read,dump,stamp,summarize

def job(issue):
    import extract
    side,mass,r=issue['side'],issue['mass'],issue['r'];D=R/side/f'm{mass:.1f}';c=read(D/'configuration.json')
    source=read(Path(issue['bestsource']).parent/'result.json')
    pd=D/('r_'+format(r,'.12g').replace('.','p'))
    q=native(OLD/side/f'm{mass:.1f}',side,mass,r,source['initialization'],pd/'lower_candidate_refit',c['rmax'],tolerance='0.001',strategy='1')
    return dict(side=side,mass=mass,r=r,status=q.get('status'),point_valid=q.get('point_valid'),D=q.get('exact_deviance'))

if __name__=='__main__':
    issues=read(R/'branch_check.json')['issues'];dump(R/'branch_runner.json',dict(status='running',started_at=stamp(),issues=issues))
    with concurrent.futures.ProcessPoolExecutor(max_workers=4) as pool:results=list(pool.map(job,issues))
    for profile_path in R.glob('*/m*/profile.json'):
        from extract import extract
        from fit_helpers import dump as safe_dump
        oldprofile=read(profile_path);D=profile_path.parent;c=read(D/'configuration.json')
        if not (D/'profile_before_branch_refit.json').exists():dump(D/'profile_before_branch_refit.json',oldprofile)
        points=[]
        for old in oldprofile['points']:
            pd=D/('r_'+format(old['r'],'.12g').replace('.','p'))
            trials=[read(p) for p in pd.glob('*/result.json')]
            valid=[q for q in trials if q.get('point_valid')]
            feasible=[q for q in trials if q.get('exact_deviance') is not None and math.isfinite(q['exact_deviance']) and q.get('minimum_bin',0)>0 and min(q.get('continuum_minima',{'x':-1}).values())>0 and q.get('continuity_passed')]
            q=dict(min(valid or feasible,key=lambda q:q['exact_deviance']))
            q['lower_feasible_candidate_gap']=q['exact_deviance']-min(x['exact_deviance'] for x in feasible)
            q['minimum_check_passed']=q['lower_feasible_candidate_gap']<=.02
            q['profile_point_valid']=q['point_valid'] and q['minimum_check_passed']
            if not q['minimum_check_passed']:q['flags']=q['flags']+['unresolved_lower_deviance_candidate']
            q['trials']=[dict(source=t.get('source'),point_valid=t.get('point_valid'),exact_deviance=t.get('exact_deviance'),flags=t['flags'],anchor=t['anchor']) for t in trials]
            points.append(q)
        profile=summarize(c,points,complete=True)
        profile['minimum_check_tolerance']=.02
        for k in ['reverse_pass','reference_refinement']:
            if k in oldprofile:profile[k]=oldprofile[k]
        safe_dump(profile_path,profile)
    dump(R/'branch_runner.json',dict(status='completed',finished_at=stamp(),results=results,one_attempt_per_issue=True))
