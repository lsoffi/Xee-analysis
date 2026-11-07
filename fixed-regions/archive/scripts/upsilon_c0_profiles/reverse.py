"""One uniform reverse continuation pass to check numerical branch stability."""
from pathlib import Path
import argparse,concurrent.futures,json,math
from run import P,S,R,OLD,native,read,dump,stamp,summarize

def job(side,mass):
    import extract
    from fit_helpers import dump as safe_dump
    D=R/side/f'm{mass:.1f}';config=read(D/'configuration.json');profile=read(D/'profile.json')
    saved=D/'profile_before_reverse.json'
    if not saved.exists():dump(saved,profile)
    points={float(x['r']):x for x in profile['points']}
    valid=[x for x in points.values() if x.get('point_valid')]
    start=max(valid,key=lambda x:x['r']);init=read(Path(start['source']).parent/'result.json')['initialization']
    for r in reversed(config['grid']):
        pd=D/('r_'+format(r,'.12g').replace('.','p'))
        candidate=native(OLD/side/f'm{mass:.1f}',side,mass,r,init,pd/'reverse_continuation',config['rmax'],tolerance='0.001')
        original=points[r]
        usable=[x for x in [original,candidate] if x.get('point_valid')]
        finite=[x for x in [original,candidate] if x.get('exact_deviance') is not None and math.isfinite(x['exact_deviance'])]
        chosen=dict(min(usable or finite,key=lambda x:x['exact_deviance'])) if finite else original
        if candidate.get('point_valid'):init=candidate['initialization']
        chosen['trials']=original.get('trials',[])+[dict(source=candidate.get('source'),anchor=candidate['anchor'],point_valid=candidate.get('point_valid'),exact_deviance=candidate.get('exact_deviance'),flags=candidate['flags'])]
        points[r]=chosen
        q=summarize(config,[points[x] for x in sorted(points)],complete=False)
        q['reverse_pass']=dict(status='running',tolerance=.001,reason='Uniform reverse continuation, same model, to audit failed points and competing minima.')
        safe_dump(D/'profile.json',q)
        print(side,mass,r,chosen.get('exact_deviance'),chosen.get('point_valid'),flush=True)
    q=summarize(config,[points[x] for x in sorted(points)],complete=True)
    q['reverse_pass']=dict(status='completed',tolerance=.001,reason='Uniform reverse continuation under unchanged model; all old attempts retained.')
    safe_dump(D/'profile.json',q)
    return dict(side=side,mass=mass,points=len(points),valid=sum(x.get('point_valid',False) for x in points.values()))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--workers',type=int,default=6);a=p.parse_args()
    assert read(R/'runner.json')['status']=='completed'
    cases=[(x['side'],x['mass']) for x in read(R/'manifest.json')['cases']]
    dump(R/'reverse_runner.json',dict(status='running',started_at=stamp(),workers=a.workers))
    with concurrent.futures.ProcessPoolExecutor(max_workers=a.workers) as pool:
        out=list(pool.map(job,[x[0] for x in cases],[x[1] for x in cases]))
    dump(R/'reverse_runner.json',dict(status='completed',finished_at=stamp(),results=out))
