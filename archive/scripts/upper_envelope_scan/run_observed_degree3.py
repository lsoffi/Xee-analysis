"""Degree-3 extension of the existing upper-SR scan; earlier degrees untouched."""
from pathlib import Path
import argparse,concurrent.futures,json
from run_observed import execute,job,P,PY,S,R

def prepare(mass):
    d=R/'degree3'/f'm{mass}'
    for region in ['signal','control']:
        execute([PY,str(S/'prepare.py'),mass,'3',region],d,'prepare_'+region,True)
    execute([PY,str(S/'build.py'),mass,'3'],d,'build',True)
    execute([PY,str(S/'seed.py'),mass,'3'],d,'seed',True)

if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--prepare-only',action='store_true')
    parser.add_argument('--masses',nargs='+')
    parser.add_argument('--workers',type=int,default=3)
    args=parser.parse_args()
    masses=args.masses or ['10.8']+[f'{m/10:.1f}' for m in range(105,121) if m!=108]
    assert 1<=args.workers<=3
    assert all(10.5<=float(m)<=12 for m in masses)
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as pool:
        if args.prepare_only:list(pool.map(prepare,masses))
        else:list(pool.map(job,[(3,m) for m in masses]))
