"""Build authorized separate-region models only; never launch statistical fits."""
from pathlib import Path
import argparse,json,time,sys,subprocess
from concurrent.futures import ThreadPoolExecutor,as_completed
import model

def main(include_cr=False,workers=4):
    points=[('low','7.5'),('high','11.0')]
    points += [('low',f'{k/10:.1f}') for k in range(60,91) if k!=75]
    points += [('high',f'{k/10:.1f}') for k in range(106,121) if k!=110]
    if include_cr:points += [(side,f'{k/10:.1f}') for k in range(91,106) for side in ['low','high']]
    rows=[];start=time.time()
    pending=[(side,mass) for side,mass in points if not (model.R/side/f'm{mass}/workspace.root').exists()]
    def prepare(point):
        side,mass=point
        result=subprocess.run([sys.executable,str(Path(model.__file__).resolve()),side,mass],capture_output=True,text=True)
        if result.returncode:raise RuntimeError(f'{side} {mass}: {result.stdout}\n{result.stderr}')
        return side,mass
    with ThreadPoolExecutor(max_workers=workers) as pool:
        futures=[pool.submit(prepare,p) for p in pending]
        for future in as_completed(futures):
            side,mass=future.result();print('BUILT',side,mass,round(time.time()-start,2),flush=True)
    for side,mass in points:
        path=model.R/side/f'm{mass}/workspace.root'
        assert path.exists()
        audit=json.loads(path.with_suffix('.audit.json').read_text())
        assert model.sha(path)==audit['workspace_sha256']
        assert model.sha(audit['source'])==audit['source_sha256']
        assert not audit['means_fixed'] and audit['discrete_categories']==0
        rows.append(dict(side=side,mass=float(mass),workspace=str(path),workspace_sha256=model.sha(path),source=audit['source'],source_sha256=audit['source_sha256']))
        print('READY',side,mass,round(time.time()-start,2),flush=True)
    manifest=dict(count=len(rows),includes_CR_mass_tests_in_both_models=include_cr,workspaces=rows,builder=str(Path(model.__file__).resolve()),builder_sha256=model.sha(model.__file__),seconds=time.time()-start)
    path=model.R/'workspace_manifest.json'
    if path.exists():
        old=json.loads(path.read_text());assert {(q['side'],q['mass']):q['workspace_sha256'] for q in old['workspaces']}=={(q['side'],q['mass']):q['workspace_sha256'] for q in rows}
    else:path.write_text(json.dumps(manifest,indent=2)+'\n')
    print('BUILD_COMPLETE',len(rows),round(time.time()-start,2),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--include-cr',action='store_true');p.add_argument('--workers',type=int,default=4);a=p.parse_args();main(a.include_cr,a.workers)
