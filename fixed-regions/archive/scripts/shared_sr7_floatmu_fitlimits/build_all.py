"""Prepare the 61 degree-seven workspaces; never launches statistical fits."""
from pathlib import Path
import hashlib
import json
import time
import model


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    rows=[];start=time.time()
    masses=['7.5','9.0']+[f'{k/10:.1f}' for k in range(60,121) if k not in [75,90]]
    for mass in masses:
        output=model.R/f'observed/m{mass}/joint/workspace_seeded.root'
        if not output.exists():
            model.build(mass,degree=7,family='envelope',output=output,freeze_means=False)
        audit=json.loads(output.with_suffix('.audit.json').read_text())
        assert audit['degree']==7 and not audit['means_fixed']
        assert sha(audit['source'])==audit['source_sha256']
        rows.append(dict(mass=float(mass),workspace=str(output),workspace_sha256=sha(output),
                         source=audit['source'],source_sha256=audit['source_sha256']))
        print('READY',mass,'elapsed',round(time.time()-start,1),flush=True)
    path=model.R/'workspace_manifest.json'
    if path.exists():
        previous=json.loads(path.read_text())
        old={q['mass']:q['workspace_sha256'] for q in previous['workspaces']}
        assert old=={q['mass']:q['workspace_sha256'] for q in rows}
        print('Existing immutable manifest verified',flush=True)
        return
    manifest=dict(workspaces=rows,count=len(rows),builder=str(Path(model.__file__).resolve()),
                  builder_sha256=sha(model.__file__),source_preservation_verified=True,
                  degree=7,floating_means=True,seconds=time.time()-start)
    path.write_text(json.dumps(manifest,indent=2)+'\n')
    print('ALL61DONE',round(time.time()-start,1),flush=True)


if __name__=='__main__':main()
