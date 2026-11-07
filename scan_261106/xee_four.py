from pathlib import Path
import subprocess
old=Path('low_logpoly3_totalC0')
total=(old/'build.py').read_text()
tails=total[total.index('            tails=[]'):total.index('        else:\n            scale=',total.index('            tails=[]'))]
for degree in (2,3,6):
 p=Path(f'low_cheb{degree}_totalC0');p.mkdir();(p/'logs').mkdir()
 s=Path('build.py').read_text().replace('BASE=Path(__file__).resolve().parent','BASE=Path(__file__).resolve().parent.parent').replace('in (3,4,5,6)','in (2,3,4,5,6)').replace('seed=np.pad(seed,(0,size-len(seed)))','seed=np.pad(seed[:size],(0,max(0,size-len(seed))))')
 mark='        # Original signal functions'
 s=s.replace(mark,'''        if contacts:
            assert len(contacts)==1 and rank==1
            g,edge=contacts[0];d=norms[g];density=targets[0][1]
'''+tails+mark)
 (p/'build.py').write_text(s)
 for name in ('worker.sh','draw.py','jobs.sub'):
  data=(old/name).read_text().replace(old.name,p.name).replace('XEE_DEGREE=3',f'XEE_DEGREE={degree}').replace('exp(pol. grado 3)',f'Chebyshev {degree}')
  (p/name).write_text(data)
 (p/'worker.sh').chmod(0o755);(p/'masses.txt').write_text('0.6\n' if degree==6 else '0.4\n0.5\n0.6\n')
 subprocess.run(['condor_submit','jobs.sub'],cwd=p,check=True)
