from pathlib import Path
import subprocess
base=Path.cwd()
stage=base/'early_fd'
stage.mkdir()
(stage/'logs').mkdir()
masses=[m for m in Path('masses.txt').read_text().split() if not (base/'results'/('m'+m)/'AsymptoticLimits.done').exists()]
(stage/'masses.txt').write_text('\n'.join(masses)+'\n')
s=Path('worker.sh').read_text().replace('out="$run/results/m$m"', '''if [[ -e "$run/results/m$m/FitDiagnostics.log" ]]; then exit 0; fi
out="$run/early_fd/results/m$m"''').replace('in Significance AsymptoticLimits FitDiagnostics','in FitDiagnostics')
s += ''
s=s.replace('exit "$rc"', '''if [[ "$rc" == 0 ]]; then
  cd "$run"
  python3 early_fd/draw_one.py "$m" "$out"
fi
exit "$rc"''')
(stage/'worker.sh').write_text(s);(stage/'worker.sh').chmod(0o755)
draw=Path('fit_gallery.py').read_text().split('for _ in range(720):')[0]
(stage/'draw_one.py').write_text(draw+'''
import sys,shutil,os
m=sys.argv[1];d=Path(sys.argv[2]);dest=out/('m'+m+'.png')
if not dest.exists():
 rec=next(r for r in records if r['mass']==float(m))
 local=Path(os.environ['_CONDOR_SCRATCH_DIR'])/'early.png'
 draw(rec,next(d.glob('fitDiagnostics*.root')),local)
 tmp=out/('.early_m'+m+'.png');shutil.copyfile(local,tmp)
 if not dest.exists():tmp.replace(dest)
 else:tmp.unlink()
''')
(stage/'jobs.sub').write_text(Path('jobs.sub').read_text())
print(len(masses),flush=True)
subprocess.run(['condor_submit','jobs.sub'],cwd=stage,check=True)
