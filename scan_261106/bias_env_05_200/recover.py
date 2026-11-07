from pathlib import Path
import json,subprocess
p=Path.cwd()
jobs=[line for n in ('pilot','remaining') for line in Path(n+'.txt').read_text().splitlines()]
bad=[]
for line in jobs:
 seed=line.split()[-1];f=Path('results')/seed/'result.json'
 if not f.exists() or not all(x['ok'] for x in json.loads(f.read_text())['branches'].values()):bad.append(line)
print('retry',len(bad),flush=True)
if bad:
 s=Path('score.py').read_text().replace('pars=fit.floatParsFinal()',"""if not fit:
  out[name]=dict(r=None,D=None,penalty=4*degree,score=None,status=-999,cov=-1,edm=None,ok=False)
  f.Close();fw.Close();continue
 pars=fit.floatParsFinal()""")
 Path('score_retry.py').write_text(s)
 s=Path('worker.sh').read_text().replace('out="$run/results/$seed"','out="$run/retry_results/$seed"').replace('MinimizerStrategy 0','MinimizerStrategy 1').replace('$run/score.py','$run/score_retry.py')
 Path('retry.sh').write_text(s);Path('retry.sh').chmod(0o755)
 Path('retry.txt').write_text('\n'.join(bad)+'\n');Path('retry.sub').write_text(Path('remaining.sub').read_text().replace('worker.sh','retry.sh').replace('remaining','retry').replace('logs/$(seed)','logs/retry_$(seed)'))
 subprocess.run(['condor_submit','retry.sub'],check=True)
