from pathlib import Path
import sys,json
P=Path.cwd();sys.path.insert(0,str(P/'scripts/sr_bernstein_selection'))
from fit import *
R=P/'runs/20261030_sr_bernstein7_bias';G=json.loads((R/'generators.json').read_text())
def job(a):
 k,ei=a;g=next(f for f in G if f['era']==eras[ei]);rng=np.random.default_rng(819221+k*10+ei);r=fit((7,ei,rng.poisson(g['mu']),g['p']));r['toy']=k;return r
if __name__=='__main__':
 with ProcessPoolExecutor(max_workers=2) as pool,(R/'gof.jsonl').open('w') as out:
  for j,r in enumerate(pool.map(job,[(k,ei) for k in range(100) for ei in range(4)])):
   out.write(json.dumps(r)+'\n');out.flush()
   if j%40==39:print(j+1,'/400',flush=True)
 print('DONE',flush=True)
