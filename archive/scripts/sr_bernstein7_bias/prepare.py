from pathlib import Path
import sys,json
P=Path.cwd();sys.path.insert(0,str(P/'scripts/sr_degree_selection'))
from fit import *
R=P/'runs/20261030_sr_bernstein7_bias'
if __name__=='__main__':
 G=json.loads((P/'runs/20261030_sr_bernstein_selection/generator7.json').read_text());assert len(G)==4 and all(g['valid'] for g in G)
 for g in G:g['degree']='bern7'
 (R/'generators.json').write_text(json.dumps(G));obs=json.loads((P/'runs/20261030_sr_degree_selection/observed.json').read_text());out=[]
 for deg in range(3,8):
  jobs=[]
  for ei,e in enumerate(eras):
   g=next(f for f in G if f['era']==e);o=next(f for f in obs if f['era']==e and f['degree']==deg);jobs.append((deg,ei,g['mu'],o['p']))
  with ProcessPoolExecutor(max_workers=3) as pool:
   for f in pool.map(fit,jobs):f['generator_degree']='bern7';out.append(f);print(deg,f['era'],f['valid'],flush=True)
 (R/'asimov_anchors.json').write_text(json.dumps(out));(R/'signal_rates.json').write_text((P/'runs/20261030_sr_degree_selection/signal_rates.json').read_text())
