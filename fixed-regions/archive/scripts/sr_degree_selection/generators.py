from fit import *
if __name__=='__main__':
 obs=json.loads((R/'observed.json').read_text());out=[f for f in obs if f['degree']==7]
 with ProcessPoolExecutor(max_workers=3) as pool:
  jobs=[]
  for ei,e in enumerate(eras):
   f=next(f for f in out if f['era']==e);seed=f['p'][:8]+[0.]+f['p'][8:];jobs.append((8,ei,None,seed))
  for f in pool.map(fit,jobs):out.append(f);print('generator8',f['era'],f['valid'],f['D'],flush=True)
 (R/'generators.json').write_text(json.dumps(out))
 anchors=[]
 for gen in [7,8]:
  for deg in [4,5,6,7]:
   jobs=[]
   for ei,e in enumerate(eras):
    ob=next(f for f in obs if f['degree']==deg and f['era']==e);g=next(f for f in out if f['degree']==gen and f['era']==e);jobs.append((deg,ei,g['mu'],ob['p']))
   with ProcessPoolExecutor(max_workers=3) as pool:
    for f in pool.map(fit,jobs):f['generator_degree']=gen;anchors.append(f);print('anchor',gen,deg,f['era'],f['valid'],flush=True)
 (R/'asimov_anchors.json').write_text(json.dumps(anchors))
