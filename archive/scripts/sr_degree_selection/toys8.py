from fit import *
def job(a):
 gen,k,ei,data,seed=a;r=fit((8,ei,data,seed));r.update(generator_degree=gen,toy=k);return r
if __name__=='__main__':
 G=json.loads((R/'generators.json').read_text());old=[json.loads(l) for l in (R/'toys.jsonl').read_text().splitlines()];jobs=[]
 for gen in [7,8]:
  for k in range(100):
   for ei,e in enumerate(eras):
    if gen==7:
     a=next(f for f in old if f['degree']==7 and f['generator_degree']==7 and f['toy']==k and f['era']==e);data=a['counts'];seed=a['p'][:8]+[0.]+a['p'][8:]
    else:
     ob=next(f for f in G if f['degree']==8 and f['era']==e);rng=np.random.default_rng(917110+k*100+8*4+ei);data=rng.poisson(ob['mu']);seed=ob['p']
    jobs.append((gen,k,ei,data,seed))
 with ProcessPoolExecutor(max_workers=3) as pool,(R/'toys8.jsonl').open('w') as out:
  for j,r in enumerate(pool.map(job,jobs)):
   out.write(json.dumps(r)+'\n');out.flush()
   if j%80==79:print(j+1,'/800',flush=True)
 print('DONE',flush=True)
