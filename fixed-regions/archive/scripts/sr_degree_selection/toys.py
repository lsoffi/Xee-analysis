from fit import *
N=100
OBS=json.loads((R/'observed.json').read_text())
def job(arg):
 deg,k,ei=arg;ob=next(f for f in OBS if f['degree']==deg and f['era']==eras[ei]);rng=np.random.default_rng(917110+k*100+deg*4+ei);data=rng.poisson(ob['mu']);a=fit((deg,ei,data,ob['p']));a.update(generator_degree=deg,toy=k);out=[a]
 if deg<7:
  seed=a['p'][:deg+1]+[0.]+a['p'][deg+1:];b=fit((deg+1,ei,data,seed));b.update(generator_degree=deg,toy=k);out.append(b)
 return out
if __name__=='__main__':
 jobs=[(deg,k,ei) for deg in range(3,8) for k in range(N) for ei in range(4)]
 with ProcessPoolExecutor(max_workers=3) as pool,(R/'toys.jsonl').open('w') as out:
  for j,rs in enumerate(pool.map(job,jobs)):
   for r in rs:out.write(json.dumps(r)+'\n')
   out.flush()
   if j%80==79:print(j+1,'/',len(jobs),flush=True)
 print('DONE',flush=True)
