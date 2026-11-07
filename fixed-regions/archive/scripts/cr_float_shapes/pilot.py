from run import *
N=50;IDS=[2,32]
def job(arg):
 kind,k,idx,ei=arg;F=json.loads((R/'fits.json').read_text());ob=next(f for f in F if f['idx']==idx and f['era']==eras[ei]);li=round((ob['lo']-6)*75);ui=round((ob['hi']-6)*75)
 if kind in ['bias','asimov']:
  g=json.loads((old/'bernstein_generator.json').read_text());rng=np.random.default_rng(987120+k);dd=[rng.poisson(f['mu']) for f in g['fits']] if kind=='bias' else [np.array(f['mu']) for f in g['fits']];d=dd[ei][li-232:ui-232]
 else:
  rng=np.random.default_rng(922120+k*100+idx*4+ei);d=rng.poisson(ob['mu'])
 r=fit((idx,ei,d,ob['p']));r.update(kind=kind,toy=k);return r
if __name__=='__main__':
 jobs=[(kind,k,idx,ei) for kind in ['gof','bias'] for k in range(N) for idx in IDS for ei in range(4)]+[('asimov',0,idx,ei) for idx in IDS for ei in range(4)]
 with ProcessPoolExecutor(max_workers=3) as pool,(R/'pilot.jsonl').open('w') as out:
  for j,r in enumerate(pool.map(job,jobs)):
   out.write(json.dumps(r)+'\n');out.flush()
   if j%40==39:print(j+1,'/',len(jobs),flush=True)
 print('DONE',flush=True)
