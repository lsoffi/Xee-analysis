from run import *
N=50;WINS=[0,5,10]
def job(arg):
 kind,k,wi,ei=arg;F=json.loads((R/'fits.json').read_text());o3=next(f for f in F if f['idx']==2*wi and f['era']==eras[ei]);o4=next(f for f in F if f['idx']==2*wi+1 and f['era']==eras[ei]);li=round((o3['lo']-6)*75);ui=round((o3['hi']-6)*75)
 if kind in ['bias4','bias5','asimov4','asimov5']:
  deg=int(kind[-1]);g=json.loads(((old/'bernstein_generator.json') if deg==4 else (R/'generator5.json')).read_text());rng=np.random.default_rng((987120 if deg==4 else 1487120)+k);dd=[rng.poisson(f['mu']) for f in g['fits']] if 'bias' in kind else [np.array(f['mu']) for f in g['fits']];data=dd[ei][li-232:ui-232]
 else:
  ob=o3 if kind=='gof3' else o4;rng=np.random.default_rng(2092120+k*100+wi*4+ei+(0 if kind=='gof3' else 100000));data=rng.poisson(ob['mu'])
 out=[];r3=None
 if kind!='gof4':
  r3=fit((2*wi,ei,data,o3['p']));r3.update(kind=kind,toy=k);out.append(r3)
 seed=o4['p'] if r3 is None else r3['p'][:4]+[0.]+r3['p'][4:]
 r4=fit((2*wi+1,ei,data,seed));r4.update(kind=kind,toy=k);out.append(r4);return out
if __name__=='__main__':
 jobs=[(kind,k,wi,ei) for kind in ['gof3','gof4','bias4','bias5'] for k in range(N) for wi in WINS for ei in range(4)]+[(kind,0,wi,ei) for kind in ['asimov4','asimov5'] for wi in WINS for ei in range(4)]
 with ProcessPoolExecutor(max_workers=3) as pool,(R/'pilot.jsonl').open('w') as out:
  for j,rs in enumerate(pool.map(job,jobs)):
   for r in rs:out.write(json.dumps(r)+'\n')
   out.flush()
   if j%60==59:print(j+1,'/',len(jobs),flush=True)
 print('DONE',flush=True)
