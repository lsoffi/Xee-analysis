from fit_stable import *
R=P/'runs/20261030_sr_bernstein7_bias';G=json.loads((R/'generators.json').read_text());A=json.loads((R/'asimov_anchors.json').read_text())
def job(args):
 toy,deg,ei=args;rng=np.random.default_rng(765421+toy);data=[rng.poisson(g['mu']) for g in G];assert [g['era'] for g in G]==eras;seed=next(f['p'] for f in A if f['degree']==deg and f['era']==eras[ei]);r=fit((deg,ei,data[ei],seed));r['toy']=toy;return r
if __name__=='__main__':
 with ProcessPoolExecutor(max_workers=6) as pool,(R/'toy_anchors.jsonl').open('w') as out:
  for j,r in enumerate(pool.map(job,[(t,d,ei) for t in range(30) for d in range(3,8) for ei in range(4)])):
   out.write(json.dumps(r)+'\n');out.flush()
   if j%40==39:print(j+1,'/600',flush=True)
 print('DONE',flush=True)
