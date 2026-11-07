from study import *
from concurrent.futures import ProcessPoolExecutor
S=json.loads((R/'scan.json').read_text());base=S[2];G=json.loads((R/'bernstein_generator.json').read_text())
def calc_truth(g):
 return [sum(np.array(ER[i]['templates'])[:,(edges[:-1]>=a-1e-8)&(edges[1:]<=b+1e-8)].sum(axis=1)@np.array(g['fits'][i]['p'][-2:]) for i in range(4)) for a,b in [(6,9.093333333333334),(10.506666666666668,12)]]
def closure(i):
 rng=np.random.default_rng(400120+i);d=[rng.poisson(f['mu']) for f in base['fits']];r=combined(base['lo'],base['hi'],3,d,[f['p'] for f in base['fits']]);return dict(valid=r['valid'],prediction=r['prediction'],error=r['error'])
if __name__=='__main__':
 out={'truth_alternative':calc_truth(G),'truth_closure':calc_truth(base),'asimov':[]}
 for li,ui in [(232,338),(240,330),(240,323),(247,323),(255,315)]:
  for deg in [1,2,3]:
   obs=next(r for r in S if round((r['lo']-6)*75)==li and round((r['hi']-6)*75)==ui and r['degree']==deg);d=[np.array(f['mu'])[li-232:ui-232] for f in G['fits']];r=combined(obs['lo'],obs['hi'],deg,d,[f['p'] for f in obs['fits']]);out['asimov'].append({k:v for k,v in r.items() if k!='fits'})
 with ProcessPoolExecutor(max_workers=3) as pool:out['closure']=list(pool.map(closure,range(200)))
 good=[x for x in out['closure'] if x['valid']];p=np.array([x['prediction'] for x in good]);err=np.mean([x['error'] for x in good],axis=0);out['closure_summary']=dict(valid=len(good),total=200,bias_sigma=((p.mean(axis=0)-out['truth_closure'])/err).tolist(),MCerror_sigma=(p.std(axis=0,ddof=1)/np.sqrt(len(good))/err).tolist());(R/'checks.json').write_text(json.dumps(out,indent=2));print(out['closure_summary']);print('Asimov baseline',out['asimov'][2]['prediction'],'truth',out['truth_alternative'])
