from toys import *
from scipy.stats import beta
G=json.loads((R/'bernstein_generator.json').read_text())
def one(i):
 rng=np.random.default_rng(1702120+i);d=[rng.poisson(f['mu']) for f in G['fits']]
 try:return bernstein_generator(d,False)['D']
 except Exception as e:return None
if __name__=='__main__':
 with ProcessPoolExecutor(max_workers=3) as pool:vals=list(pool.map(one,range(200)))
 good=[v for v in vals if v is not None];k=sum(v>=G['D'] for v in good);res=dict(observed=G['D'],valid=len(good),total=200,exceedances=k,p=k/len(good),deviances=vals);(R/'generator_gof.json').write_text(json.dumps(res));print({k:v for k,v in res.items() if k!='deviances'})
