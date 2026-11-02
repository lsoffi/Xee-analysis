from study import *
s=json.loads((R/'scan.json').read_text());base=s[2]
for r in s:
 if r['degree']==3:print(round(r['lo'],3),round(r['hi'],3),'p',round(r['p_asym'],3),'tail',np.round(r['prediction'],1),'err',np.round(r['error'],1),'shift',np.round((np.array(r['prediction'])-base['prediction'])/base['error'],2),'bounds',r['boundary_count'])
g=combined(6+232/75,6+338/75,4);(R/'generator4.json').write_text(json.dumps(g));print('GEN4',g['valid'],g['p_asym']);
for f in g['fits']:
 # polynomial p(t) = Bernstein in u=(t+1)/2
 from numpy.polynomial import Polynomial
 from math import comb
 q=C.Chebyshev(f['p'][:5]).convert(kind=Polynomial)(Polynomial([-1,2]));b=[sum(q.coef[k]*comb(j,k)/comb(4,k) for k in range(j+1)) for j in range(5)];print('Bernstein4 coefficients',b)
