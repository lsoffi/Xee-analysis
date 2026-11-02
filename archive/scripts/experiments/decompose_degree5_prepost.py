import numpy as np,json
from pathlib import Path
rows=json.load(open('outputs/degree5_allFloat_all_eras_parameters.json'));J=np.array([0.,-2.,0.,-.4,0.]);out=[]
for era in ['2022','2022EE','2023','2023BPix']:
 pars=[next(a for a in rows if a['era']=='combined' and a['fit']=='fit_s' and a['parameter']==f't{i}_{era}_envelope') for i in range(5)];pre=np.array([a['prefit'] for a in pars]);post=np.array([a['postfit'] for a in pars]);ipre=6+pre@J;ipost=6+post@J
 rec={'era':era,'parameters':pars,'points':[]}
 for mass in [9,11,11.5,11.9]:
  T=np.polynomial.chebyshev.chebvander(np.array([(mass-9)/3]),5)[0,1:];Ppre=1+pre@T;Ppost=1+post@T;fp=Ppre/ipre;fq=Ppost/ipost;terms=(post-pre)*(T*ipre-Ppre*J)/(ipre*ipost);assert abs(terms.sum()-(fq-fp))<1e-12
  rec['points'].append(dict(mass=mass,shape_ratio=float(fq/fp),contributions_relative_to_prefit={f't{i}':float(terms[i]/fp) for i in range(5)}))
 print(era,'pre',pre,'post',post)
 for v in rec['points']:print(v)
 out.append(rec)
Path('outputs/degree5_prepost_exact_decomposition.json').write_text(json.dumps(out,indent=2))
