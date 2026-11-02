from pathlib import Path
import ROOT,json,sys,numpy as np
ROOT.gSystem.Load('libHiggsAnalysisCombinedLimit');P=Path.cwd();R=P/'runs/20261030_sr6_cr4_scan';d=R/('m'+sys.argv[1])/'joint';f=ROOT.TFile.Open(str(d/'workspace.root'));w=f.Get('w');rows=json.loads((P/'runs/20261030_sr_degree_selection/observed.json').read_text())
for a in rows:
 if a['degree']!=6:continue
 e=a['era'];p=np.array(a['p']);assert a['valid']
 for suffix,co in [('',p[:7]),('_CR',p[7:12])]:
  for j in range(1,len(co)):w.var(f'local_c{j}_{e}{suffix}').setVal(co[j]/co[0])
  integral=np.polynomial.chebyshev.chebint(co);y=(np.polynomial.chebyshev.chebval(1,integral)-np.polynomial.chebyshev.chebval(-1,integral))/2;w.var(f'local_yield_{e}{suffix}').setVal(y)
 # Convert the same Chebyshev polynomial to positive Bernstein coefficients for an identical initial shape.
 from math import comb
 from numpy.polynomial import Polynomial,Chebyshev
 power=Chebyshev(p[:7]).convert(kind=Polynomial)(Polynomial([-1,2])).coef
 bs=np.array([sum(power[k]*comb(i,k)/comb(6,k) for k in range(i+1)) for i in range(7)])
 assert min(bs)>0,(e,bs)
 for j in range(1,7):w.var(f'local_b{j}_{e}').setVal(bs[j]/bs[0])
 for j,pr in enumerate(['upsilon1s','upsilon2s']):
  w.var(f'CMS_EXO25020_bkgScale_{pr}_{e}').setVal(p[12+j]);w.var(f'Ymean_{pr}_{e}').setVal(p[14+j])
w.var('r').setVal(.01);w.writeToFile(str(d/'workspace_seeded.root'));f.Close();(d/'workspace.root').rename(d/'workspace_unseeded.root');(d/'workspace_seeded.root').rename(d/'workspace.root')
