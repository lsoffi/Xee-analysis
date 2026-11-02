from pathlib import Path
import ROOT,json,numpy as np,hashlib
ROOT.gSystem.Load('libHiggsAnalysisCombinedLimit');P=Path.cwd();R=P/'runs/20261029_upsilon_control_cheb3';q=json.loads((R/'parameters.json').read_text());out={}
for mode in ['control','joint']:
 f=ROOT.TFile.Open(str(R/mode/'workspace.root'));w=f.Get('w');checks=[]
 for stage,fit in q[mode].items():
  for a in fit['parameters']:
   if w.var(a['name']):w.var(a['name']).setVal(a['value'])
  for suffix in (['','_CR'] if mode=='joint' else ['']):
   for era in ['2022','2022EE','2023','2023BPix']:
    vs=[w.var(f'local_c{i}_{era}{suffix}') for i in [1,2,3]];c=[1]+[v.getVal() for v in vs];roots=np.polynomial.chebyshev.chebroots(np.polynomial.chebyshev.chebder(c));points=[-1,1]+[float(z.real) for z in roots if abs(z.imag)<1e-10 and -1<z.real<1];minimum=min(np.polynomial.chebyshev.chebval(points,c));assert minimum>0
    assert all(v.getMin()+1e-5<v.getVal()<v.getMax()-1e-5 for v in vs)
    checks.append(dict(stage=stage,era=era,region=suffix or 'SR' if mode=='joint' else 'CR',coefficients=c,minimum=float(minimum)))
 out[mode]=checks
 f.Close()
for part in ['control','signal']:
 a=json.loads((R/part/'audit.json').read_text());assert hashlib.sha256(Path(a['source']).read_bytes()).hexdigest()==a['source_sha256']
(R/'validation.json').write_text(json.dumps(out,indent=2));print('Verified exact extrema, coefficient boundaries and source hashes')
