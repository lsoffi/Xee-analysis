from pathlib import Path
import ROOT,json,sys,numpy as np
from numpy.polynomial import Polynomial,Chebyshev
from math import comb
ROOT.gSystem.Load('libHiggsAnalysisCombinedLimit');P=Path.cwd();mass=sys.argv[1];deg=int(sys.argv[2]);D=P/f'runs/20261030_upper_envelope_scan/degree{deg}/m{mass}/joint'
f=ROOT.TFile.Open(str(D/'workspace.root'));w=f.Get('w')
ref=json.loads((P/f'runs/20261030_upper_sr_cr_mu_fixed/degree{deg}/results.json').read_text())
for a in ref['parameters']:
    if w.var(a['name']):w.var(a['name']).setVal(a['value'])
means=json.loads((P/'runs/20261030_upper_sr_cr_mu_fixed/cr_only/results.json').read_text())['means']
for a in means:
    w.var(a['name']).setVal(a['value']);w.var(a['name']).setConstant(True)
for era in ['2022','2022EE','2023','2023BPix']:
    coeff=np.array([1]+[w.var(f'local_c{j}_{era}').getVal() for j in range(1,deg+1)])
    power=Chebyshev(coeff).convert(kind=Polynomial)(Polynomial([-1,2])).coef
    bs=np.array([sum(power[k]*comb(i,k)/comb(deg,k) for k in range(i+1)) for i in range(deg+1)])
    assert min(bs)>0,(era,bs)
    for j in range(1,deg+1):w.var(f'local_b{j}_{era}').setVal(bs[j]/bs[0])
    env=w.pdf('local_envelope_'+era);corr=[]
    for j in [0,1]:w.cat('pdfindex_'+era).setIndex(j);corr.append(env.getCorrection())
    assert abs(corr[0]-corr[1])<1e-8;w.cat('pdfindex_'+era).setIndex(0)
w.var('mass').setVal(11.2);w.var('r').setVal(0);w.var('r').setRange(-1000,1000)
w.writeToFile(str(D/'workspace_seeded.root'));f.Close()
