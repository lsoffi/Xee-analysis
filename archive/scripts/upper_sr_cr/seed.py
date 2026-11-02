from pathlib import Path
import sys, json, ROOT
ROOT.gSystem.Load('libHiggsAnalysisCombinedLimit')
P=Path.cwd();degree=int(sys.argv[1]);D=P/f'runs/20261030_upper_sr_cr/degree{degree}/joint'
f=ROOT.TFile.Open(str(D/'workspace.root'));w=f.Get('w')
# Old data fits provide only initial values, never constraints or fixed estimates.
q=ROOT.TFile.Open(str(P/'runs/20261030_sr6_cr4_scan/m7.5/joint/fitdiagnostics/fitDiagnostics_local.root'))
fit=q.Get('fit_b')
for v in fit.floatParsFinal():
    name=v.GetName()
    if ('_CR' in name or name.startswith('Ymean_') or 'bkgScale_upsilon' in name) and w.var(name):
        w.var(name).setVal(v.getVal())
for i,era in enumerate(['2022','2022EE','2023','2023BPix'],1):
    for j in range(1,degree+1):w.var(f'local_c{j}_{era}').setVal(-.5 if j==1 else .05 if j==2 else 0)
    n=w.data('data_obs').reduce(ROOT.RooFit.Cut(f'CMS_channel==CMS_channel::SR_ch{i}')).sumEntries()
    y=sum(w.function(f'local_rate_{p}_{era}').getVal() for p in ['upsilon1s','upsilon2s'])
    w.var('local_yield_'+era).setVal(max(1,n-y))
w.var('r').setVal(0)
w.writeToFile(str(D/'workspace_seeded.root'))
f.Close();q.Close()
