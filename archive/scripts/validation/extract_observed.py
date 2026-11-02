import ROOT,pathlib,json,csv,math
ROOT.gSystem.Load('libHiggsAnalysisCombinedLimit');P=pathlib.Path('/Users/imb/Documents/Codex/Analisi-statistica-Xee/runs/20261029_background_validation');rows=[];pars=[];corr=[]
for m in ['7.5','9.5']:
 for model in ['original','degree5']:
  wd=P/f'observed/m{m}_{model}';f=ROOT.TFile.Open(str(wd/'fitDiagnostics_validation.root'));b=f.Get('fit_b');s=f.Get('fit_s');t=f.Get('tree_fit_sb');t.GetEntry(0);r=s.floatParsFinal().find('r');row=dict(mass=float(m),model=model,b_status=b.status(),s_status=s.status(),b_covQual=b.covQual(),s_covQual=s.covQual(),r=r.getVal(),r_error=r.getError(),r_low=float(t.r)-abs(float(t.rLoErr)),r_high=float(t.r)+abs(float(t.rHiErr)),q0=max(0,-2*float(t.nll_nll0)) if t.r>0 else 0)
  row['Z_local_from_FD']=math.sqrt(row['q0']);row['p0_local']=.5*math.erfc(row['Z_local_from_FD']/math.sqrt(2));rows.append(row)
  for stage,fit in [('B',b),('SB',s)]:
   for v in fit.floatParsFinal():
    pars.append(dict(mass=float(m),model=model,stage=stage,parameter=v.GetName(),value=v.getVal(),error=v.getError()))
   if stage=='SB':
    for v in fit.floatParsFinal():
     if v.GetName()!='r':corr.append(dict(mass=float(m),model=model,parameter=v.GetName(),correlation_r=fit.correlation('r',v.GetName())))
for name,rs in [('observed_comparison',rows),('observed_parameters',pars),('signal_background_correlations',corr)]:
 with (P/(name+'.csv')).open('w') as o:
  w=csv.DictWriter(o,fieldnames=list(rs[0]));w.writeheader();w.writerows(rs)
print(json.dumps(rows,indent=2));(P/'observed_comparison.json').write_text(json.dumps(rows,indent=2))
