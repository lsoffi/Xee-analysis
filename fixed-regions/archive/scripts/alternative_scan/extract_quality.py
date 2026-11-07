import ROOT,pathlib,json,csv,math
ROOT.gSystem.Load('libHiggsAnalysisCombinedLimit');P=pathlib.Path('/Users/imb/Documents/Codex/Analisi-statistica-Xee');R=P/'runs/20261029_alternative_background_scan';rows=[];pars=[];correlations=[]
for a in json.loads((R/'audit_scan.json').read_text()):
 wd=R/f"fits/{a['family']}/m{a['mass_label']}/fitdiagnostics";f=ROOT.TFile.Open(str(wd/'fitDiagnostics_alternative.root'))
 if not f or f.IsZombie():rows.append(dict(family=a['family'],mass=a['mass'],error='missing ROOT'));continue
 row=dict(family=a['family'],mass=a['mass'],mass_label=a['mass_label']);bounds=[]
 for name in ['fit_b','fit_s']:
  fit=f.Get(name)
  if not fit:row[name+'_missing']=True;continue
  row[name+'_status']=fit.status();row[name+'_covQual']=fit.covQual();row[name+'_edm']=fit.edm()
  for v in fit.floatParsFinal():
   value=v.getVal();lo=v.getMin();hi=v.getMax();at_bound=math.isfinite(lo+hi) and min(value-lo,hi-value)<1e-4*(hi-lo)
   pars.append(dict(family=a['family'],mass=a['mass'],stage=name,parameter=v.GetName(),value=value,error=v.getError(),min=lo,max=hi,at_bound=at_bound))
   if at_bound:bounds.append(dict(stage=name,parameter=v.GetName(),value=value,min=lo,max=hi))
   if name=='fit_s' and v.GetName()!='r':correlations.append(dict(family=a['family'],mass=a['mass'],parameter=v.GetName(),correlation_r=fit.correlation('r',v.GetName())))
  if name=='fit_s':
   r=fit.floatParsFinal().find('r');row['r_fit_result']=r.getVal();row['r_error_hesse']=r.getError()
 row['parameters_at_bounds']=bounds;rows.append(row);f.Close()
(R/'fit_quality_detailed.json').write_text(json.dumps(rows,indent=2))
for name,data in [('fit_quality_detailed',rows),('parameters_prepost_fits',pars),('correlations_with_r',correlations)]:
 keys=list(dict.fromkeys(k for r in data for k in r))
 with (R/(name+'.csv')).open('w') as f:w=csv.DictWriter(f,fieldnames=keys);w.writeheader();w.writerows(data)
print(json.dumps(dict(fits=len(rows),nonzero_status=sum(r.get('fit_b_status')!=0 or r.get('fit_s_status')!=0 for r in rows),covariance_warnings=sum(r.get('fit_b_covQual')!=3 or r.get('fit_s_covQual')!=3 for r in rows),bounds=sum(bool(r.get('parameters_at_bounds')) for r in rows)),indent=2))
