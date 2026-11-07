from pathlib import Path
import ROOT,json,csv,math
ROOT.gSystem.Load('libHiggsAnalysisCombinedLimit');P=Path('/Users/imb/Documents/Codex/Analisi-statistica-Xee');R=P/'runs/20261029_upsilon_control_bern3';eras=['2022','2022EE','2023','2023BPix'];names=['CMS_EXO25020_bkgScale_'+proc+'_'+era for era in eras for proc in ['upsilon1s','upsilon2s']];results={}
for mode,filename,stages in [('control','fitDiagnostics_control.root',['b']),('joint','fitDiagnostics_local.root',['b','s'])]:
 f=ROOT.TFile.Open(str(R/mode/'fitdiagnostics'/filename));wf=ROOT.TFile.Open(str(R/mode/'workspace.root'));w=wf.Get('w');result={}
 for stage in stages:
  fit=f.Get('fit_'+stage);pars=fit.floatParsFinal();values=[]
  for v in pars:
   target=w.var(v.GetName())
   if target:target.setVal(v.getVal())
  for era in eras:
   rates={}
   for proc in ['upsilon1s','upsilon2s']:
    name='CMS_EXO25020_bkgScale_'+proc+'_'+era;v=pars.find(name);rates[proc]=dict(scale=v.getVal(),error=math.sqrt(fit.covarianceMatrix()[pars.index(name)][pars.index(name)]),roo_error=v.getError())
    if mode=='control':
     full=w.function(f'n_exp_binch{eras.index(era)+1}_proc_{proc}');local=w.function('local_rate_'+proc+'_'+era);rates[proc].update(full_yield=full.getVal(),full_yield_error=full.getVal()/v.getVal()*rates[proc]['error'],control_yield=local.getVal(),control_yield_error=local.getVal()/v.getVal()*rates[proc]['error'])
   values.append(dict(era=era,**rates,correlation_1s_2s=fit.correlation(names[2*eras.index(era)],names[2*eras.index(era)+1])))
  corr=[[fit.correlation(a,b) for b in names] for a in names]
  result[stage]=dict(status=fit.status(),covQual=fit.covQual(),edm=fit.edm(),normalizations=values,correlation_names=names,correlations=corr,covariance=[[float(fit.covarianceMatrix()[pars.index(a)][pars.index(b)]) for b in names] for a in names],parameters=[dict(name=v.GetName(),value=v.getVal(),error=v.getError()) for v in pars])
 results[mode]=result;f.Close();wf.Close()
(R/'parameters.json').write_text(json.dumps(results,indent=2));print(json.dumps({k:{s:{a:b for a,b in val.items() if a not in ['parameters','correlations','correlation_names']} for s,val in v.items()} for k,v in results.items()},indent=2))
