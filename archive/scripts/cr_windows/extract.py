from pathlib import Path
import ROOT,json,math,hashlib
ROOT.gSystem.Load('libHiggsAnalysisCombinedLimit');P=Path.cwd();R=P/'runs/20261029_cr_windows';configs=[dict(tag='baseline',range=[8.8,10.8],bins=150)]+json.loads((R/'configs.json').read_text());eras=['2022','2022EE','2023','2023BPix'];sr=json.loads((P/'runs/20261029_upsilon_control_cheb3/signal/audit.json').read_text());out={}
for cfg in configs:
 out[cfg['tag']]={}
 for n in [2,3]:
  rr=P/'runs'/('20261029_upsilon_control' if n==2 else '20261029_upsilon_control_cheb3') if cfg['tag']=='baseline' else R/cfg['tag']/f'degree{n}'
  f=ROOT.TFile.Open(str(rr/'control/fitdiagnostics/fitDiagnostics_control.root'));fit=f.Get('fit_b');pars=fit.floatParsFinal();cov=fit.covarianceMatrix();wf=ROOT.TFile.Open(str(rr/'control/workspace.root'));w=wf.Get('w');rows=[]
  for p in pars:
   if w.var(p.GetName()):w.var(p.GetName()).setVal(p.getVal())
  for i,era in enumerate(eras):
   row=dict(era=era);names=[f'CMS_EXO25020_bkgScale_{proc}_{era}' for proc in ['upsilon1s','upsilon2s']];grad=[];prediction=0
   for proc,name in zip(['upsilon1s','upsilon2s'],names):
    p=pars.find(name);idx=pars.index(name);err=math.sqrt(cov[idx][idx]);g=sr['eras'][i]['processes'][proc]['local_rate_r1'];grad.append(g);prediction+=g*p.getVal();cors={f'c{j}':fit.correlation(name,f'local_c{j}_{era}') for j in range(1,n+1)};cors['continuum_yield']=fit.correlation(name,'local_yield_'+era)
    row[proc]=dict(scale=p.getVal(),error=err,relative_error=err/p.getVal() if p.getVal()>1e-12 else None,SR_yield=g*p.getVal(),SR_error=g*err,continuum_correlations=cors,boundary=bool(p.getVal()<1e-4 or p.getVal()>9.9999))
   var=sum(grad[a]*grad[b]*cov[pars.index(names[a])][pars.index(names[b])] for a in range(2) for b in range(2));row.update(SR_total=prediction,SR_total_error=math.sqrt(max(0,var)),correlation_Y1_Y2=fit.correlation(*names));rows.append(row)
  record=dict(range=cfg['range'],bins=cfg['bins'],status=fit.status(),covQual=fit.covQual(),edm=fit.edm(),eras=rows,SR_sum=sum(x['SR_total'] for x in rows),SR_sum_error=math.sqrt(sum(x['SR_total_error']**2 for x in rows)))
  # Original shapes have no shared nuisance dependence in B-only; verify inter-era scale covariances negligible.
  names=[f'CMS_EXO25020_bkgScale_{proc}_{era}' for era in eras for proc in ['upsilon1s','upsilon2s']];record['cross_era_max_correlation']=max(abs(fit.correlation(a,b)) for a in names for b in names if a.split('_')[-1]!=b.split('_')[-1]);assert record['cross_era_max_correlation']<1e-5
  out[cfg['tag']][str(n)]=record
  d=R/cfg['tag']/f'degree{n}'/'gof';d.mkdir(parents=True,exist_ok=True);w.var('r').setVal(0);w.saveSnapshot('control_fit',w.allVars());w.writeToFile(str(d/'workspace.root'))
  audit=json.loads((rr/'control/audit.json').read_text());assert hashlib.sha256(Path(audit['source']).read_bytes()).hexdigest()==audit['source_sha256'];wf.Close();f.Close()
(R/'control_comparison.json').write_text(json.dumps(out,indent=2));print(json.dumps({k:{n:{a:b for a,b in v.items() if a!='eras'} for n,v in vals.items()} for k,vals in out.items()},indent=2))
