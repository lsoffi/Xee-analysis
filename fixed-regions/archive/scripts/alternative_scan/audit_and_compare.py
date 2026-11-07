import ROOT,pathlib,json,csv,hashlib,re,math
import numpy as np
ROOT.gROOT.SetBatch(True);ROOT.gSystem.Load('libHiggsAnalysisCombinedLimit');P=pathlib.Path('/Users/imb/Documents/Codex/Analisi-statistica-Xee');R=P/'runs/20261029_alternative_background_scan';audits=[];sources={};shapechecks=[]
for a in json.loads((R/'audit_scan.json').read_text()):
 sources[a['source']]=a['source_sha256'];f=ROOT.TFile.Open(a['path']);w=f.Get('w');g=ROOT.TFile.Open(str(R/f"fits/{a['family']}/m{a['mass_label']}/fitdiagnostics/fitDiagnostics_alternative.root"));obs=ROOT.RooArgSet(w.var('mass'));x=np.linspace(6,12,1201)
 assert w.data('data_obs').sumEntries()==a['events'] and w.var('mass').getBins()==450
 for stage in ['fit_b','fit_s']:
  fit=g.Get(stage);w.allVars().assignValueOnly(fit.floatParsFinal())
  for era in ['2022','2022EE','2023','2023BPix']:
   pdf=w.pdf(f"ALT_{a['family']}_{era}_pdf");values=[]
   for v in x:w.var('mass').setVal(float(v));values.append(pdf.getVal(obs))
   values=np.array(values);area=float(np.trapezoid(values,x));shapechecks.append(dict(family=a['family'],mass=a['mass'],stage=stage,era=era,min_pdf=float(values.min()),integral_grid=area,positive_finite=bool(np.isfinite(values).all() and values.min()>=0),normalization_grid_ok=bool(abs(area-1)<2e-4)))
 g.Close();f.Close()
for path,sha in sources.items():assert hashlib.sha256(pathlib.Path(path).read_bytes()).hexdigest()==sha
(R/'model_integrity_final.json').write_text(json.dumps(dict(source_files_unchanged=len(sources),workspaces=204,shape_checks=len(shapechecks),shape_anomalies=[r for r in shapechecks if not r['positive_finite'] or not r['normalization_grid_ok']],checks=shapechecks),indent=2))
comparisons=[]
for a in json.loads((R/'refinement_points.json').read_text()):
 if not a['reason']['covariance_or_status']:continue
 d=R/f"numerical_checks/{a['family']}/m{a['mass_label']}";row=dict(family=a['family'],mass=a['mass'])
 for label,path in [('primary',R/f"fits/{a['family']}/m{a['mass_label']}"),('accurate',d)]:
  name='alternative' if label=='primary' else 'accurate';ff=ROOT.TFile.Open(str(path/f'fitdiagnostics/fitDiagnostics_{name}.root'));fb=ff.Get('fit_b');fs=ff.Get('fit_s');row[label+'_b_status']=fb.status();row[label+'_s_status']=fs.status();row[label+'_b_covQual']=fb.covQual();row[label+'_s_covQual']=fs.covQual();ff.Close()
  log=(path/'significance/combine.log').read_text();m=re.search(r'Significance:\s*([\d.eE+-]+)',log);row[label+'_Z']=float(m.group(1)) if m else None
 row['Z_difference']=row['accurate_Z']-row['primary_Z'];comparisons.append(row)
(R/'numerical_comparison.json').write_text(json.dumps(comparisons,indent=2))
with (R/'numerical_comparison.csv').open('w') as f:wr=csv.DictWriter(f,fieldnames=list(comparisons[0]));wr.writeheader();wr.writerows(comparisons)
print(json.dumps(dict(source_files_unchanged=len(sources),shape_checks=len(shapechecks),shape_anomalies=sum(not a['positive_finite'] or not a['normalization_grid_ok'] for a in shapechecks),numerical_checks=len(comparisons),max_Z_difference=max(abs(r['Z_difference']) for r in comparisons)),indent=2))
