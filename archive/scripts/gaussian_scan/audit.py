import ROOT,pathlib,json,csv,hashlib,math
import numpy as np
ROOT.gROOT.SetBatch(True);ROOT.gSystem.Load('libHiggsAnalysisCombinedLimit');P=pathlib.Path('/Users/imb/Documents/Codex/Analisi-statistica-Xee');R=P/'runs/20261029_gaussian_signal_scan';ERAS=['2022','2022EE','2023','2023BPix'];sources={};checks=[];post=[];definitions=json.loads((R/'signal_shape_definitions.json').read_text())
for a in json.loads((R/'audit_scan.json').read_text()):
 sources[a['source']]=a['source_sha256'];f=ROOT.TFile.Open(a['path']);w=f.Get('w');s=ROOT.TFile.Open(a['source']);orig=s.Get('w');obs=ROOT.RooArgSet(w.var('mass'));obs0=ROOT.RooArgSet(orig.var('mass'));w.var('r').setVal(1);orig.var('r').setVal(1);assert w.data('data_obs').sumEntries()==orig.data('data_obs').sumEntries()==318924
 for i,era in enumerate(ERAS,1):
  dcb=w.pdf(f'shapeSig_Zd_ch{i}');g=w.pdf('gaussian_signal_'+era);assert g.getMean().GetName()==dcb.x0().GetName();assert g.getSigma().GetName()==dcb.sigmaL().GetName()==dcb.sigmaR().GetName();assert abs(w.function(f'n_exp_binch{i}_proc_Zd').getVal()-orig.function(f'n_exp_binch{i}_proc_Zd').getVal())<1e-9
  mu=g.getMean().getVal();sigma=g.getSigma().getVal();assert sigma>0;norm=.5*(math.erf((12-mu)/(sigma*math.sqrt(2)))-math.erf((6-mu)/(sigma*math.sqrt(2))));errs=[]
  for x in np.clip(mu+sigma*np.array([-2,-1,0,1,2]),6,12):
   w.var('mass').setVal(float(x));manual=math.exp(-.5*((x-mu)/sigma)**2)/(sigma*math.sqrt(2*math.pi)*norm);errs.append(abs(g.getVal(obs)-manual)/manual)
  for proc in ['upsilon1s','upsilon2s']:
   p=w.pdf(f'shapeBkg_{proc}_ch{i}');p0=orig.pdf(f'shapeBkg_{proc}_ch{i}')
   for x in [6.2,7.5,9.46,10.02,11.8]:w.var('mass').setVal(x);orig.var('mass').setVal(x);assert abs(p.getVal(obs)-p0.getVal(obs0))<1e-12
  checks.append(dict(model=a['model'],mass=a['mass'],era=era,mean=mu,sigma=sigma,gaussian_analytic_error=max(errs),same_mean_sigma_nodes=True,signal_rate_unchanged=True,upsilon_shapes_unchanged=True));assert max(errs)<1e-10
  if a['model']=='gaussian' and a['mass'] in [7.5,9.5]:
   x=np.linspace(6,12,1801);dc=[];gg=[]
   for value in x:w.var('mass').setVal(float(value));dc.append(dcb.getVal(obs));gg.append(g.getVal(obs))
   np.savetxt(R/f'nominal_shapes_m{a["mass_label"]}_{era}.csv',np.c_[x,dc,gg],delimiter=',',header='mass_GeV,dcb_density,gaussian_density',comments='')
 result=ROOT.TFile.Open(str(R/f'fits/{a["model"]}/m{a["mass_label"]}/fitdiagnostics/fitDiagnostics_alternative.root'))
 for stage in ['fit_b','fit_s']:
  fit=result.Get(stage);w.allVars().assignValueOnly(fit.floatParsFinal())
  for i,era in enumerate(ERAS,1):
   pdf=w.pdf(f'shapeBkg_dy_ch{i}').getPdf(1);xx=np.linspace(6,12,601);yy=[]
   for x in xx:w.var('mass').setVal(float(x));yy.append(pdf.getVal(obs))
   yy=np.array(yy);post.append(dict(model=a['model'],mass=a['mass'],era=era,stage=stage,min_pdf=float(np.nanmin(yy)),finite=bool(np.isfinite(yy).all()),positive=bool(np.isfinite(yy).all() and yy.min()>=0),integral_grid=float(np.trapezoid(yy,xx))))
 result.Close();s.Close();f.Close()
for path,sha in sources.items():assert hashlib.sha256(pathlib.Path(path).read_bytes()).hexdigest()==sha
out=dict(source_files_unchanged=len(sources),workspaces=102,nominal_checks=len(checks),checks=checks,postfit_continuum_checks=post,postfit_continuum_anomalies=[r for r in post if not r['positive']]);(R/'integrity_audit.json').write_text(json.dumps(out,indent=2))
with (R/'signal_widths.csv').open('w') as f:
 rows=[r for r in definitions if r['model']=='gaussian'];wr=csv.DictWriter(f,fieldnames=list(rows[0]));wr.writeheader();wr.writerows(rows)
print(json.dumps(dict(source_files_unchanged=len(sources),mean_sigma_rate_upsilon_checks=len(checks),postfit_continuum_anomalies=len(out['postfit_continuum_anomalies'])),indent=2))
