from pathlib import Path
import ROOT,json,math
ROOT.gROOT.SetBatch(True);ROOT.gSystem.Load('libHiggsAnalysisCombinedLimit')
R=Path('runs/20261030_upsilon_mu_sigma_cheb3');out=[]
for a in json.loads((R/'configs.json').read_text()):
 d=R/a['tag']
 if not (d/'multidimfit_shape.root').exists():print('MISSING',a['tag']);continue
 a.update(json.loads((d/'audit.json').read_text()));f=ROOT.TFile.Open(str(d/'multidimfit_shape.root'));r=f.Get('fit_mdf');fw=ROOT.TFile.Open(str(d/'workspace.root'));w=fw.Get('w');pars=r.floatParsFinal();a['status']=r.status();a['covQual']=r.covQual();a['parameters']={};a['boundaries']=[]
 for p in pars:
  w.var(p.GetName()).setVal(p.getVal());a['parameters'][p.GetName()]={'value':p.getVal(),'lo':p.getErrorLo(),'hi':p.getErrorHi()}
  if min(p.getVal()-p.getMin(),p.getMax()-p.getVal())<1e-4*(p.getMax()-p.getMin()):a['boundaries'].append(p.GetName())
 log=(d/'combine.log').read_text();a['profile_failed']='Minos error calculation failed' in log or (d/'profile_missing').exists();a['reliable_errors']=a['status']==0 and a['covQual']==3 and not a['profile_failed'];x=w.var('mass');obs=ROOT.RooArgSet(x);bw=(a['range'][1]-a['range'][0])/a['bins'];curves={k:[] for k in ['Y1','Y2','continuum','total']}
 for point in a['data']:
  x.setVal(point['mass']);y1=w.pdf('upsilon1s_bin').getVal(obs)*bw*w.var('N_upsilon1s').getVal();y2=w.pdf('upsilon2s_bin').getVal(obs)*bw*w.function('N_upsilon2s').getVal()
  bg=w.pdf('continuum_shape').getVal(obs)*bw*w.var('N_continuum').getVal() if a['background_degree']==3 else sum(w.pdf('continuum_'+side).getVal(obs)*bw*w.var('N_continuum_'+side).getVal() for side in ['left','right'])
  for k,v in zip(curves,[y1,y2,bg,y1+y2+bg]):curves[k].append(v)
 a['curves']=curves;a['expected_sum']=sum(curves['total']);out.append(a);fw.Close();f.Close()
(R/'results.json').write_text(json.dumps(out,indent=2))
print('fits',len(out),'status0cov3',sum(a['status']==0 and a['covQual']==3 for a in out),'profile_problem',[(a['tag'],a['status'],a['covQual']) for a in out if not a['reliable_errors']])
