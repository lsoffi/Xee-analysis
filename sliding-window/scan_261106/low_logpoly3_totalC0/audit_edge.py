from pathlib import Path
import json,math
import ROOT
from scipy.integrate import quad
ROOT.gSystem.Load('libHiggsAnalysisCombinedLimit')
base=Path.cwd()
for variant in ('low_logpoly3','low_logpoly3_totalC0'):
 d=base/variant/'results/m0.6'; wf=ROOT.TFile.Open(str(d/'workspace.root'));w=wf.Get('w');ff=ROOT.TFile.Open(str(next(d.glob('fitDiagnostics*.root'))));md=json.loads((d/'model.json').read_text())
 for stage in ('b','s'):
  fit=ff.Get('fit_'+stage)
  for par in fit.floatParsFinal():w.var(par.GetName()).setVal(par.getVal())
  chi=0;diffs=[]
  for em in md['eras']:
   era=em['era'];g=next(x for x in em['CRs'] if x['group']=='omega_phi')
   def get(k):return (w.arg(k)).getVal()
   edge=.65
   unit=g['slope']/math.expm1(g['slope']*(g['range'][1]-edge)) if g['slope'] else 1/(g['range'][1]-edge)
   extra=0
   for name in ('jpsi','psi2s'):
    mu=get('mu_'+name+'_'+era);sig,al,nl,ar,nr=[get(name+'_'+k+'_'+era) for k in ('sigma','alphaL','nL','alphaR','nR')]
    def density(x):
     t=(x-mu)/sig
     if t < -al:return math.exp(-al*al/2-nl*math.log1p(al*(-al-t)/nl))
     if t > ar:return math.exp(-ar*ar/2-nr*math.log1p(ar*(t-ar)/nr))
     return math.exp(-t*t/2)
    area=quad(density,2.8,4.,points=[x for x in [mu-sig*al,mu+sig*ar] if 2.8<x<4],epsabs=1e-13)[0]
    extra+=get('N_'+name+'_'+era)*density(edge)/area
   diffs.append((get('logpoly_anchor_'+era)+extra-get(g['continuum_norm'])*unit)*.005)
   ch=next(c for c in md['channels'] if c['label']=='SR0_'+era);h=ff.Get('shapes_fit_'+stage+'/'+ch['label']+'/total')
   for i,y in enumerate(ch['counts'],1):
    v=h.GetBinContent(i);chi+=(y-v)**2/v
  print(variant,stage,fit.status(),fit.covQual(),round(chi,3),diffs,flush=True)
  if variant.endswith('totalC0'):assert max(map(abs,diffs))<1e-5
 ff.Close();wf.Close()
