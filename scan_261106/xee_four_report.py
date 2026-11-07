from pathlib import Path
import json,csv,math
import numpy as np
import ROOT
from scipy.integrate import quad
ROOT.gSystem.Load('libHiggsAnalysisCombinedLimit')
web=Path('/eos/user/s/soffi/www/Xee/261106/scan_SR_CR')
out=web/'confronto_4_modelli_low';out.mkdir(exist_ok=True)
rows=[]
for m in ('0.4','0.5','0.6'):
 md=json.loads((Path('low_logpoly3/results')/('m'+m)/'model.json').read_text())
 for label in ('Cheb2','Cheb3','Cheb6','Log3'):
  name='low_logpoly3'+('_totalC0' if m=='0.6' else '') if label=='Log3' else ('low_'+label.lower()+'_totalC0' if label!='Cheb6' or m=='0.6' else '')
  d=Path(name)/'results'/('m'+m)
  f=ROOT.TFile.Open(str(next(d.glob('fitDiagnostics*.root'))))
  if m=='0.6':
   wf=ROOT.TFile.Open(str(d/'workspace.root'));w=wf.Get('w');meta=json.loads((d/'model.json').read_text())
  for stage in ('b','s'):
   fit=f.Get('fit_'+stage);pars=fit.floatParsFinal()
   chi=0;bins=0
   for ch in md['channels']:
    if ch['region']!='SR0':continue
    h=f.Get('shapes_fit_'+stage+'/'+ch['label']+'/total')
    for i,y in enumerate(ch['counts'],1):
     v=h.GetBinContent(i);assert v>0
     chi+=(y-v)**2/v;bins+=1
   jump=0.
   if m=='0.6':
    for par in pars:w.var(par.GetName()).setVal(par.getVal())
    def get(k):return w.arg(k).getVal()
    for era in meta['eras']:
     e=era['era'];g=next(q for q in era['CRs'] if q['group']=='omega_phi')
     edge=.65
     amp=get('logpoly_anchor_'+e) if label=='Log3' else sum(get(n) for n in era['continuum_coefficients'])
     other=0.
     for res in ('jpsi','psi2s'):
      mu=get('mu_'+res+'_'+e);sig,al,nl,ar,nr=[get(res+'_'+k+'_'+e) for k in ('sigma','alphaL','nL','alphaR','nR')]
      def pdf(x):
       t=(x-mu)/sig
       if t < -al:return math.exp(-al*al/2-nl*math.log1p(al*(-al-t)/nl))
       if t > ar:return math.exp(-ar*ar/2-nr*math.log1p(ar*(t-ar)/nr))
       return math.exp(-t*t/2)
      area=quad(pdf,2.8,4,points=[x for x in [mu-sig*al,mu+sig*ar] if 2.8<x<4],epsabs=1e-13)[0]
      other+=get('N_'+res+'_'+e)*pdf(edge)/area
     slope=g['slope'];delta=g['range'][1]-edge
     unit=slope/math.expm1(slope*delta) if slope else 1/delta
     jump=max(jump,abs(amp+other-get(g['continuum_norm'])*unit)*.005)
    assert jump<1e-5,(label,jump)
   r=pars.find('r')
   rows.append([m,label,stage,round(chi,2),bins,fit.status(),fit.covQual(),round(r.getVal(),4) if r else 0,jump if m=='0.6' else ''])
  if m=='0.6':wf.Close()
  f.Close()
with (out/'confronto.csv').open('w') as fp:
 wr=csv.writer(fp);wr.writerow(['massa','modello','fit','chi2_Pearson_SR','bin_SR','status','covQual','r','max_salto_eventi_per_5MeV']);wr.writerows(rows)
for row in rows:print(*row)
