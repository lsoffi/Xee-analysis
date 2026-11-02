from pathlib import Path
import json,csv,math
import numpy as np
from scipy.special import erf
from scipy.optimize import minimize_scalar,brentq
from scipy.integrate import quad
P=Path.cwd();R=P/'runs/20261029_upsilon_peak_windows';p=json.loads((R/'parameters.json').read_text());ys=p['upsilon'];signals=p['signal'];eras=['2022','2022EE','2023','2023BPix'];fractions=[.68,.80,.90]
def primitive(x,r):
 scalar=np.ndim(x)==0;x=np.atleast_1d(np.asarray(x,dtype=float));sig=np.where(x<r['mean'],r['sigmaL'],r['sigmaR']);t=(x-r['mean'])/sig;v=np.sqrt(np.pi/2)*erf(t/np.sqrt(2))
 for side,sgn in [('L',-1),('R',1)]:
  a=abs(r['alpha'+side]);n=r['n'+side];z=sgn*t;mask=z>a;B=n/a-a;A=(n/a)**n*np.exp(-a*a/2)
  if abs(n-1)<1e-8:tail=A*np.log((B+z[mask])/(B+a))
  else:tail=A*((B+z[mask])**(1-n)-(B+a)**(1-n))/(1-n)
  v[mask]=sgn*(np.sqrt(np.pi/2)*erf(a/np.sqrt(2))+tail)
 v*=sig;return float(v[0]) if scalar else v

def cdf(x,r):return (primitive(x,r)-primitive(6,r))/(primitive(12,r)-primitive(6,r))
def integral(lo,hi,r):return float(cdf(hi,r)-cdf(lo,r))
def density(x,r):
 sig=r['sigmaL'] if x<r['mean'] else r['sigmaR'];t=(x-r['mean'])/sig
 if t<-abs(r['alphaL']):a=abs(r['alphaL']);n=r['nL'];return np.exp(-a*a/2)*(1+(-t-a)*a/n)**(-n)
 if t>abs(r['alphaR']):a=abs(r['alphaR']);n=r['nR'];return np.exp(-a*a/2)*(1+(t-a)*a/n)**(-n)
 return np.exp(-t*t/2)

def interval(r,f):
 grid=np.linspace(6,12,120001);c=cdf(grid,r);us=np.linspace(0,1-f,10001);width=np.interp(us+f,c,grid)-np.interp(us,c,grid);i=int(np.argmin(width));a=us[max(0,i-2)];b=us[min(len(us)-1,i+2)]
 def quant(u):return 6. if u<=0 else (12. if u>=1 else brentq(lambda x:cdf(x,r)-u,6,12,xtol=1e-12))
 def objective(u):return quant(u+f)-quant(u)
 opt=minimize_scalar(objective,bounds=(a,b),method='bounded',options={'xatol':1e-12});u=min([0,1-f,float(opt.x)],key=objective);lo=quant(u);hi=quant(u+f);assert abs(integral(lo,hi,r)-f)<1e-9;return lo,hi

def merge(intervals):
 out=[]
 for lo,hi in sorted(intervals):
  if out and lo<=out[-1][1]+1e-10:out[-1][1]=max(hi,out[-1][1])
  else:out.append([lo,hi])
 return out

def fraction_union(intervals,r):return sum(integral(lo,hi,r) for lo,hi in merge(intervals))
# Independent ROOT point-density checks and numerical quadrature.
validation=dict(max_root_density_difference=0.,max_quad_integral_difference=0.)
for r in ys+signals:
 norm=primitive(12,r)-primitive(6,r)
 for x,value in r['checks']:validation['max_root_density_difference']=max(validation['max_root_density_difference'],abs(density(x,r)/norm-value))
 points=[r['mean']-abs(r['alphaL'])*r['sigmaL'],r['mean'],r['mean']+abs(r['alphaR'])*r['sigmaR']]
 numerical=quad(density,6,12,args=(r,),points=[x for x in points if 6<x<12],epsabs=1e-11,epsrel=1e-10)[0];validation['max_quad_integral_difference']=max(validation['max_quad_integral_difference'],abs(numerical/norm-1))
assert validation['max_root_density_difference']<1e-8 and validation['max_quad_integral_difference']<1e-8
perera=[];common=[]
for f in fractions:
 pieces=[]
 for proc in ['upsilon1s','upsilon2s']:
  rows=[]
  for r in [r for r in ys if r['process']==proc]:
   lo,hi=interval(r,f);row=dict(target_fraction=f,process=proc,era=r['era'],low=lo,high=hi,width=hi-lo,mean=r['mean'],sigma_core=r['sigmaL']);rows.append(row);perera.append(row)
  lo=min(v['low'] for v in rows);hi=max(v['high'] for v in rows);blo=6+math.floor((lo-6)*75+1e-9)/75;bhi=6+math.ceil((hi-6)*75-1e-9)/75
  coverage={r['era']:integral(blo,bhi,r) for r in ys if r['process']==proc};assert min(coverage.values())>=f-1e-9
  pieces.append(dict(process=proc,low=blo,high=bhi,width=bhi-blo,coverage=coverage,exact_envelope=[lo,hi]))
 unions=merge([[x['low'],x['high']] for x in pieces]);width=sum(hi-lo for lo,hi in unions);common.append(dict(target_fraction=f,pieces=pieces,union=unions,union_width=width,gap=max(0,pieces[1]['low']-pieces[0]['high']),overlap=max(0,pieces[0]['high']-pieces[1]['low']),recovered_vs_current_CR=(10.506666666666668-9.093333333333334)-width,coverage_union={r['process']+'_'+r['era']:fraction_union(unions,r) for r in ys}))
cores=[]
for scale in [1.,1.5]:
 pieces=[]
 for proc in ['upsilon1s','upsilon2s']:
  shapes=[r for r in ys if r['process']==proc];lo=min(r['mean']-scale*r['sigmaL'] for r in shapes);hi=max(r['mean']+scale*r['sigmaR'] for r in shapes);lo=6+math.floor((lo-6)*75+1e-9)/75;hi=6+math.ceil((hi-6)*75-1e-9)/75;pieces.append(dict(process=proc,low=lo,high=hi,width=hi-lo,coverage={r['era']:integral(lo,hi,r) for r in shapes}))
 unions=merge([[x['low'],x['high']] for x in pieces]);width=sum(b-a for a,b in unions);cores.append(dict(label='core'+str(scale),sigma_core_multiplier=scale,pieces=pieces,union=unions,union_width=width,gap=max(0,pieces[1]['low']-pieces[0]['high']),overlap=max(0,pieces[0]['high']-pieces[1]['low']),recovered_vs_current_CR=(10.506666666666668-9.093333333333334)-width,coverage_union={r['process']+'_'+r['era']:fraction_union(unions,r) for r in ys}))
(R/'core_windows.json').write_text(json.dumps(cores,indent=2))
# Include current contiguous CR as reference; fractions are conditional on original6–12 support.
candidates=[dict(label='current',union=[[9.093333333333334,10.506666666666668]])]+[dict(label=str(round(100*r['target_fraction'])),union=r['union']) for r in common]+[dict(label=r['label'],union=r['union']) for r in cores]
rows=[]
for r in signals:
 for cand in candidates:
  frac=fraction_union(cand['union'],r);rows.append(dict(mass=r['mass'],era=r['era'],candidate=cand['label'],signal_fraction_in_CR=frac,signal_fraction_outside_CR=1-frac))
# Consequences for existing common95 search windows, including retained signal in SR minus CR.
old=json.loads((P/'runs/20261029_signal_widths/common95.json').read_text())['windows'];srrows=[]
for window in old:
 m=window['mass_GeV'];lo=6+math.floor((window['low95_common_GeV']-6)*75+1e-9)/75;hi=6+math.ceil((window['high95_common_GeV']-6)*75-1e-9)/75
 for cand in candidates:
  overlap=merge([[max(a,lo),min(b,hi)] for a,b in cand['union'] if max(a,lo)<min(b,hi)]);inside=[];outside=[]
  for r in [r for r in signals if r['mass']==m]:
   total=integral(lo,hi,r);removed=fraction_union(overlap,r);inside.append(total);outside.append(total-removed)
  srrows.append(dict(mass=m,candidate=cand['label'],original_SR=[lo,hi],overlap_width=sum(b-a for a,b in overlap),overlap_intervals=overlap,disjoint=not overlap,retained_signal_SR_min=min(outside),retained_signal_SR_max=max(outside),original_signal_SR_min=min(inside)))
(R/'per_era_windows.json').write_text(json.dumps(perera,indent=2));(R/'common_windows.json').write_text(json.dumps(common,indent=2));(R/'signal_fractions.json').write_text(json.dumps(rows,indent=2));(R/'SR_overlap.json').write_text(json.dumps(srrows,indent=2));(R/'validation.json').write_text(json.dumps(validation,indent=2))
for filename,data in [('per_era_windows.csv',perera),('signal_fractions.csv',rows),('SR_overlap.csv',srrows)]:
 with (R/filename).open('w') as h:w=csv.DictWriter(h,fieldnames=list(data[0]));w.writeheader();w.writerows(data)
print(json.dumps(common,indent=2));print('VALIDATION',validation)
