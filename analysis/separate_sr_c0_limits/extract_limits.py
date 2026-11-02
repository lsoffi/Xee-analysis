from pathlib import Path
import argparse,json,math,re,os
import ROOT
ROOT.gROOT.SetBatch(True)
P=Path(__file__).resolve().parents[2];R=P/'runs/20261031_separate_sr_c0_limits'
a=argparse.ArgumentParser();a.add_argument('region');a.add_argument('mass');args=a.parse_args();m=f'{float(args.mass):.1f}'
D=R/args.region/('m'+m);J=D/'limits';receipt=json.loads((J/'combine.receipt.json').read_text());rows=[]
files=sorted(J.glob('higgsCombine*.root'))
if files:
 f=ROOT.TFile.Open(str(files[-1]))
 if f and not f.IsZombie():
  tree=f.Get('limit')
  if tree:rows=[dict(limit=float(v.limit),q=float(v.quantileExpected)) for v in tree]
  f.Close()
L={}
for key,q in [('observed',-1),('expected_2p5',.025),('expected_16',.16),('expected_50',.5),('expected_84',.84),('expected_97p5',.975)]:
 near=[v for v in rows if abs(v['q']-q)<.004]
 if near and math.isfinite(near[0]['limit']):L[key]=near[0]['limit']
flags=[]
if receipt['exit_code']!=0:flags.append('native_nonzero_exit')
if len(L)!=6:flags.append('incomplete_limit_quantiles')
text=(J/'combine.log').read_text(errors='replace');pairs=re.findall(r'At r =\s*([\d.eE+-]+):[^\n]*CLs\s*=\s*([\d.eE+-]+)',text)
check={'verified':False}
if 'observed' in L and pairs:
 rr,cl=min(pairs,key=lambda v:abs(float(v[0])-L['observed']));check=dict(r=float(rr),CLs=float(cl),verified=abs(float(cl)-.05)<.002)
if not check['verified']:flags.append('observed_limit_crossing_unverified')
bracket=[]
if len(L)==6:
 exp=[L[k] for k in ['expected_2p5','expected_16','expected_50','expected_84','expected_97p5']]
 if not all(0<x<=y for x,y in zip(exp,exp[1:])):flags.append('unordered_expected_quantiles')
 sig=L['expected_50']/1.959963984540054
 if sig>0:bracket=[k for k,n in [('expected_84',2),('expected_97p5',4)] if abs(L[k]/(L['expected_50']+n*sig)-1)<.003]
 if bracket:flags.append('expected_bands_bracket_hit')
warnings=[s.strip() for s in text.splitlines() if 'WARNING' in s or 'Best fit of asimov dataset' in s]
if any('Best fit of asimov dataset' in s for s in warnings):flags.append('nonzero_asimov_bestfit_warning')
a=dict(region=args.region,mass=float(m),degree=6 if args.region=='low' else 4,limits=L,flags=flags,receipt=receipt,
 observed_limit_check=check,expected_bracket_hits=bracket,native_warnings=warnings,expected_validated=False,
 status='complete' if len(L)==6 and receipt['exit_code']==0 else 'partial',
 note='Only native AsymptoticLimits requested; no independent FitDiagnostics minimum/covariance audit. Shared CR makes the two analyses correlated.')
p=D/'results.json';tmp=p.with_suffix('.tmp.'+str(os.getpid()));tmp.write_text(json.dumps(a,indent=2));tmp.replace(p)
print(args.region,m,a['status'],L,flags)
