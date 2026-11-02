from pathlib import Path
import argparse,json,math,re,fcntl
import ROOT
from observed_helpers import dump,read_stage,read_tree
ROOT.gROOT.SetBatch(True);ROOT.gSystem.Load('libHiggsAnalysisCombinedLimit')
P=Path(__file__).resolve().parents[2];R=P/'runs/20261031_shared_sr7_floatmu_fitlimits'
parser=argparse.ArgumentParser();parser.add_argument('mass');a=parser.parse_args();m=f'{float(a.mass):.1f}'
D=R/'observed'/('m'+m);J=D/'joint';workspace=J/'workspace_seeded.root'
lockfile=(D/'extraction.lock').open('a');fcntl.flock(lockfile,fcntl.LOCK_EX)
result=dict(mass=float(m),degree=7,fit={},flags=[],expected_validated=False,status='partial')
curves={};candidates={}
for stage in ['b','s']:
    trial=[]
    for directory in sorted(J.glob('fitdiagnostics_anchor*')):
        f=directory/'fitDiagnostics_shared.root'
        if f.exists() and (directory/'combine.receipt.json').exists():
            x,bins=read_stage(workspace,f,stage)
            if x.get('readable'):trial.append(x)
    if not trial:continue
    usable=[x for x in trial if x['point_valid']]
    finite=[x for x in trial if math.isfinite(x['exact_deviance'])]
    selected=min(usable or finite or trial,key=lambda x:x['exact_deviance'])
    value,bins=read_stage(workspace,selected['source'],stage,J/f'workspace_postfit_{stage}.root')
    result['fit'][stage]=value;curves[stage]=bins;candidates[stage]=trial
    if not value['numerical_ok']:result['flags'].append('fit_'+stage+'_numerical_flag')
    if not value['edm_valid']:result['flags'].append('fit_'+stage+'_high_or_nonfinite_EDM')
    if value['mean_boundary']:result['flags'].append('fit_'+stage+'_upsilon_mean_boundary')
    if stage=='s':result.update(r=value['r'],r_error=value['error'])
if len(result['fit'])==2:
    delta=result['fit']['b']['exact_deviance']-result['fit']['s']['exact_deviance']
    result['Z_FD']=math.sqrt(max(0,delta)) if result['r']>0 else 0
    result['delta_deviance_B_minus_SB']=delta
    if delta<-.01:result['flags'].append('S_minimum_worse_than_B')
files=sorted((J/'limits').glob('higgsCombine*.root')) if (J/'limits/combine.receipt.json').exists() else []
if files:
    values=read_tree(files[-1]);L={}
    for key,q in [('observed',-1),('expected_2p5',.025),('expected_16',.16),('expected_50',.5),('expected_84',.84),('expected_97p5',.975)]:
        near=[v for v in values if abs(v['quantile']-q)<.004]
        if near:L[key]=near[0]['limit']
    if L:result['limits']=L
    logfile=J/'limits/combine.log';text=logfile.read_text(errors='replace') if logfile.exists() else ''
    pairs=re.findall(r'At r =\s*([\d.eE+-]+):[^\n]*CLs\s*=\s*([\d.eE+-]+)',text)
    if 'observed' in L and pairs:
        rr,cl=min(pairs,key=lambda v:abs(float(v[0])-L['observed']))
        result['observed_limit_check']=dict(r=float(rr),CLs=float(cl),verified=abs(float(cl)-.05)<.002)
    else:result['observed_limit_check']=dict(verified=False)
    if not result['observed_limit_check']['verified']:result['flags'].append('observed_limit_crossing_unverified')
    result['expected_warnings']=[s.strip() for s in text.splitlines() if 'Best fit of asimov dataset' in s or 'WARNING' in s]
    if len(L)!=6:result['flags'].append('incomplete_limit_quantiles')
    else:
        exp=[L[k] for k in ['expected_2p5','expected_16','expected_50','expected_84','expected_97p5']]
        if not all(0<x<=y for x,y in zip(exp,exp[1:])):result['flags'].append('unordered_expected_quantiles')
        sig=L['expected_50']/1.959963984540054
        result['expected_bracket_hits']=[key for key,n in [('expected_84',2),('expected_97p5',4)] if abs(L[key]/(L['expected_50']+n*sig)-1)<.003]
        if result['expected_bracket_hits']:result['flags'].append('expected_bands_bracket_hit')
result['method_receipts']={str(p.parent.relative_to(J)):json.loads(p.read_text()) for p in J.glob('*/combine.receipt.json')}
for name,receipt in result['method_receipts'].items():
    if receipt.get('exit_code')!=0:result['flags'].append(name+'_nonzero_exit')
result['numerical_ok']=len(result['fit'])==2 and all(v['numerical_ok'] for v in result['fit'].values()) and result.get('observed_limit_check',{}).get('verified',False)
result['status']='complete' if len(result.get('limits',{}))==6 and len(result['fit'])==2 and result['method_receipts'].get('limits',{}).get('exit_code')==0 else 'partial'
f=ROOT.TFile.Open(str(workspace));w=f.Get('w')
result['mu_prefit']={v.GetName():float(v.getVal()) for v in w.allVars() if v.GetName().startswith('Ymean_')}
f.Close()
dump(D/'results.json',result);dump(D/'bins.json',curves);dump(D/'fit_trials.json',candidates)
print(m,result['status'],result.get('r'),result['flags'],flush=True)
