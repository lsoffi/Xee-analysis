"""Read the native significance; FD is a check, never a replacement result."""
from pathlib import Path
import argparse,json,math,os,re
import ROOT

ROOT.gROOT.SetBatch(True)
P=Path(__file__).resolve().parents[2];R=P/'runs/20261031_separate_sr_c0_limits'

def log_probability(Z,p0):
    if p0>0:return math.log(p0)
    # erfc underflows for very large positive Z. The alternating Mills-ratio
    # expansion is accurate to floating precision here (Z>~38); stop before
    # the asymptotic terms grow or once their contribution is negligible.
    assert Z>30
    series=1.;term=1.
    for k in range(1,100):
        nxt=-term*(2*k-1)/(Z*Z)
        if abs(nxt)>=abs(term):break
        series+=nxt;term=nxt
        if abs(term)<1e-17*abs(series):break
    return -.5*Z*Z-math.log(Z)-.5*math.log(2*math.pi)+math.log(series)

def extract(side,mass):
    mass=float(mass);D=R/side/f'm{mass:.1f}';J=D/'significance'
    receipt=json.loads((J/'combine.receipt.json').read_text())
    fd=json.loads((D/'fit_results.json').read_text());flags=[];values=[];output=None
    files=sorted(J.glob('higgsCombine*.Significance.*.root'))
    if len(files)==1:
        output=str(files[0]);f=ROOT.TFile.Open(output)
        if f and not f.IsZombie():
            tree=f.Get('limit')
            if tree:values=[float(row.limit) for row in tree]
        if f:f.Close()
    if receipt.get('exit_code')!=0:flags.append('native_nonzero_exit')
    Z=values[0] if len(values)==1 and math.isfinite(values[0]) else None
    if Z is None:flags.append('missing_or_invalid_native_significance')
    elif Z<0:flags.append('negative_native_significance')
    p0=.5*math.erfc(Z/math.sqrt(2.)) if Z is not None else None
    log_p0=log_probability(Z,p0) if Z is not None else None
    underflow=bool(Z is not None and p0==0.)
    if underflow:flags.append('p0_floating_point_underflow')
    ref_valid=bool(fd['point_valid']);rhat=fd['fit']['s']['r'];delta=fd['fit']['b']['exact_deviance']-fd['fit']['s']['exact_deviance']
    ref_Z=math.sqrt(max(0.,delta)) if rhat>0 else 0.
    if not ref_valid:flags.append('FD_reference_point_invalid')
    if any(fd['fit'][stage]['covQual']!=3 for stage in ['b','s']):flags.append('FD_covariance_warning')
    if any(fd['fit'][stage]['mean_boundary'] for stage in ['b','s']):flags.append('FD_Upsilon_mean_boundary')
    agreement=dict(checked=ref_valid and Z is not None,tolerance=.1)
    if agreement['checked']:
        agreement.update(absolute_Z_difference=abs(Z-ref_Z),agrees=abs(Z-ref_Z)<=.1)
        if not agreement['agrees']:flags.append('native_FD_significance_mismatch')
    else:agreement.update(absolute_Z_difference=None,agrees=None)
    text=(J/'combine.log').read_text(errors='replace');printed=re.findall(r'Significance:\s*([\d.eE+-]+)',text)
    if printed and Z is not None and abs(float(printed[-1])-Z)>1e-4:flags.append('native_log_root_mismatch')
    complete=receipt.get('exit_code')==0 and Z is not None
    result=dict(side=side,mass=mass,method='Significance',Z=Z,p0=p0,log10_p0=log_p0/math.log(10) if log_p0 is not None else None,
                p0_underflow=underflow,p0_underflow_note='Numerical underflow in erfc, not a probability exactly zero; log10_p0 remains finite.' if underflow else None,
                status='complete' if complete else 'failed',receipt=receipt,flags=flags,native_output=output,
                reference_FD=dict(point_valid=ref_valid,rhat=rhat,Z=ref_Z,delta_deviance=delta,
                    status={s:fd['fit'][s]['status'] for s in ['b','s']},covQual={s:fd['fit'][s]['covQual'] for s in ['b','s']}),
                FD_agreement=agreement,numerical_ok=bool(complete and ref_valid and agreement['agrees']),
                covariance_ok=all(fd['fit'][s]['covQual']==3 for s in ['b','s']),
                interpretation='Native one-sided local significance and p0; FD only checks the minimum. No look-elsewhere correction. Correlated low/high analyses remain separate.')
    path=D/'significance_results.json';temp=path.with_suffix('.tmp.'+str(os.getpid()));temp.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n');temp.replace(path)
    print(side,f'{mass:.1f}',Z,p0,flags,flush=True);return result

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('side',choices=['low','high']);p.add_argument('mass');a=p.parse_args();extract(a.side,a.mass)
