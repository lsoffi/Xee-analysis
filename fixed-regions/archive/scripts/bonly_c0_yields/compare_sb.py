"""Read original free-r fit_s snapshots and compare physical CR yields.

No new minimization. Covariance propagation includes r and all floating
parameters, unlike covariance from a conditional fixed-r likelihood point.
"""
from pathlib import Path
import csv
import hashlib
import json
import math
import sys
import numpy as np
import ROOT

P = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(P / 'scripts/separate_sr_c0_limits'))
from fit_helpers import ERAS, add_pdf, apply_fit, dump, inspect_model

ROOT.gROOT.SetBatch(True)
OUT = P / 'runs/20261102_c0_background_yields'
PROCS = ['upsilon1s', 'upsilon2s', 'Zd', 'dy']


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def extract(side, mass):
    d = P / 'runs/20261031_separate_sr_c0_limits' / side / f'm{mass:.1f}'
    saved = json.loads((d / 'fit_results.json').read_text())['fit']['s']
    ws_path, fit_path = d/'workspace.root', Path(saved['source'])
    hashes = {str(p): sha(p) for p in [ws_path, fit_path]}
    wf, ff = ROOT.TFile.Open(str(ws_path)), ROOT.TFile.Open(str(fit_path))
    w, fit = wf.Get('w'), ff.Get('fit_s')
    apply_fit(w, fit)
    assert fit.status() == 0 and fit.covQual() == 3
    assert fit.floatParsFinal().find('r') and not fit.constPars().find('r')
    metric, bins = inspect_model(w, 's')
    assert abs(metric['exact_deviance']-saved['exact_deviance']) < 1e-7
    pars = list(fit.floatParsFinal())
    vc = fit.covarianceMatrix()
    cov = np.array([[vc[i][j] for j in range(len(pars))] for i in range(len(pars))])
    funcs, keep = {}, []
    for i, era in enumerate(ERAS, 1):
        pdf = add_pdf(w, 's', f'CR_ch{i}')
        for shape, rate in zip(pdf.pdfList(), pdf.coefList()):
            name = shape.GetName()
            proc = 'dy' if name.startswith('local_continuum_') else next(p for p in PROCS[:-1] if name.startswith('shape_'+p+'_'))
            funcs[(era,proc)] = rate
    for proc in PROCS:
        args = ROOT.RooArgList(*[funcs[(era,proc)] for era in ERAS])
        obj = ROOT.RooAddition('combined_'+proc, '', args)
        keep.extend([args,obj]);funcs[('combined',proc)] = obj
    keys = list(funcs)

    def values():
        return np.array([funcs[k].getVal() for k in keys])

    central = values()
    jacobians, derivatives = [], []
    for scale in [1e-3,5e-4]:
        jac = np.zeros((len(keys),len(pars)))
        methods = {}
        for j, par in enumerate(pars):
            var = w.var(par.GetName())
            if not var or par.getError() <= 0:
                continue
            x = par.getVal(); h = scale*par.getError()
            low, high = x-var.getMin(), var.getMax()-x
            if min(low,high) >= h:
                var.setVal(x+h);plus=values()
                var.setVal(x-h);minus=values()
                deriv=(plus-minus)/(2*h);method='central'
            else:
                direction=1 if high>=low else -1
                room=high if direction==1 else low
                h=min(h,room*.25)
                assert h>0,(par.GetName(),x,room)
                var.setVal(x+direction*h);one=values()
                var.setVal(x+direction*2*h);two=values()
                deriv=direction*(-3*central+4*one-two)/(2*h)
                method='forward_second_order' if direction==1 else 'backward_second_order'
            var.setVal(x);jac[:,j]=deriv;methods[par.GetName()]=method
        jacobians.append(jac);derivatives.append(methods)
    propagated=jacobians[0]@cov@jacobians[0].T
    errors=np.sqrt(np.maximum(np.diag(propagated),0))
    half=np.sqrt(np.maximum(np.diag(jacobians[1]@cov@jacobians[1].T),0))
    stability=np.abs(errors-half)/np.maximum(errors,1e-10)
    assert stability.max()<1e-5,stability.max()
    assert np.allclose(values(),central,rtol=1e-12,atol=1e-8)
    rows=[]
    for i,(era,proc) in enumerate(keys):
        rows.append(dict(side=side,mass=mass,era=era,region='CR',process=proc,yield_events=float(central[i]),hesse_error=float(errors[i])))
    baseline=json.loads((OUT/'yields.json').read_text())['results'][side]
    by={(x['era'],x['region'],x['process']):x for x in baseline['rows']}
    for row in rows:
        if row['process'] in PROCS[:2]:
            b=by[(row['era'],'CR',row['process'])]
            row.update(bonly_yield=b['yield_events'],bonly_error=b['hesse_error'],shift_from_bonly=row['yield_events']-b['yield_events'],relative_shift_from_bonly=row['yield_events']/b['yield_events']-1)
    profile=json.loads((P/'runs/20261102_upsilon_c0_profiles'/side/f'm{mass:.1f}'/'profile.json').read_text())
    best=min((x for x in profile['points'] if x.get('profile_point_valid')),key=lambda x:x['exact_deviance'])
    indices=[keys.index(('combined',proc)) for proc in PROCS]
    vc4=propagated[np.ix_(indices,indices)]
    sd4=np.sqrt(np.diag(vc4));corr=vc4/sd4[:,None]/sd4[None,:]
    norm_boundaries=[]
    for par in pars:
        if 'bkgScale_upsilon' in par.GetName():
            var=w.var(par.GetName())
            if min(var.getVal()-var.getMin(),var.getMax()-var.getVal())<1e-4*(var.getMax()-var.getMin()):
                norm_boundaries.append(par.GetName())
    result=dict(side=side,mass=mass,stage='free-r S+B',source_workspace=str(ws_path),source_fit=str(fit_path),source_hashes=hashes,
        status=fit.status(),covQual=fit.covQual(),edm=fit.edm(),r=saved['r'],r_error=saved['error'],mean_boundaries=saved['mean_boundary'],
        normalization_boundaries=norm_boundaries,
        error_method='Local Jacobian times complete free-r native covariance; second-order one-sided derivative near parameter boundaries.',
        caveat='Symmetric Hesse approximations only. Upsilon means at parameter boundaries and weak identifiability can make Gaussian intervals misleading; these are not validated profile intervals.',
        error_step_halving_max_relative=float(stability.max()),one_sided_derivatives={k:v for k,v in derivatives[0].items() if v!='central'},
        CR_range=[w.var('mass_CR').getMin(),w.var('mass_CR').getMax()],rows=rows,
        combined_order=PROCS,combined_covariance=vc4.tolist(),combined_correlation=corr.tolist(),
        covariance_order=[list(k) for k in keys],yield_covariance=propagated.tolist(),
        profile_best=dict(r=best['r'],deviance=best['exact_deviance'],covQual=best['covQual'],source=best['source'],CR_yields=best['selected_yields_by_region']['CR'],
                          source_is_fixed_r=True,covariance_not_used=True),**metric)
    ff.Close();wf.Close()
    assert hashes=={p:sha(p) for p in hashes}
    result['source_hashes_unchanged']=True
    return result


def main():
    data=dict(description='CR-only four-era Upsilon yields from original free-r S+B fits, compared with original B-only fits; no new minimizations.',results=[extract(s,m) for s in ['low','high'] for m in [9.5,10.0]])
    dump(OUT/'sb_comparison.json',data)
    rows=[r for q in data['results'] for r in q['rows']]
    fields=list(dict.fromkeys(k for r in rows for k in r))
    with (OUT/'sb_comparison.csv').open('w',newline='') as f:
        wr=csv.DictWriter(f,fieldnames=fields);wr.writeheader();wr.writerows(rows)
    for q in data['results']:
        print(q['side'],q['mass'],'r',q['r'],'+/-',q['r_error'],'status/cov',q['status'],q['covQual'],'boundary',q['mean_boundaries'])
        for r in q['rows']:
            if r['era']=='combined':print(r)
        print('correlation',q['combined_correlation'],'stability',q['error_step_halving_max_relative'])


if __name__=='__main__':
    main()
