"""Read a native conditional FitDiagnostics fit; never minimize or change inputs.

FitDiagnostics --customStartingPoint --skipSBFit saves the fixed-r fit to
fit_b although the minimized PDF is ModelConfig's full signal+background PDF.
"""
from pathlib import Path
import argparse, math, sys
import ROOT

P = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(P / 'scripts/separate_sr_c0_limits'))
from fit_helpers import apply_fit, inspect_model, polynomial_minimum, continuity, dump, ERAS

ROOT.gROOT.SetBatch(True)
ROOT.gSystem.Load('libHiggsAnalysisCombinedLimit')
ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.WARNING)
PROCESSES = ('dy', 'upsilon1s', 'upsilon2s', 'Zd')


def process(name):
    if name.startswith('local_continuum_'):
        return 'dy'
    return next(p for p in PROCESSES[1:] if name.startswith('shape_' + p + '_'))


def extract(workspace, fit_path, fixed_r, output=None):
    workspace, fit_path = Path(workspace), Path(fit_path)
    result = dict(readable=False, workspace=str(workspace), source=str(fit_path),
                  requested_r=float(fixed_r), point_valid=False, numerical_ok=False,
                  stored_fit_name='fit_b', evaluated_model='ModelConfig', flags=[])
    sf = ff = None
    try:
        sf = ROOT.TFile.Open(str(workspace))
        if not sf or sf.IsZombie() or not sf.Get('w'):
            raise RuntimeError('Missing readable source workspace')
        w = sf.Get('w')
        ff = ROOT.TFile.Open(str(fit_path))
        if not ff or ff.IsZombie() or not ff.Get('fit_b'):
            raise RuntimeError('Missing native conditional fit_b')
        fit = ff.Get('fit_b')
        floats, consts = fit.floatParsFinal(), fit.constPars()
        r = consts.find('r')
        assert r and not floats.find('r'), 'r must be in constPars, absent from floatParsFinal'
        value = float(r.getVal())
        assert abs(value-float(fixed_r)) <= 1e-9*max(1., abs(float(fixed_r))), 'Native r differs from requested grid point'
        assert r.isConstant(), 'Native r is not marked constant'
        mean_names = [f'Ymean_{peak}_{era}' for peak in ('upsilon1s','upsilon2s') for era in ERAS]
        norm_names = [f'CMS_EXO25020_bkgScale_{peak}_{era}' for peak in ('upsilon1s','upsilon2s') for era in ERAS]
        norm_names += [f'local_yield_{era}_CR' for era in ERAS]
        for name in mean_names + norm_names:
            assert floats.find(name) and not consts.find(name), 'Required floating parameter missing: ' + name
            assert not floats.find(name).isConstant(), 'Required parameter marked constant: ' + name
        apply_fit(w, fit)
        w.var('r').setConstant(True)
        parameters = []
        boundary_parameters = []
        gaussian_names={g.GetName()[:-3] for g in w.obj('ModelConfig').GetGlobalObservables()}
        native_unbounded_gaussian_parameters=[]
        for v in floats:
            if not v.InheritsFrom('RooRealVar'):
                continue
            p = dict(name=v.GetName(), value=float(v.getVal()), error=float(v.getError()),
                     low=float(v.getMin()), high=float(v.getMax()))
            assert w.var(v.GetName()), 'Native floating parameter absent from source workspace'
            source_var = w.var(v.GetName())
            p['source_low']=float(source_var.getMin());p['source_high']=float(source_var.getMax())
            native_unbounded=(p['low']==-math.inf and p['high']==math.inf)
            if native_unbounded and p['name'] in gaussian_names:
                # Native Combine removes bounds from some Gaussian-constrained
                # nuisances, also in the previously saved unconditional fits.
                native_unbounded_gaussian_parameters.append(p['name'])
                assert p['source_low']<=p['value']<=p['source_high'], 'Native nuisance value outside source evaluation range: '+p['name']
            else:
                assert abs(p['low']-source_var.getMin()) < 1e-8*max(1.,abs(p['low'])), 'Changed parameter lower bound: '+p['name']
                assert abs(p['high']-source_var.getMax()) < 1e-8*max(1.,abs(p['high'])), 'Changed parameter upper bound: '+p['name']
            if math.isfinite(p['low']) and math.isfinite(p['high']):
                if min(abs(p['value']-p['low']),abs(p['value']-p['high'])) <= 1e-4*max(1.,p['high']-p['low']):
                    boundary_parameters.append(p['name'])
            parameters.append(p)
        metric, bins = inspect_model(w, 's')
        minima = polynomial_minimum(w)
        c0 = continuity(w, 's')
        by_era = {era:{p:0. for p in PROCESSES} for era in ERAS}
        by_region = {region:{p:0. for p in PROCESSES} for region in ('SR','CR')}
        by_channel = {}
        for channel,b in bins.items():
            amounts={p:0. for p in PROCESSES}
            for name,y in b['component_yields'].items():
                amounts[process(name)]+=float(y)
            for p,y in amounts.items():
                by_era[b['era']][p]+=y
                by_region[b['region']][p]+=y
            by_channel[channel]=dict(era=b['era'],region=b['region'],yields=amounts,
                expected_events=b['expected_events'],observed_events=b['observed_events'],deviance=b['deviance'])
        selected={p:sum(v[p] for v in by_era.values()) for p in PROCESSES}
        edm=float(fit.edm());status=int(fit.status());cov=int(fit.covQual())
        valid_edm=math.isfinite(edm) and 0 <= edm < .01
        valid=status==0 and valid_edm and min(minima.values())>0 and metric['minimum_bin']>0 and math.isfinite(metric['exact_deviance'])
        c0passed=c0['maximum_relative_difference']<1e-8
        if status != 0:result['flags'].append('native_fit_status_nonzero')
        if not valid_edm:result['flags'].append('invalid_EDM')
        if cov != 3:result['flags'].append('covariance_warning')
        if min(minima.values())<=0 or metric['minimum_bin']<=0:result['flags'].append('nonpositive_model')
        if not math.isfinite(metric['exact_deviance']):result['flags'].append('nonfinite_deviance')
        if not c0passed:result['flags'].append('C0_mismatch')
        result.update(readable=True,r=value,r_fixed=True,status=status,covQual=cov,edm=edm,
            status_history=[dict(label=fit.statusLabelHistory(i),code=int(fit.statusCodeHistory(i))) for i in range(fit.numStatusHistory())],
            num_invalid_NLL=int(fit.numInvalidNLL()),
            point_valid=bool(valid and c0passed),numerical_ok=bool(valid and c0passed and cov==3),
            parameters=parameters,initialization={q['name']:q['value'] for q in parameters},
            floating_means={n:float(w.var(n).getVal()) for n in mean_names},
            mean_boundaries=[n for n in mean_names if n in boundary_parameters],
            normalization_boundaries=[n for n in norm_names if n in boundary_parameters],
            parameter_boundaries=boundary_parameters,continuum_minima=minima,
            native_unbounded_gaussian_parameters=native_unbounded_gaussian_parameters,
            continuity_passed=c0passed,C0_max_relative_difference=c0['maximum_relative_difference'],
            selected_yields=selected,selected_yields_by_era=by_era,selected_yields_by_region=by_region,
            channels=by_channel,**metric)
    except Exception as error:
        result['failure']=type(error).__name__+': '+str(error)
        result['flags'].append('extraction_failed')
    finally:
        if ff:ff.Close()
        if sf:sf.Close()
    if output:dump(output,result)
    return result


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('workspace',type=Path)
    parser.add_argument('fit_path',type=Path)
    parser.add_argument('fixed_r',type=float)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    q=extract(args.workspace,args.fit_path,args.fixed_r,args.output)
    print(q.get('r',args.fixed_r),q.get('status'),q.get('edm'),q.get('exact_deviance'),q['flags'],flush=True)
    raise SystemExit(0 if q['readable'] else 1)
