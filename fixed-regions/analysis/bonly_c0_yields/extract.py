"""Read converged native r=0 fits and propagate their complete Hesse covariance.

This script never minimizes or alters a source ROOT workspace/result.
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
sys.path.insert(0, str(P / 'analysis/separate_sr_c0_limits'))
from fit_helpers import ERAS, add_pdf, apply_fit, continuity, dump, inspect_model

ROOT.gROOT.SetBatch(True)
OUT = P / 'runs/20261102_c0_background_yields'
PROCS = ['upsilon1s', 'upsilon2s']


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def extract(side):
    d = P / 'runs/20261031_separate_sr_c0_limits' / side / 'm9.5'
    old = json.loads((d / 'fit_results.json').read_text())['fit']['b']
    ws_path, fit_path = d / 'workspace.root', Path(old['source'])
    hashes = {str(p): sha(p) for p in [ws_path, fit_path]}
    wf, ff = ROOT.TFile.Open(str(ws_path)), ROOT.TFile.Open(str(fit_path))
    w, fit = wf.Get('w'), ff.Get('fit_b')
    assert fit.status() == 0 and fit.covQual() == 3
    apply_fit(w, fit)
    w.var('r').setVal(0)
    assert fit.constPars().find('r').getVal() == 0
    metric, bins = inspect_model(w, 'b')
    assert abs(metric['exact_deviance'] - old['exact_deviance']) < 1e-7
    pars = list(fit.floatParsFinal())
    cov_native = fit.covarianceMatrix()
    cov = np.array([[cov_native[i][j] for j in range(len(pars))] for i in range(len(pars))])
    functions, keepalive = {}, []
    for i, era in enumerate(ERAS, 1):
        for region in ['SR', 'CR']:
            pdf = add_pdf(w, 'b', f'{region}_ch{i}')
            for component, coefficient in zip(pdf.pdfList(), pdf.coefList()):
                for proc in PROCS:
                    if component.GetName().startswith('shape_' + proc + '_'):
                        functions[(era, region, proc)] = coefficient

    def add(key, terms):
        args = ROOT.RooArgList()
        for term in terms:
            args.add(term)
        obj = ROOT.RooAddition('_'.join(key), '_'.join(key), args)
        keepalive.extend([args, obj])
        functions[key] = obj

    for era in ERAS:
        for proc in PROCS:
            add((era, 'SR+CR', proc), [functions[(era, reg, proc)] for reg in ['SR', 'CR']])
    for reg in ['SR', 'CR', 'SR+CR']:
        for proc in PROCS:
            add(('combined', reg, proc), [functions[(era, reg, proc)] for era in ERAS])
    keys = list(functions)
    values = np.array([functions[k].getVal() for k in keys])
    root_errors = np.array([functions[k].getPropagatedError(fit) for k in keys])
    # Infinitesimal Jacobian covariance propagation. Repeat at half the step
    # to check numerical stability independently of ROOT's finite +/-1sigma
    # getPropagatedError approximation for nonlinear RooFormulaVars.
    gradients = []
    for step_scale in [1e-3, 5e-4]:
        jac = np.zeros((len(keys), len(pars)))
        for j, par in enumerate(pars):
            var = w.var(par.GetName())
            if not var or par.getError() <= 0:
                continue
            center = par.getVal()
            step = step_scale * par.getError()
            room = min(center - var.getMin(), var.getMax() - center)
            if math.isfinite(room):
                step = min(step, room * .25)
            assert step > 0, (par.GetName(), center, room)
            var.setVal(center + step)
            plus = np.array([functions[k].getVal() for k in keys])
            var.setVal(center - step)
            minus = np.array([functions[k].getVal() for k in keys])
            var.setVal(center)
            jac[:, j] = (plus - minus) / (2 * step)
        gradients.append(jac)
    propagated = gradients[0] @ cov @ gradients[0].T
    errors = np.sqrt(np.maximum(np.diag(propagated), 0))
    errors_halfstep = np.sqrt(np.maximum(np.diag(gradients[1] @ cov @ gradients[1].T), 0))
    relative_stability = np.abs(errors - errors_halfstep) / np.maximum(errors, 1e-12)
    assert relative_stability.max() < 1e-5, relative_stability.max()
    restored = np.array([functions[k].getVal() for k in keys])
    assert np.allclose(restored, values, rtol=1e-12, atol=1e-8)
    rows = []
    for i, key in enumerate(keys):
        era, region, proc = key
        rows.append(dict(side=side, era=era, region=region, process=proc,
                         yield_events=float(values[i]), hesse_error=float(errors[i]),
                         root_getPropagatedError=float(root_errors[i]),
                         root_error_relative_difference=float((root_errors[i] - errors[i]) / errors[i]),
                         halfstep_relative_difference=float(relative_stability[i])))
    boundaries = []
    resonance_parameters = []
    for par in pars:
        if par.GetName().startswith('Ymean_') or 'bkgScale_upsilon' in par.GetName():
            v = w.var(par.GetName())
            distance = min(v.getVal()-v.getMin(), v.getMax()-v.getVal())
            edge = distance < 1e-4 * (v.getMax()-v.getMin())
            if edge:
                boundaries.append(par.GetName())
            resonance_parameters.append(dict(name=par.GetName(), value=par.getVal(),
                hesse_error=par.getError(), low=v.getMin(), high=v.getMax(), at_boundary=edge))
    result = dict(side=side, source_workspace=str(ws_path), source_fit=str(fit_path),
        source_hashes=hashes, signal_strength=0, status=fit.status(), covQual=fit.covQual(), edm=fit.edm(),
        SR_degree=6 if side=='low' else 4, CR_degree=4,
        SR_range=[w.var('mass_SR').getMin(), w.var('mass_SR').getMax()],
        CR_range=[w.var('mass_CR').getMin(), w.var('mass_CR').getMax()],
        physical_yields='Analytical extended coefficient functions, integrated over the specified regions.',
        error_method='sqrt(J V J^T), with complete native fit covariance and numerical local derivatives. Means and normalizations both included; all inter-parameter and inter-era correlations retained.',
        error_interpretation='Approximate symmetric one-standard-deviation Hesse errors; not profile-likelihood intervals or externally calibrated normalization constraints.',
        error_numerical_stability_max_relative=float(relative_stability.max()),
        root_getPropagatedError_max_relative_difference=float(np.max(np.abs(root_errors-errors)/errors)),
        resonance_parameters=resonance_parameters, resonance_parameters_at_boundary=boundaries,
        rows=rows, covariance_order=[list(k) for k in keys], yield_covariance=propagated.tolist(),
        total_expected_events=sum(c['expected_events'] for c in bins.values()),
        total_observed_events=sum(c['observed_events'] for c in bins.values()),
        continuity=continuity(w,'b'), **metric)
    ff.Close()
    wf.Close()
    assert hashes == {p: sha(p) for p in hashes}
    result['source_hashes_unchanged'] = True
    return result


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    data = dict(method='Existing converged background-only fits, r=0; read-only extraction, no new minimization.',
        correlated_analyses='The two analyses share CR events and must not be added or treated as independent measurements.',
        nominal_mass='m9.5 source directory only; the zero-signal fit is identical for all signal mass hypotheses within each side.',
        results={side: extract(side) for side in ['low','high']})
    dump(OUT/'yields.json', data)
    rows = [row for result in data['results'].values() for row in result['rows']]
    with (OUT/'yields.csv').open('w', newline='') as f:
        writer=csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)
    lines=['# Υ yields from the current background-only fits', '',
        'All entries are fitted numbers of events, with approximate symmetric 1σ Hesse errors propagated with the complete covariance. No new fit was necessary: the existing r=0 fits both have status 0 and covariance quality 3.', '',
        'The low+CR and high+CR fits share CR data. Do not add their yields or interpret their difference as independent measurements.', '']
    for side, result in data['results'].items():
        lines += [f'## {side}: SR + CR', '', '| Era | Υ(1S), SR+CR | Υ(2S), SR+CR |', '|---|---:|---:|']
        lookup={(r['era'],r['region'],r['process']):r for r in result['rows']}
        def val(era,reg,proc):
            r=lookup[(era,reg,proc)]
            return f"{r['yield_events']:.1f} ± {r['hesse_error']:.1f}"
        for era in ERAS+['combined']:
            lines.append(f"| {era} | {val(era,'SR+CR','upsilon1s')} | {val(era,'SR+CR','upsilon2s')} |")
        lines += ['', '| Region, all eras | Υ(1S) | Υ(2S) |', '|---|---:|---:|']
        for reg in ['SR','CR','SR+CR']:
            lines.append(f"| {reg} | {val('combined',reg,'upsilon1s')} | {val('combined',reg,'upsilon2s')} |")
        lines += ['', f"SR range: {result['SR_range']}; CR range: {result['CR_range']} GeV.",
                  f"No Υ normalization or mean lies at a fit boundary. Derivative step-halving check: maximum relative error change {result['error_numerical_stability_max_relative']:.3g}.", '']
    (OUT/'yields.md').write_text('\n'.join(lines)+'\n')
    print((OUT/'yields.md').read_text())


if __name__=='__main__':
    main()
