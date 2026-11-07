"""ROOT-side helpers for the upper-SR Cheb5 -> envelope3 toy experiment."""
from pathlib import Path
import json, math, os, signal, subprocess, time
import ROOT
import numpy as np

ERAS = ['2022', '2022EE', '2023', '2023BPix']


def atomic_json(path, value):
    path = Path(path)
    temp = path.with_suffix(path.suffix + '.tmp')
    def clean(item):
        if isinstance(item, dict):
            return {k: clean(v) for k, v in item.items()}
        if isinstance(item, (tuple, list)):
            return [clean(v) for v in item]
        if isinstance(item, float) and not math.isfinite(item):
            return None
        return item
    temp.write_text(json.dumps(clean(value), indent=2, allow_nan=False))
    temp.replace(path)


def execute(command, directory, timeout=180):
    """Retain every receipt and kill a whole subprocess group on timeout."""
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    receipt = directory / 'receipt.json'
    if receipt.exists():
        return json.loads(receipt.read_text())
    atomic_json(directory / 'command.json', command)
    start = time.time()
    with (directory / 'combine.log').open('w') as log:
        proc = subprocess.Popen(command, cwd=directory, stdout=log,
                                stderr=subprocess.STDOUT, start_new_session=True)
        try:
            code = proc.wait(timeout=timeout)
        except subprocess.TimeoutExpired:
            os.killpg(proc.pid, signal.SIGTERM)
            try:
                proc.wait(timeout=5)
            except subprocess.TimeoutExpired:
                os.killpg(proc.pid, signal.SIGKILL)
                proc.wait()
            code = 124
    value = dict(exit_code=code, seconds=time.time() - start, timeout_seconds=timeout)
    atomic_json(receipt, value)
    return value


def apply_globals(workspace, snapshot):
    globals_ = workspace.obj('ModelConfig').GetGlobalObservables()
    values = {}
    for observable in globals_:
        value = snapshot.find(observable.GetName())
        if not value:
            raise ValueError('Missing auxiliary observation ' + observable.GetName())
        observable.setVal(value.getVal())
        values[observable.GetName()] = value.getVal()
    return values


def apply_result(workspace, fit):
    for item in list(fit.floatParsFinal()) + list(fit.constPars()):
        name = item.GetName()
        if item.InheritsFrom('RooAbsCategory'):
            if workspace.cat(name):
                workspace.cat(name).setIndex(item.getIndex())
        elif workspace.var(name):
            workspace.var(name).setVal(item.getVal())


def exact_deviance(workspace, signal=True):
    """Extended-bin deviance plus original unit-Gaussian auxiliary penalties.

    All comparisons use these unshifted values, never offset minNll values.
    Gaussian global-observable pairs are verified by their model constraints.
    """
    mc = workspace.obj('ModelConfig' if signal else 'ModelConfig_bonly')
    dataset = workspace.data('toy_data') or workspace.data('data_obs')
    regions = [('SR_', 'mass'), ('CR_', 'mass_CR')] if workspace.var('mass_CR') else [('', 'mass')]
    value, minimum_bin = 0., math.inf
    for prefix, observable in regions:
        for index in range(1, 5):
            channel = prefix + f'ch{index}'
            cp = mc.GetPdf().getPdf(channel)
            pdf = next(p for p in cp.getComponents() if p.InheritsFrom('RooAddPdf'))
            x = workspace.var(observable)
            obs = ROOT.RooArgSet(x)
            bw = (x.getMax() - x.getMin()) / x.getBins()
            total = pdf.expectedEvents(obs)
            data = dataset.reduce(ROOT.RooFit.Cut(f'CMS_channel==CMS_channel::{channel}'))
            ns, term = 0., 0.
            for k in range(data.numEntries()):
                row = data.get(k)
                x.setVal(row.getRealValue(observable))
                n = data.weight()
                mu = pdf.getVal(obs) * total * bw
                minimum_bin = min(minimum_bin, mu)
                if mu <= 0 or not math.isfinite(mu):
                    return None, minimum_bin
                ns += n
                if n:
                    term += n * math.log(n / mu)
            value += 2 * (total - ns + term)
    for glob in mc.GetGlobalObservables():
        name = glob.GetName()
        if not name.endswith('_In') or not workspace.var(name[:-3]):
            raise ValueError('Unknown auxiliary constraint ' + name)
        value += (workspace.var(name[:-3]).getVal() - glob.getVal()) ** 2
    return value, minimum_bin


def continuum_positivity(workspace):
    minima = {}
    for era in ERAS:
        for region in ['SR', 'CR'] if workspace.var('mass_CR') else ['CRonly']:
            suffix = '_CR' if region == 'CR' else ''
            degree = 3 if region == 'SR' else 4
            cat = workspace.cat('pdfindex_' + era) if region == 'SR' else None
            if cat and cat.getIndex() == 1:
                # Positive Bernstein coefficients ensure positive shape everywhere.
                minima[f'{region}_{era}'] = min([1.] + [workspace.var(f'local_b{j}_{era}').getVal() for j in range(1, 4)])
            else:
                coeff = [1.] + [workspace.var(f'local_c{j}_{era}{suffix}').getVal() for j in range(1, degree + 1)]
                poly = np.polynomial.Chebyshev(coeff)
                roots = poly.deriv().roots()
                points = [-1., 1.] + [float(q.real) for q in roots if abs(q.imag) < 1e-9 and -1 < q.real < 1]
                minima[f'{region}_{era}'] = float(min(poly(q) for q in points))
    return minima


def read_fit(directory, signal, interval_method):
    directory = Path(directory)
    filename = directory / 'fitDiagnostics_upper.root'
    if not filename.exists():
        return dict(readable=False, failure='missing_fit_result')
    handle = ROOT.TFile.Open(str(filename))
    if not handle or handle.IsZombie():
        return dict(readable=False, failure='unreadable_fit_result')
    fit = handle.Get('fit_s' if signal else 'fit_b')
    if not fit:
        handle.Close()
        return dict(readable=False, failure='missing_fit_s' if signal else 'missing_fit_b')
    params = []
    for item in fit.floatParsFinal():
        params.append(dict(name=item.GetName(), value=item.getVal(), error=item.getError(),
                           low=item.getMin(), high=item.getMax()))
    result = dict(readable=True, status=fit.status(), covQual=fit.covQual(),
                  edm=fit.edm(), minNll_offset=fit.minNll(), parameters=params,
                  pdf_indices={v.GetName(): v.getIndex() for v in fit.constPars() if v.InheritsFrom('RooAbsCategory') and v.GetName().startswith('pdfindex_')})
    source = ROOT.TFile.Open(str(directory / 'workspace.root'))
    w = source.Get('w')
    apply_result(w, fit)
    if not signal:
        w.var('r').setVal(0.)
    result['deviance_with_aux'], result['minimum_bin_expectation'] = exact_deviance(w, signal)
    result['continuum_minima'] = continuum_positivity(w)
    result['positive'] = (result['minimum_bin_expectation'] > 0 and min(result['continuum_minima'].values()) > 0)
    result['boundary_parameters'] = []
    for a in params:
        if not all(math.isfinite(a[k]) for k in ['value', 'low', 'high']):
            continue
        if abs(a['value'] - a['low']) < 1e-4 * max(1., a['high'] - a['low']) or abs(a['value'] - a['high']) < 1e-4 * max(1., a['high'] - a['low']):
            result['boundary_parameters'].append(a['name'])
    if signal:
        r = fit.floatParsFinal().find('r')
        result.update(r=r.getVal(), error=r.getError(), interval_method=interval_method)
        if interval_method == 'profile68':
            low = r.getMin('err68') if r.hasRange('err68') else r.getVal() + r.getAsymErrorLo()
            high = r.getMax('err68') if r.hasRange('err68') else r.getVal() + r.getAsymErrorHi()
        else:
            low, high = r.getVal() - r.getError(), r.getVal() + r.getError()
        result.update(interval_low=low, interval_high=high,
                      error_lo=r.getVal() - low, error_hi=high - r.getVal(),
                      r_boundary='r' in result['boundary_parameters'],
                      interval_boundary=(low <= r.getMin() + 1e-6 or high >= r.getMax() - 1e-6))
        result['interval_valid'] = (all(math.isfinite(v) for v in [low, high, r.getVal()]) and low < r.getVal() < high and not result['interval_boundary'])
    result['fit_valid'] = (result['status'] == 0 and result['positive'] and result['deviance_with_aux'] is not None and not result.get('r_boundary', False))
    source.Close()
    handle.Close()
    return result
