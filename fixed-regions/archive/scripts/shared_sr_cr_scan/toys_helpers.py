"""ROOT-side helpers for the global shared-SR Cheb7 -> envelope6 toy experiment."""
from pathlib import Path
import json, math, os, resource, signal, subprocess, time
import ROOT
import numpy as np

ERAS = ['2022', '2022EE', '2023', '2023BPix']


def process_cpu_usage():
    own=resource.getrusage(resource.RUSAGE_SELF)
    children=resource.getrusage(resource.RUSAGE_CHILDREN)
    return dict(worker_seconds=own.ru_utime+own.ru_stime,
                native_children_seconds=children.ru_utime+children.ru_stime,
                total_seconds=own.ru_utime+own.ru_stime+children.ru_utime+children.ru_stime)


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


def channel_layout(workspace):
    mc = workspace.obj('ModelConfig')
    pdf = mc.GetPdf()
    regions = [('SRL', 'mass_SRL'), ('CR', 'mass_CR'), ('SRH', 'mass_SRH')]
    result = []
    for region, observable in regions:
        for index in range(1, 5):
            channel = f'{region}_ch{index}'
            cp = pdf.getPdf(channel)
            if cp:
                result.append((channel, observable))
    if not result:
        raise ValueError('No expected shared-SR channels')
    return result


def exact_deviance(workspace, signal=True):
    """Unshifted extended-bin deviance plus Gaussian auxiliary penalties.

    Combine evaluates the native binned likelihood at bin centres. The yield
    integrals themselves use the PDF's exact analytic normalizations. Equal
    discrete-model penalties are omitted from minimum comparisons only.
    """
    mc = workspace.obj('ModelConfig' if signal else 'ModelConfig_bonly')
    dataset = workspace.data('toy_data') or workspace.data('data_obs')
    value, minimum_bin = 0., math.inf
    for channel, observable in channel_layout(workspace):
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
    active_regions = set(c.split('_')[0] for c, _ in channel_layout(workspace))
    for era in ERAS:
        if 'CR' in active_regions:
            coeff = [1.] + [workspace.var(f'local_c{j}_{era}_CR').getVal() for j in range(1, 5)]
            poly = np.polynomial.Chebyshev(coeff)
            roots = poly.deriv().roots()
            points = [-1., 1.] + [float(q.real) for q in roots if abs(q.imag) < 1e-9 and -1 < q.real < 1]
            minima[f'CR_{era}'] = float(min(poly(q) for q in points))
        if not ({'SRL', 'SRH'} & active_regions):
            continue
        cat = workspace.cat('pdfindex_' + era)
        if cat and cat.getIndex() == 1 and workspace.var(f'local_b1_{era}'):
            co = [1.]
            j = 1
            while workspace.var(f'local_b{j}_{era}'):
                co.append(workspace.var(f'local_b{j}_{era}').getVal()); j += 1
            minima[f'SR_{era}'] = min(co)
        else:
            co = [1.]
            j = 1
            while workspace.var(f'local_c{j}_{era}'):
                co.append(workspace.var(f'local_c{j}_{era}').getVal()); j += 1
            poly = np.polynomial.Chebyshev(co)
            roots = poly.deriv().roots()
            # Builder imposes positivity on full reference6–12, including the gap.
            # Exact stationary-point audit strengthens its construction grid.
            points = [-1., 1.] + [float(q.real) for q in roots if abs(q.imag)<1e-9 and -1<q.real<1]
            for region in ['SRL', 'SRH']:
                if region not in active_regions:
                    continue
                x = workspace.var('mass_' + region)
                lo, hi = (x.getMin() - 9) / 3, (x.getMax() - 9) / 3
                points += [lo, hi] + [float(q.real) for q in roots if abs(q.imag) < 1e-9 and lo < q.real < hi]
            minima[f'SR_{era}'] = float(min(poly(q) for q in points))
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
    discrete = w.genobj('discreteParams')
    if discrete and discrete.getSize() == 4 and len(result['pdf_indices']) != 4:
        source.Close(); handle.Close()
        return dict(readable=False, failure='missing_envelope_indices_in_saved_fit')
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
    result['edm_valid'] = bool(math.isfinite(result['edm']) and 0 <= result['edm'] < .01)
    result['numerical_validity_version'] = 'status0_edm_lt_0p01_v1'
    result['fit_valid'] = (result['status'] == 0 and result['edm_valid'] and result['positive'] and result['deviance_with_aux'] is not None and not result.get('r_boundary', False))
    source.Close()
    handle.Close()
    return result
