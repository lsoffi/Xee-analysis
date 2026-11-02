"""Conditional frequentist toys for the shared low/high SR study.

Bias r=0 uses one common 50-dataset ensemble across masses (paired fits).
Injected ensembles are independent by mass. GoF uses 100 separate common
B-only datasets from the fitted envelope model, never the alternate bias truth.
"""
from pathlib import Path
import argparse, fcntl, hashlib, json, math, struct, sys, time, shutil
import ROOT
from toys_helpers import atomic_json, execute, read_fit, apply_result, channel_layout
from toys_worker import set_observable_values, fit_source

ROOT.gROOT.SetBatch(True)
ROOT.gSystem.Load('libHiggsAnalysisCombinedLimit')
P = Path(__file__).resolve().parents[2]
R = P / 'runs/20261031_shared_sr_cr_scan'
C = str(P / 'bin/combine')
MASSES = [f'{k / 10:.1f}' for k in range(60, 121)]
CHANNELS = ','.join(f'{reg}_ch{i}' for reg in ['SRL', 'CR', 'SRH'] for i in range(1, 5))


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load(path):
    return json.loads(Path(path).read_text())


def prepare_control(source):
    directory = R / 'control'
    directory.mkdir(parents=True, exist_ok=True)
    with (directory / 'prepare.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        if (directory / 'reference.json').exists():
            return
        from model import make_cr_workspace
        make_cr_workspace(source, directory / 'workspace.root', means_floating=True)
        handle = ROOT.TFile.Open(str(directory / 'workspace.root'))
        w = handle.Get('w')
        params = [dict(name=v.GetName(), value=v.getVal()) for v in w.allVars()
                  if v.GetName().startswith(('local_', 'Ymean_', 'CMS_EXO25020_bkgScale_'))]
        assert len([p for p in params if p['name'].startswith('Ymean_')]) == 8
        atomic_json(directory / 'reference.json', dict(parameters=params, source=str(source),
                    source_sha256=sha(source), means_floating=True,
                    channels=[c for c, _ in channel_layout(w)]))
        handle.Close()


def prepare_generator(mass='7.5'):
    from model import build
    directory = R / 'bias/generator'
    directory.mkdir(parents=True, exist_ok=True)
    with (directory / 'prepare.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        if (directory / 'results.json').exists():
            return load(directory / 'results.json')
        attempts = []
        for variant in range(3):
            wd = directory / f'start{variant}'
            wd.mkdir(exist_ok=True)
            path = wd / 'workspace.root'
            build(mass, degree=7, family='cheb', output=path, freeze_means=True)
            f = ROOT.TFile.Open(str(path)); w = f.Get('w')
            if variant:
                for era in ['2022', '2022EE', '2023', '2023BPix']:
                    for j in [6, 7]:
                        v = w.var(f'local_c{j}_{era}')
                        v.setVal(v.getVal() + (1 if variant == 1 else -1) * .0005)
            set_observable_values(w); w.writeToFile(str(path)); f.Close()
            command = [C, '-M', 'FitDiagnostics', str(path), '-m', mass, '-n', '_upper',
                       '--skipSBFit', '--minos', 'none', '--robustFit', '0', '--keepFailures',
                       '--cminDefaultMinimizerStrategy', '2', '--cminDefaultMinimizerTolerance', '.0001',
                       '--cminDefaultMinimizerPrecision', '1e-12', '--X-rtd',
                       'MINIMIZER_freezeDisassociatedParams', '-v', '0']
            receipt = execute(command, wd, 900)
            q = read_fit(wd, False, 'none')
            q.update(variant=variant, directory=str(wd), receipt=receipt)
            q['valid'] = bool(q.get('fit_valid') and q.get('covQual') == 3 and receipt['exit_code'] == 0)
            attempts.append(q)
            atomic_json(wd / 'results.json', q)
            print('GENERATOR', variant, q.get('status'), q.get('covQual'), q.get('deviance_with_aux'), flush=True)
        valid = [q for q in attempts if q['valid']]
        if len(valid) < 2:
            atomic_json(directory / 'attempts.json', attempts)
            raise RuntimeError('Cheb7 generator requires at least two valid B-only starts')
        best = min(valid, key=lambda q: q['deviance_with_aux'])
        f = ROOT.TFile.Open(str(Path(best['directory']) / 'workspace.root')); w = f.Get('w')
        ff = ROOT.TFile.Open(str(Path(best['directory']) / 'fitDiagnostics_upper.root'))
        apply_result(w, ff.Get('fit_b'))
        means = [dict(name=v.GetName(), value=v.getVal(), constant=bool(v.isConstant()))
                 for v in w.allVars() if v.GetName().startswith('Ymean_')]
        assert len(means) == 8 and all(v['constant'] for v in means)
        w.var('r').setVal(0); set_observable_values(w)
        w.writeToFile(str(directory / 'workspace_postfit.root'))
        ff.Close(); f.Close()
        result = dict(best, mass=float(mass), generator_degree=7, fit_degree=6, CR_degree=4,
                      means=means, attempts=attempts,
                      max_valid_deviance_difference=max(q['deviance_with_aux'] for q in valid)-min(q['deviance_with_aux'] for q in valid),
                      interpretation='Alternative global Chebyshev7 shared SR intensity, CR Cheb4 independent, B-only observed fit; fixed means from CR-only observed fit')
        atomic_json(directory / 'results.json', result)
        prepare_control(directory / 'workspace_postfit.root')
        return result


def expected_source(mass):
    return R / f'observed/m{float(mass):.1f}/results.json'


def expected_anchor(mass):
    return R / f'bias/anchors/m{float(mass):.1f}.json'


def expected_value(mass):
    path = expected_source(mass)
    value = float(load(path)['limits']['expected_50'])
    if not math.isfinite(value) or value <= 0:
        raise ValueError('Expected median must be positive and finite: ' + str(path))
    anchor = expected_anchor(mass); anchor.parent.mkdir(parents=True, exist_ok=True)
    with anchor.with_suffix('.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        if anchor.exists():
            old = load(anchor)
            if abs(old['expected_median'] - value) > 1e-12 * max(1., abs(value)):
                raise RuntimeError('Expected median changed after toy anchor was fixed: ' + str(path))
            return float(old['expected_median'])
        atomic_json(anchor, dict(mass=float(mass), expected_median=value,
                    source=str(path), source_sha256_at_capture=sha(path),
                    field='limits.expected_50',
                    interpretation='Immutable injection/initialization anchor; source results may later gain unrelated GoF or audit fields'))
    return value


def audit_toys(path, ntoys):
    handle = ROOT.TFile.Open(str(path))
    if not handle or handle.IsZombie() or handle.TestBit(ROOT.TFile.kRecovered):
        raise RuntimeError('Unusable generated ROOT file ' + str(path))
    rows = []
    for i in range(1, ntoys + 1):
        toy = handle.Get(f'toys/toy_{i}'); snapshot = handle.Get(f'toys/toy_{i}_snapshot')
        assert toy and snapshot
        digest = hashlib.sha256(); counts = {}
        for k in range(toy.numEntries()):
            row = toy.get(k); channel = row.find('CMS_channel').getLabel(); n = toy.weight()
            assert math.isfinite(n) and n >= 0 and abs(n-round(n)) < 1e-8
            counts[channel] = counts.get(channel, 0) + n
            digest.update(channel.encode())
            digest.update(struct.pack('dddd', row.getRealValue('mass_SRL'), row.getRealValue('mass_CR'), row.getRealValue('mass_SRH'), n))
        assert len(counts) == 12 and toy.numEntries() == 1800, (len(counts), toy.numEntries())
        rows.append(dict(toy=i, events=toy.sumEntries(), entries=toy.numEntries(), counts=counts,
                         data_sha256=digest.hexdigest(), snapshot_parameters=snapshot.getSize()))
    handle.Close()
    assert len({q['data_sha256'] for q in rows}) == ntoys
    return rows


def generate_file(directory, source, seed, ntoys, rtrue, mass, extra):
    directory.mkdir(parents=True, exist_ok=True)
    if (directory / 'truth.json').exists():
        old = load(directory / 'truth.json')
        assert sha(old['toy_file']) == old['toy_file_sha256']
        assert sha(source) == old['source_sha256']
        return old
    f = ROOT.TFile.Open(str(source)); w = f.Get('w')
    w.var('r').setRange(-1e6, 1e6); w.var('r').setVal(rtrue)
    for v in w.allVars():
        if v.GetName().startswith(('local_', 'Ymean_', 'CMS_EXO25020_bkgScale_')):
            v.setConstant(True)
    for cat in w.allCats():
        if cat.GetName().startswith('pdfindex_'):
            cat.setConstant(True)
    if w.genobj('discreteParams'):
        w.genobj('discreteParams').removeAll()
    set_observable_values(w); w.saveSnapshot('generation', w.allVars())
    generation_workspace = directory / 'generator.root'
    w.writeToFile(str(generation_workspace)); f.Close()
    command = [C, '-M', 'GenerateOnly', str(generation_workspace), '--snapshotName', 'generation',
               '--saveToys', '--toysFrequentist', '--bypassFrequentistFit', '--expectSignal', str(rtrue),
               '--genBinnedChannels', CHANNELS, '-t', str(ntoys), '-s', str(seed), '-m', mass, '-n', '_generated']
    receipt = execute(command, directory, 900)
    assert receipt['exit_code'] == 0, receipt
    files = list(directory.glob('higgsCombine*.root')); assert len(files) == 1
    rows = audit_toys(files[0], ntoys)
    truth = dict(mass=float(mass), rtrue=rtrue, seed=seed, ntoys=ntoys, source=str(source),
                 source_sha256=sha(source), generator_workspace_sha256=sha(generation_workspace),
                 toy_file=str(files[0]), toy_file_sha256=sha(files[0]), toy_audit=rows,
                 generation='Frequentist binned Poisson data and auxiliary observations; conditional fitted truth',
                 signal_in_all_regions=True, regions={'SRL':[6,9.093333333333334],
                 'CR':[9.093333333333334,10.506666666666668], 'SRH':[10.506666666666668,12]})
    truth.update(extra); atomic_json(directory / 'truth.json', truth)
    return truth


def generate_bias(mass, label):
    mass = f'{float(mass):.1f}'
    directory = R / f'bias/generation/m{mass}/{label}'
    directory.mkdir(parents=True, exist_ok=True)
    with (directory / 'generate.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        if (directory / 'truth.json').exists():
            return load(directory / 'truth.json')
        expected = expected_value(mass)
        generator = prepare_generator()
        fit_hash = sha(R / 'bias/generator/results.json')
        if label == 'r0':
            common = R / 'bias/generation/common_r0'
            common.mkdir(exist_ok=True)
            with (common / 'generate.lock').open('a') as lock0:
                fcntl.flock(lock0, fcntl.LOCK_EX)
                truth = generate_file(common, R / 'bias/generator/workspace_postfit.root', 913000, 50, 0., '7.5',
                                      dict(generator_degree=7, fit_degree=6, label='r0', paired_across_masses=True,
                                           generator_fit_sha256=fit_hash))
            truth = dict(truth, mass=float(mass), expected_median=expected,
                         expected_source=str(expected_anchor(mass)), expected_source_sha256=sha(expected_anchor(mass)))
            atomic_json(directory / 'truth.json', truth)
            return truth
        from model import build
        gen_source = directory / 'source.root'
        build(mass, degree=7, family='cheb', output=gen_source, freeze_means=True)
        f = ROOT.TFile.Open(str(gen_source)); w = f.Get('w')
        for q in generator['parameters'] + generator['means']:
            if w.var(q['name']):
                w.var(q['name']).setVal(q['value'])
        set_observable_values(w); w.writeToFile(str(gen_source)); f.Close()
        return generate_file(directory, gen_source, 913001 + round(float(mass)*10)*10, 50, expected, mass,
                             dict(generator_degree=7, fit_degree=6, label='rExpected', paired_across_masses=False,
                                  expected_median=expected, generator_fit_sha256=fit_hash,
                                  expected_source=str(expected_anchor(mass)), expected_source_sha256=sha(expected_anchor(mass))))


def generate_gof(mass='7.5'):
    observed_source = R / f'observed/m{mass}/joint/workspace_postfit_b.root'
    directory = R / 'gof/generation'; directory.mkdir(parents=True, exist_ok=True)
    source = directory / 'source_snapshot.root'
    with (directory / 'generate.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        if not source.exists():
            shutil.copy2(observed_source, source)
            atomic_json(directory / 'source_snapshot.json', dict(source=str(observed_source), source_sha256_at_capture=sha(observed_source), immutable_copy=str(source), sha256=sha(source)))
        prepare_control(source)
        return generate_file(directory, source, 914001, 100, 0., mass,
                             dict(generator_degree=6, fit_degree=6, label='GoF_Bonly', paired_across_masses=True,
                                  meaning='100 common B-only datasets; 6100 r-free saturated scores across61 masses'))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--prepare-only', action='store_true')
    parser.add_argument('--gof', action='store_true')
    parser.add_argument('--mass', nargs='+', default=MASSES)
    parser.add_argument('--label', nargs='+', choices=['r0', 'rExpected'], default=['r0', 'rExpected'])
    args = parser.parse_args()
    if args.prepare_only:
        prepare_generator()
    elif args.gof:
        generate_gof()
    else:
        for mass in args.mass:
            for label in args.label:
                print('GENERATE', mass, label, flush=True); generate_bias(mass, label)
