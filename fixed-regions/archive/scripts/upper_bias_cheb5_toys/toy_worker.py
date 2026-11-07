"""One full two-stage pseudo-experiment; no ensemble launches on import."""
from pathlib import Path
import argparse, fcntl, hashlib, json, math, os, shutil, sys, time, traceback
import ROOT
from worker_helpers import ERAS, atomic_json, execute, apply_globals, read_fit

ROOT.gROOT.SetBatch(True)
ROOT.gSystem.Load('libHiggsAnalysisCombinedLimit')
P = Path(__file__).resolve().parents[2]
R = P / 'runs/20261030_upper_bias_cheb5_toys'
OLD = P / 'runs/20261030_upper_envelope_scan'
C = str(P / 'bin/combine')
PROTOCOL = 'two_fixed_anchors_v2'


def parser():
    p = argparse.ArgumentParser()
    p.add_argument('--mass', required=True)
    p.add_argument('--injection', choices=['r0', 'rExpected'], required=True)
    p.add_argument('--toy', type=int, required=True)
    p.add_argument('--interval', choices=['profile68', 'hesse68'], default=None,
                   help='Default policy: profile only preselected toy1; other toys retain Hesse diagnostics')
    p.add_argument('--fit-timeout', type=int, default=600)
    p.add_argument('--output-root', type=Path, default=R / 'fits')
    p.add_argument('--generation-root', type=Path, default=R / 'generation')
    p.add_argument('--upgrade-multistart', action='store_true',
                   help='Audit/upgrade a completed legacy result; retain all old artifacts')
    p.add_argument('--recover-control-point', action='store_true',
                   help='Read-only recovery of toy1 point estimates from existing valid anchors; no fits launched')
    return p


def write_cr_workspace(directory, toy, snapshot, reference, variant):
    directory.mkdir(parents=True, exist_ok=True)
    handle = ROOT.TFile.Open(str(OLD / 'degree3/m10.8/control/workspace.root'))
    w = handle.Get('w')
    for a in reference:
        if w.var(a['name']):
            w.var(a['name']).setVal(a['value'])
    if variant == 1:
        for v in w.allVars():
            if v.GetName().startswith('Ymean_'):
                v.setVal((v.getMin() + v.getMax()) / 2 - .07)
            if 'bkgScale_upsilon2s' in v.GetName():
                v.setVal(1.5)
    weight = ROOT.RooRealVar('toy_weight', '', 0, 1e12)
    obs = ROOT.RooArgSet(w.var('mass'), w.cat('CMS_channel'), weight)
    data = ROOT.RooDataSet('toy_data', '', obs, ROOT.RooFit.WeightVar(weight))
    count = 0.
    for k in range(toy.numEntries()):
        row = toy.get(k)
        channel = row.find('CMS_channel').getLabel()
        if not channel.startswith('CR_'):
            continue
        w.cat('CMS_channel').setLabel(channel[3:])
        w.var('mass').setVal(row.getRealValue('mass_CR'))
        weight.setVal(toy.weight())
        data.add(obs, toy.weight())
        count += toy.weight()
    getattr(w, 'import')(data)
    globals_ = apply_globals(w, snapshot)
    w.var('r').setVal(0.)
    w.var('mass').setVal(9.8)
    w.writeToFile(str(directory / 'workspace.root'))
    handle.Close()
    return count, globals_


def common(directory, mass, tolerance):
    return [str(directory / 'workspace.root'), '-D', 'toy_data', '-m', mass,
            '-n', '_upper', '--cminDefaultMinimizerStrategy', '2',
            '--cminDefaultMinimizerTolerance', str(tolerance),
            '--cminDefaultMinimizerPrecision', '1e-12',
            '--X-rtd', 'MINIMIZER_freezeDisassociatedParams', '-v', '0']


def run_cr(args, directory, toy, snapshot):
    final = directory / 'cr_result.json'
    if final.exists():
        return json.loads(final.read_text())
    reference = json.loads((P / 'runs/20261030_upper_sr_cr_mu_fixed/cr_only/results.json').read_text())['selected']['parameters']
    attempts = []
    for variant in [0, 1]:
        wd = directory / f'cr_start{variant}'
        events, globals_ = write_cr_workspace(wd, toy, snapshot, reference, variant)
        command = [C, '-M', 'FitDiagnostics'] + common(wd, args.mass, .001) + ['--skipSBFit', '--minos', 'none', '--robustFit', '0', '--keepFailures']
        receipt = execute(command, wd, args.fit_timeout)
        result = read_fit(wd, False, 'none')
        result.update(variant=variant, receipt=receipt, directory=str(wd))
        attempts.append(result)
    readable = [a for a in attempts if a.get('readable') and a.get('deviance_with_aux') is not None]
    if not readable:
        raise RuntimeError('No readable CR attempt')
    good = [a for a in readable if a['fit_valid'] and a['receipt']['exit_code'] == 0]
    selected = min(good or readable, key=lambda a: a['deviance_with_aux'])
    # One continuation for numerical failures, independent of r/pull values.
    if not good:
        wd = directory / 'cr_retry'
        write_cr_workspace(wd, toy, snapshot, selected['parameters'], 0)
        command = [C, '-M', 'FitDiagnostics'] + common(wd, args.mass, .00001) + ['--skipSBFit', '--minos', 'none', '--robustFit', '0', '--keepFailures']
        receipt = execute(command, wd, args.fit_timeout)
        result = read_fit(wd, False, 'none')
        result.update(variant=2, receipt=receipt, directory=str(wd))
        attempts.append(result)
        if result.get('fit_valid') and receipt['exit_code'] == 0:
            selected = result
    means = [a for a in selected['parameters'] if a['name'].startswith('Ymean_')]
    if len(means) != 8:
        raise RuntimeError('Expected eight fitted resonance means')
    result = dict(selected=selected, means=means, attempts=attempts, events=events,
                  global_observables=globals_, valid=selected['fit_valid'] and selected['receipt']['exit_code'] == 0)
    atomic_json(final, result)
    return result


def write_joint_workspace(directory, mass, toy, snapshot, cr, rtrue, rmax, previous=None):
    directory.mkdir(parents=True, exist_ok=True)
    source = OLD / f'degree3/m{mass}/joint/workspace_seeded.root'
    handle = ROOT.TFile.Open(str(source))
    w = handle.Get('w')
    getattr(w, 'import')(toy.Clone('toy_data'))
    global_values = apply_globals(w, snapshot)
    for a in cr['selected']['parameters']:
        name = a['name'] + '_CR' if a['name'].startswith('local_') else a['name']
        if w.var(name) and not name.startswith('Ymean_'):
            w.var(name).setVal(a['value'])
    for a in cr['means']:
        w.var(a['name']).setVal(a['value'])
        w.var(a['name']).setConstant(True)
    if previous:
        for a in previous['parameters']:
            if w.var(a['name']) and a['name'] != 'r':
                w.var(a['name']).setVal(a['value'])
        for name, value in previous['pdf_indices'].items():
            w.cat(name).setIndex(value)
    w.var('r').setRange(-rmax, rmax)
    w.var('r').setVal(previous['r'] if previous else rtrue)
    w.var('mass').setVal(11.2)
    w.var('mass_CR').setVal(9.8)
    means_audit = {a['name']: dict(value=w.var(a['name']).getVal(),
                                 CR_value=a['value'], constant=bool(w.var(a['name']).isConstant()))
                   for a in cr['means']}
    assert len(means_audit) == 8 and all(a['constant'] and a['value'] == a['CR_value'] for a in means_audit.values())
    atomic_json(directory / 'workspace_audit.json', dict(means=means_audit,
                global_observables=global_values, toy_events=toy.sumEntries(),
                CR_signal_included=True, discrete_indices=4,
                source=str(source), source_sha256=hashlib.sha256(source.read_bytes()).hexdigest()))
    w.writeToFile(str(directory / 'workspace.root'))
    handle.Close()
    return str(source)


def point_valid(result):
    value = result.get('deviance_with_aux')
    return bool(result.get('readable') and result.get('fit_valid') and
                result.get('receipt', {}).get('exit_code') == 0 and
                value is not None and math.isfinite(value))


def run_anchors(args, directory, toy, snapshot, cr, expected, rmax):
    anchors = []
    for label, anchor in [('zero', 0.), ('expected', expected)]:
        attempts = []
        previous = None
        for attempt in range(2):
            wd = directory / (f'anchor_v2_{label}' if not attempt else f'anchor_v2_{label}_retry')
            # Each primary is independent, sharing only the fixed observed model/CR seeds.
            write_joint_workspace(wd, args.mass, toy, snapshot, cr, anchor, rmax, previous)
            command = [C, '-M', 'FitDiagnostics'] + common(wd, args.mass, .001 if not attempt else .00001)
            command += ['--rMin', str(-rmax), '--rMax', str(rmax),
                        '--cminRunAllDiscreteCombinations', '--skipBOnlyFit',
                        '--keepFailures', '--preFitValue', str(previous['r'] if previous else anchor),
                        '--minos', 'none', '--robustFit', '0']
            receipt = execute(command, wd, args.fit_timeout)
            result = read_fit(wd, True, 'hesse68')
            result.update(attempt=attempt, directory=str(wd), receipt=receipt,
                          anchor_label=label, anchor_value=anchor)
            attempts.append(result)
            if point_valid(result):
                break
            if result.get('readable'):
                previous = result
        readable = [q for q in attempts if q.get('readable') and q.get('deviance_with_aux') is not None and math.isfinite(q['deviance_with_aux'])]
        valid = [q for q in readable if point_valid(q)]
        selected = min(valid or readable, key=lambda q: q['deviance_with_aux']) if readable else attempts[-1]
        anchor_result = dict(label=label, initial_r=anchor, attempts=attempts,
                             selected=selected, valid=point_valid(selected))
        anchors.append(anchor_result)
    readable = [a['selected'] for a in anchors if a['selected'].get('deviance_with_aux') is not None and math.isfinite(a['selected']['deviance_with_aux'])]
    valid = [a['selected'] for a in anchors if a['valid']]
    if not readable:
        atomic_json(directory / 'multistart_v2.json', dict(protocol=PROTOCOL,
                    anchors=anchors, valid=False, error='Neither anchor produced a readable minimum'))
        raise RuntimeError('Neither fixed anchor produced a readable simultaneous minimum')
    selected = min(valid or readable, key=lambda a: a['deviance_with_aux'])
    meta = dict(protocol=PROTOCOL, anchors=anchors, selected_anchor=selected['anchor_label'],
                selected_deviance=selected['deviance_with_aux'],
                valid=bool(valid), both_anchors_valid=all(a['valid'] for a in anchors),
                expected_median_anchor=expected, r_range=[-rmax, rmax],
                anchor_policy='Same fixed r=0 and r=expectedMedian starts for BOTH injection ensembles')
    atomic_json(directory / 'multistart_v2.json', meta)
    return selected, meta


def legacy_profile_candidate(old):
    if not old or not old.get('directory'):
        return None
    try:
        selected = read_fit(Path(old['directory']), True, old.get('interval_method', 'profile68'))
        selected.update(directory=old['directory'], receipt=old.get('receipt', {}),
                        attempt=old.get('selected_attempt', old.get('attempt', 0)))
        return selected
    except Exception:
        return None


def retain_point_without_profile(args, directory, best_anchor, multistart, attempts, reason):
    """A failed uncertainty calculation cannot erase a converged point estimate."""
    candidates = [best_anchor] if point_valid(best_anchor) else []
    if not candidates:
        raise RuntimeError('No valid point estimate available after profile failure')
    selected = min(candidates, key=lambda q: q['deviance_with_aux'])
    result = dict(selected)
    result.update(source=str(OLD / f'degree3/m{args.mass}/joint/workspace_seeded.root'),
                  attempts=attempts, selected_attempt=selected.get('attempt', 0),
                  profile_reused=False, profile_computed=False, profile_attempted=True,
                  profile_failed_point_retained=True, profile_failure_reason=reason,
                  point_from_anchor=selected.get('anchor_label') is not None,
                  multistart=multistart, multistart_valid=multistart['valid'],
                  r_range=multistart['r_range'],
                  profile_anchor_deviance_difference=selected['deviance_with_aux'] - best_anchor['deviance_with_aux'],
                  interval_method='profile68', interval_valid=False,
                  interval_low=None, interval_high=None, error_lo=None, error_hi=None,
                  interval_boundary=False)
    return result


def recover_control_point(args, directory):
    """Use already-saved fits only: this function never calls Combine."""
    multistart = json.loads((directory / 'multistart_v2.json').read_text())
    anchors = [a['selected'] for a in multistart['anchors'] if point_valid(a['selected'])]
    if not anchors:
        raise RuntimeError('Control recovery requires an existing valid anchor')
    best_anchor = min(anchors, key=lambda q: q['deviance_with_aux'])
    attempts = []
    for i, name in enumerate(['profile_v2', 'profile_v2_retry']):
        wd = directory / name
        if not (wd / 'receipt.json').exists():
            continue
        receipt = json.loads((wd / 'receipt.json').read_text())
        try:
            result = read_fit(wd, True, 'profile68')
        except Exception as error:
            result = dict(readable=False, failure=repr(error))
        result.update(attempt=i, directory=str(wd), receipt=receipt)
        attempts.append(result)
    return retain_point_without_profile(args, directory, best_anchor, multistart, attempts,
                'Recovery from stored anchors after unsuccessful control profile; no new fits')


def run_joint(args, directory, toy, snapshot, cr, expected, old=None):
    if args.recover_control_point:
        return recover_control_point(args, directory)
    # The numerical range, like both seeds, is identical for both injection ensembles.
    rmax = max(100., 5 * expected)
    best_anchor, multistart = run_anchors(args, directory, toy, snapshot, cr, expected, rmax)
    legacy = legacy_profile_candidate(old)
    reuse_point = (legacy and point_valid(legacy) and
                   legacy['deviance_with_aux'] <= best_anchor['deviance_with_aux'] + 1e-4)
    reuse_interval = (reuse_point and legacy.get('interval_valid') and
                      legacy.get('interval_method') == args.interval and
                      -rmax < legacy['interval_low'] < legacy['interval_high'] < rmax)
    if reuse_point and (args.interval == 'hesse68' or reuse_interval):
        result = dict(legacy)
        result.update(source=old.get('source'), attempts=old.get('attempts', [legacy]),
                      selected_attempt=legacy['attempt'],
                      profile_reused=legacy.get('interval_method') == 'profile68',
                      legacy_point_reused=True,
                      multistart=multistart, multistart_valid=multistart['valid'],
                      profile_anchor_deviance_difference=legacy['deviance_with_aux'] - best_anchor['deviance_with_aux'],
                      r_range=[-rmax, rmax])
        return result
    if args.interval == 'hesse68':
        # The user-selected fast study needs unbiased point estimates for all toys.
        # Both full envelope anchors already provide that fit and Hesse diagnostics.
        # Do not repeat an identical fit or invent profile-based coverage.
        result = dict(best_anchor)
        result.update(source=str(OLD / f'degree3/m{args.mass}/joint/workspace_seeded.root'),
                      attempts=[best_anchor], selected_attempt=best_anchor['attempt'],
                      profile_reused=False, legacy_point_reused=False,
                      point_from_anchor=True, profile_computed=False,
                      multistart=multistart, multistart_valid=multistart['valid'],
                      profile_anchor_deviance_difference=0.,
                      minimum_consistent_with_anchors=True, r_range=[-rmax, rmax])
        return result
    attempts = []
    previous = best_anchor
    if legacy and point_valid(legacy) and legacy['deviance_with_aux'] < best_anchor['deviance_with_aux']:
        previous = legacy
    for attempt in range(2):
        wd = directory / ('profile_v2' if attempt == 0 else 'profile_v2_retry')
        source = write_joint_workspace(wd, args.mass, toy, snapshot, cr, previous['r'], rmax, previous)
        command = [C, '-M', 'FitDiagnostics'] + common(wd, args.mass, .001 if attempt == 0 else .00001)
        command += ['--rMin', str(-rmax), '--rMax', str(rmax), '--cminRunAllDiscreteCombinations',
                    '--skipBOnlyFit', '--keepFailures', '--preFitValue', str(previous['r'])]
        if args.interval == 'profile68':
            command += ['--minos', 'poi', '--robustFit', '1', '--stepSize', '.01',
                        '--setRobustFitStrategy', '1', '--setRobustFitTolerance', '.001']
        else:
            command += ['--minos', 'none', '--robustFit', '0']
        receipt = execute(command, wd, args.fit_timeout)
        result = read_fit(wd, True, args.interval)
        result.update(attempt=attempt, directory=str(wd), receipt=receipt)
        attempts.append(result)
        if result.get('readable') and result.get('deviance_with_aux') is not None:
            previous = result
        hesse_ok = result.get('covQual') == 3 and result.get('error', 0.) > 0
        minimum_consistent = (result.get('deviance_with_aux') is not None and
                              result['deviance_with_aux'] <= best_anchor['deviance_with_aux'] + 1e-4)
        result['minimum_consistent_with_anchors'] = minimum_consistent
        good = point_valid(result) and result.get('interval_valid') and minimum_consistent
        if good and (args.interval == 'profile68' or hesse_ok):
            break
    readable = [a for a in attempts if a.get('readable') and a.get('deviance_with_aux') is not None]
    if not readable:
        return retain_point_without_profile(args, directory, best_anchor, multistart, attempts,
                    'Profile attempts unreadable or timed out')
    good = [a for a in readable if point_valid(a) and a.get('interval_valid') and a.get('minimum_consistent_with_anchors')]
    if not good:
        return retain_point_without_profile(args, directory, best_anchor, multistart, attempts,
                    'No valid profile interval at a minimum consistent with both fixed anchors')
    selected = min(good or readable, key=lambda a: a['deviance_with_aux'])
    result = dict(selected)
    result.update(source=source, attempts=attempts, selected_attempt=selected['attempt'],
                  profile_reused=False, multistart=multistart,
                  multistart_valid=multistart['valid'], r_range=[-rmax, rmax],
                  profile_anchor_deviance_difference=selected['deviance_with_aux'] - best_anchor['deviance_with_aux'])
    if not selected.get('minimum_consistent_with_anchors'):
        result['fit_valid'] = False
        result['interval_valid'] = False
    return result


def main(args):
    if not 1 <= args.toy <= 100:
        raise ValueError('Toy must be 1..100')
    interval_requested = args.interval
    args.interval = args.interval or ('profile68' if args.toy == 1 else 'hesse68')
    args.mass = f'{float(args.mass):.1f}'
    args.output_root = args.output_root.resolve()
    args.generation_root = args.generation_root.resolve()
    directory = args.output_root / f'm{args.mass}/{args.injection}/toy{args.toy:03d}'
    directory.mkdir(parents=True, exist_ok=True)
    lock = (directory / 'worker.lock').open('a')
    fcntl.flock(lock, fcntl.LOCK_EX)
    result_path = directory / 'result.json'
    old = json.loads(result_path.read_text()) if result_path.exists() else None
    if args.recover_control_point:
        if args.toy != 1 or not old:
            raise ValueError('Control-point recovery is restricted to completed toy1 results')
        if old.get('fit_valid'):
            print('POINT ALREADY USABLE', result_path, flush=True)
            return
        saved = json.loads((directory / 'multistart_v2.json').read_text())
        if not any(point_valid(a['selected']) for a in saved.get('anchors', [])):
            raise ValueError('No saved valid anchor: refusing read-only point recovery')
        if not (directory / 'cr_result.json').exists():
            raise ValueError('Missing saved CR fit: recovery must not launch fits')
    elif old and (old.get('protocol') == PROTOCOL or not args.upgrade_multistart):
        print('EXISTS', result_path, old.get('protocol', 'legacy_v1'), flush=True)
        return
    if old:
        archive = directory / ('before_control_point_recovery' if args.recover_control_point else 'legacy_v1')
        archive.mkdir(exist_ok=True)
        if not (archive / 'result.json').exists():
            shutil.copy2(result_path, archive / 'result.json')
    start = time.time()
    result = dict(mass=float(args.mass), label=args.injection, injection=args.injection,
                  toy=args.toy, generator_degree=5, fit_degree=3,
                  fit_model='Chebyshev3/Bernstein3 discrete envelope', interval_method=args.interval,
                  valid=False, fit_valid=False, interval_valid=False, failure_flags=[],
                  protocol=PROTOCOL, upgraded_from_legacy=bool(old),
                  recovered_control_point=bool(args.recover_control_point),
                  interval_policy='profile_controls_toy1_else_hesse_diagnostic',
                  maininterval_policy='profile_controls_toy1_else_hesse_diagnostic',
                  interval_requested=interval_requested or 'automatic_by_preselected_toy_id',
                  preselected_profile_control=args.toy == 1,
                  ensemble_coverage_evaluated=False,
                  worker_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    try:
        generation = args.generation_root / f'm{args.mass}/{args.injection}'
        truth = json.loads((generation / 'truth.json').read_text())
        result.update(rtrue=float(truth['rtrue']), generation_seed=truth['seed'],
                      truth_source=str(generation / 'truth.json'))
        files = list(generation.glob('higgsCombine*.root'))
        if len(files) != 1:
            raise ValueError(f'Expected one generated toy file, got {files}')
        handle = ROOT.TFile.Open(str(files[0]))
        toy = handle.Get(f'toys/toy_{args.toy}')
        snapshot = handle.Get(f'toys/toy_{args.toy}_snapshot')
        if not toy or not snapshot:
            raise ValueError('Missing toy or auxiliary snapshot')
        result['toy_events'] = toy.sumEntries()
        cr = run_cr(args, directory, toy, snapshot)
        expected = float(truth['expected_median'])
        result['expected_median_anchor'] = expected
        joint = run_joint(args, directory, toy, snapshot, cr, expected, old)
        result.update(joint)
        result.update(CR_valid=cr['valid'], CR_status=cr['selected']['status'],
                      CR_covQual=cr['selected']['covQual'], CR_boundary_parameters=cr['selected']['boundary_parameters'])
        result['joint_fit_valid'] = result['fit_valid']
        result['joint_interval_valid'] = result['interval_valid']
        result['fit_valid'] = bool(result['joint_fit_valid'] and cr['valid'] and result['receipt']['exit_code'] == 0 and result['multistart_valid'])
        result['interval_valid'] = bool(result['joint_interval_valid'] and result['fit_valid'])
        result['bias'] = result['r'] - result['rtrue']
        result['bias_sigma'] = result['bias'] / result['error'] if result['error'] > 0 else None
        result['hesse_valid'] = result['covQual'] == 3 and result['error'] > 0 and math.isfinite(result['error'])
        if result['interval_method'] == 'hesse68':
            result['interval_valid'] = bool(result['interval_valid'] and result['hesse_valid'])
        toward_truth = result['error_lo'] if result['bias'] >= 0 else result['error_hi']
        result['pull'] = (result['bias'] / toward_truth if result['interval_valid'] and
                          toward_truth is not None and math.isfinite(toward_truth) and toward_truth > 0 else None)
        result['covered68'] = (result['interval_low'] <= result['rtrue'] <= result['interval_high']
                               if result['interval_valid'] else None)
        result['valid'] = bool(result['fit_valid'] and result['interval_valid'] and cr['valid'] and result['receipt']['exit_code'] == 0 and (result['interval_method'] == 'profile68' or result['hesse_valid']))
        if not cr['valid']:
            result['failure_flags'].append('CR_fit_invalid')
        if result.get('profile_failed_point_retained'):
            result['failure_flags'].append('profile_failed_point_retained')
        if not result['multistart_valid']:
            result['failure_flags'].append('multistart_no_converged_anchor')
        elif not result['multistart']['both_anchors_valid']:
            result['failure_flags'].append('one_anchor_numerically_invalid')
        if result.get('profile_anchor_deviance_difference', 0) > 1e-4:
            result['failure_flags'].append('profile_minimum_worse_than_anchor')
        if not result['fit_valid']:
            result['failure_flags'].append('simultaneous_fit_invalid')
        if not result['interval_valid']:
            result['failure_flags'].append('interval_invalid')
        if not result['hesse_valid']:
            result['failure_flags'].append('Hesse_covariance_not_full_quality')
        if result['r_boundary']:
            result['failure_flags'].append('r_at_boundary')
        if result['interval_boundary']:
            result['failure_flags'].append('interval_at_boundary')
        if result['receipt']['exit_code']:
            result['failure_flags'].append('combine_nonzero_exit')
        handle.Close()
    except Exception as error:
        result['exception'] = repr(error)
        result['failure_flags'].append('worker_exception')
        (directory / 'exception.txt').write_text(traceback.format_exc())
    result['seconds'] = time.time() - start
    atomic_json(result_path, result)
    print(json.dumps({k: result.get(k) for k in ['mass', 'label', 'toy', 'rtrue', 'r', 'error', 'error_lo', 'error_hi', 'pull', 'covered68', 'valid', 'failure_flags', 'seconds']}), flush=True)


if __name__ == '__main__':
    main(parser().parse_args())
