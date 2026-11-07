"""R-free saturated GoF on common B-only toys, with CR means re-estimated."""
from pathlib import Path
import argparse, fcntl, json, math, re, time, traceback
from types import SimpleNamespace
import ROOT
from toys_helpers import atomic_json, execute, process_cpu_usage
from toys_worker import RUN, C, run_cr, write_joint_workspace, fit_source
from generator_toys import expected_value
from gof_diagnostics import DIAGNOSTIC_COMBINE, RTD_FLAG, read_nominal_diagnostics


def main(args):
    previous_cpu=dict(worker_seconds=0.,native_children_seconds=0.,total_seconds=0.)
    base = RUN / f'gof/toy{args.toy:03d}'
    base.mkdir(parents=True, exist_ok=True)
    truth = json.loads((RUN / 'gof/generation/truth.json').read_text())
    handle = ROOT.TFile.Open(truth['toy_file'])
    toy = handle.Get(f'toys/toy_{args.toy}'); snapshot = handle.Get(f'toys/toy_{args.toy}_snapshot')
    if not toy or not snapshot:
        raise RuntimeError('Missing GoF toy or auxiliary snapshot')
    with (base / 'cr.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        crargs = SimpleNamespace(mass=args.mass[0], fit_timeout=args.fit_timeout)
        cr = run_cr(crargs, base, toy, snapshot)
    for mass in args.mass:
        mass = f'{float(mass):.1f}'
        directory = base / f'm{mass}'; directory.mkdir(exist_ok=True)
        with (directory / 'worker.lock').open('a') as lock:
            fcntl.flock(lock, fcntl.LOCK_EX)
            if (directory / 'result.json').exists():
                continue
            start = time.time()
            result = dict(mass=float(mass), toy=args.toy, score=None, valid=False,
                          CR_valid=cr['valid'], CR_status=cr['selected']['status'],
                          CR_covQual=cr['selected']['covQual'], CR_edm=cr['selected'].get('edm'), CR_edm_valid=cr['selected'].get('edm_valid',False), joint_EDM_available=False, nominal_continuum_positivity_checked=False, CR_boundary_parameters=cr['selected']['boundary_parameters'],
                          CR_mu_boundary=any(n.startswith('Ymean_') for n in cr['selected']['boundary_parameters']), protocol='shared_sr_gof_rfree_v1',
                          paired_across_masses=True, generation_seed=truth['seed'])
            try:
                expected = expected_value(mass); rmax = max(100., 5 * expected)
                # Uniform observed B-only seed (then toy-specific CR seeds/means).
                source = fit_source(mass)
                previous = None
                if source.exists():
                    sf = ROOT.TFile.Open(str(source)); sw = sf.Get('w')
                    previous = dict(parameters=[dict(name=v.GetName(), value=v.getVal()) for v in sw.allVars()
                                               if v.GetName().startswith(('local_', 'CMS_EXO25020_bkgScale_'))], r=0.,
                                    pdf_indices={v.GetName():v.getIndex() for v in sw.allCats() if v.GetName().startswith('pdfindex_')})
                    sf.Close()
                write_joint_workspace(directory, mass, toy, snapshot, cr, 0., rmax, previous)
                # Re-apply toy CR seed values after the observed B-only seed.
                sf = ROOT.TFile.Open(str(directory / 'workspace.root')); sw = sf.Get('w')
                for item in cr['selected']['parameters']:
                    if sw.var(item['name']):
                        sw.var(item['name']).setVal(item['value'])
                for item in cr['means']:
                    sw.var(item['name']).setVal(item['value']); sw.var(item['name']).setConstant(True)
                sw.var('r').setVal(0.)
                sw.writeToFile(str(directory / 'workspace.root')); sf.Close()
                command = [str(DIAGNOSTIC_COMBINE), '-M', 'GoodnessOfFit', str(directory / 'workspace.root'),
                           '-D', 'toy_data', '-m', mass, '-n', '_gof', '--algo', 'saturated',
                           '--trackParameters', 'r', '--rMin', str(-rmax), '--rMax', str(rmax),
                           '--cminDefaultMinimizerStrategy', '0', '--cminDefaultMinimizerTolerance', '.001',
                           '--cminDefaultMinimizerPrecision', '1e-12', '--X-rtd',
                           'MINIMIZER_freezeDisassociatedParams', '--X-rtd', RTD_FLAG, '--cminRunAllDiscreteCombinations', '-v', '0']
                receipt = execute(command, directory, args.fit_timeout)
                files = list(directory.glob('higgsCombine*.root'))
                if receipt['exit_code'] or len(files) != 1:
                    raise RuntimeError('GoF execution failed: ' + str(receipt))
                ff = ROOT.TFile.Open(str(files[0])); tree = ff.Get('limit')
                assert tree and tree.GetEntries() == 1
                tree.GetEntry(0); score = float(tree.limit); rhat = float(tree.trackedParam_r); ff.Close()
                result.update(score=score, r=rhat, r_range=[-rmax,rmax], r_boundary=abs(rhat)>=.99*rmax)
                text = (directory / 'combine.log').read_text(errors='replace')
                flags = [line[:500] for line in text.splitlines() if re.search(r'(Fit failed|Minimization failed|Invalid NLL|\bnan\b)', line, re.I)]
                try:
                    diagnostic=read_nominal_diagnostics(files[0], directory/'workspace.root')
                except Exception as error:
                    diagnostic=dict(available=False,attempted=True,outcome_complete=True,valid=False,flags=['nominal_diagnostic_parse_exception'],exception=repr(error))
                result.update(nominal_diagnostics=diagnostic, joint_EDM_available=diagnostic.get('available',False), nominal_continuum_positivity_checked=diagnostic.get('available',False),
                              nominal_status=diagnostic.get('status'), nominal_edm=diagnostic.get('edm'), nominal_valid=diagnostic.get('valid',False), diagnostic_runtime=str(DIAGNOSTIC_COMBINE))
                result.update(score=score, r=rhat, r_range=[-rmax,rmax],
                              r_boundary=abs(rhat) >= .99*rmax, receipt=receipt,
                              expected_median_anchor=expected, numerical_log_flags=flags,
                              valid=bool(math.isfinite(score) and cr['valid'] and abs(rhat)<.99*rmax and not flags and diagnostic.get('valid',False)))
            except Exception as error:
                result['exception'] = repr(error)
                result.setdefault('nominal_diagnostics',dict(available=False,attempted=False,outcome_complete=True,valid=False,flags=['native_GoF_result_unavailable'],exception=repr(error)))
                (directory / 'exception.txt').write_text(traceback.format_exc())
            result['seconds'] = time.time() - start
            current_cpu=process_cpu_usage()
            result['process_cpu_usage']={k:current_cpu[k]-previous_cpu[k] for k in current_cpu}
            result['process_cpu_usage']['includes_CR_stage']=mass==f'{float(args.mass[0]):.1f}'
            previous_cpu=current_cpu
            atomic_json(directory / 'result.json', result)
            print('GOF', args.toy, mass, result.get('score'), result['valid'], flush=True)
    handle.Close()


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--toy', type=int, required=True)
    p.add_argument('--mass', nargs='+', default=[f'{k/10:.1f}' for k in range(60,121)])
    p.add_argument('--fit-timeout', type=int, default=600)
    args=p.parse_args()
    try:
        main(args)
    except Exception as error:
        # Includes toy-file and CR-stage failures, which happen before the per-mass try.
        for mass in args.mass:
            directory=RUN/f'gof/toy{args.toy:03d}/m{float(mass):.1f}'
            directory.mkdir(parents=True,exist_ok=True)
            with (directory/'worker.lock').open('a') as lock:
                fcntl.flock(lock,fcntl.LOCK_EX)
                if (directory/'result.json').exists():continue
                atomic_json(directory/'result.json',dict(mass=float(mass),toy=args.toy,score=None,valid=False,CR_valid=False,nominal_valid=False,
                    exception=repr(error),failure_flags=['GoF_preparation_or_CR_failure'],
                    nominal_diagnostics=dict(available=False,attempted=False,outcome_complete=True,valid=False,flags=['GoF_preparation_or_CR_failure'],exception=repr(error))))
                (directory/'exception.txt').write_text(traceback.format_exc())
        raise
