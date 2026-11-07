"""Read preserved GoF inputs; compare isolated diagnostic scores to original scores."""
from pathlib import Path
import hashlib
import json
import subprocess
import time
import ROOT
from toys_helpers import continuum_positivity, exact_deviance

PROJECT = Path(__file__).resolve().parents[2]
RUN = PROJECT / 'runs/20261031_shared_sr_cr_scan'
RUNTIME = PROJECT / 'work/shared_gof_diagnostic_runtime/v2'


def row(directory):
    files = list(directory.glob('higgsCombine*.root'))
    assert len(files) == 1, files
    handle = ROOT.TFile.Open(str(files[0]))
    tree = handle.Get('limit'); assert tree and tree.GetEntries() == 1
    tree.GetEntry(0)
    result = {b.GetName(): float(getattr(tree, b.GetName())) for b in tree.GetListOfBranches()}
    leaf_types = {b.GetName(): b.GetLeaf(b.GetName()).GetTypeName() for b in tree.GetListOfBranches()}
    handle.Close()
    return result, leaf_types


def main():
    cases = [
        ('toy007_6.0', RUN / 'gof/toy007/m6.0', False),
    ]
    results = []
    for name, original, run_off in cases:
        original_values, _ = row(original)
        command_file = original / 'combine.command.json'
        if not command_file.exists():
            command_file = original / 'command.json'
        original_command = json.loads(command_file.read_text())
        source = next(Path(x) for x in original_command if x.endswith('.root'))
        source_hash = hashlib.sha256(source.read_bytes()).hexdigest()
        case = dict(case=name, original_directory=str(original), original_values=original_values,
                    input=str(source), input_sha256=source_hash, attempts={})
        for enabled in ([False, True] if run_off else [True]):
            variant = 'on' if enabled else 'off'
            directory = RUNTIME / 'validation' / name / variant
            directory.mkdir(parents=True, exist_ok=True)
            if list(directory.glob('higgsCombine*.root')):
                raise RuntimeError('Refusing to overwrite completed parity job ' + str(directory))
            command = [str(RUNTIME / 'combine')] + original_command[1:]
            if enabled:
                command.extend(['--X-rtd', 'GOF_TRACK_NOMINAL'])
            (directory / 'command.json').write_text(json.dumps(command, indent=2))
            start = time.time()
            with (directory / 'combine.log').open('w') as log:
                result = subprocess.run(command, cwd=directory, stdout=log, stderr=subprocess.STDOUT, timeout=600)
            assert result.returncode == 0, (name, result.returncode)
            values, types = row(directory)
            entry = dict(seconds=time.time()-start, values=values,
                         score_difference_from_original=values['limit']-original_values['limit'],
                         native_tracked_r_difference=values['trackedParam_r']-original_values['trackedParam_r'])
            assert abs(entry['score_difference_from_original']) < 1e-8, entry
            assert abs(entry['native_tracked_r_difference']) < 1e-8, entry
            if enabled:
                diagnostics = {key.removeprefix('gof_nominal_'): value for key,value in values.items() if key.startswith('gof_nominal_')}
                assert all(types['gof_nominal_'+key] == 'Double_t' for key in diagnostics)
                assert sum(key.startswith('pdfindex_') for key in diagnostics) == 4
                assert diagnostics['fitresult_available'] == 1
                assert diagnostics['fitresult_state_max_abs_difference'] < 1e-12
                handle = ROOT.TFile.Open(str(source)); workspace = handle.Get('w')
                for key,value in diagnostics.items():
                    if workspace.var(key):
                        workspace.var(key).setVal(value)
                    elif workspace.cat(key):
                        workspace.cat(key).setIndex(int(value))
                mc = workspace.obj('ModelConfig')
                parameters = mc.GetPdf().getParameters(mc.GetObservables())
                floating = [p.GetName() for p in parameters if p.InheritsFrom('RooRealVar') and not p.isConstant()]
                entry['uncaptured_floating_parameters'] = [name for name in floating if name not in diagnostics]
                assert not entry['uncaptured_floating_parameters'], entry['uncaptured_floating_parameters']
                entry['all_modelconfig_nuisances_captured'] = all(p.GetName() in diagnostics for p in mc.GetNuisanceParameters())
                assert entry['all_modelconfig_nuisances_captured']
                assert 'lumi_1' in diagnostics and 'lumi_2' in diagnostics
                deviance, min_bin = exact_deviance(workspace)
                multipdfs = [p for p in mc.GetPdf().getComponents() if p.InheritsFrom('RooMultiPdf')]
                assert len(multipdfs) == 8
                penalty = 2 * sum(p.getCorrection() for p in multipdfs)
                entry.update(exact_deviance=deviance, minimum_bin_expectation=min_bin,
                             envelope_penalty=penalty,
                             score_minus_reconstructed=values['limit']-deviance-penalty)
                assert abs(entry['score_minus_reconstructed']) < 1e-6, entry['score_minus_reconstructed']
                entry['exact_continuum_minima'] = continuum_positivity(workspace)
                entry['diagnostic_count'] = len(diagnostics)
                entry['selected_indices'] = {k:v for k,v in diagnostics.items() if k.startswith('pdfindex_')}
                handle.Close()
            case['attempts'][variant] = entry
            print(name, variant, 'D', values['limit'], 'delta', entry['score_difference_from_original'], flush=True)
        assert hashlib.sha256(source.read_bytes()).hexdigest() == source_hash
        results.append(case)
        (RUNTIME / 'validation/parity.json').write_text(json.dumps(dict(all_score_parities_passed=True, cases=results),indent=2)+'\n')
    print('ALL PARITY CHECKS PASSED', flush=True)


if __name__ == '__main__':
    main()
