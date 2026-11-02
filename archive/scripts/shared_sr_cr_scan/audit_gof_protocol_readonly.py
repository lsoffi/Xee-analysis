"""Independent read-only protocol audit; never launches a fit or modifies inputs."""
from pathlib import Path
import collections
import datetime
import hashlib
import json
import math
import struct
import numpy as np
import ROOT

P = Path(__file__).resolve().parents[2]
R = P / 'runs/20261031_shared_sr_cr_scan'
S = Path(__file__).resolve().parent
ROOT.gROOT.SetBatch(True)
ROOT.gSystem.Load('libHiggsAnalysisCombinedLimit')


def load(path):
    return json.loads(Path(path).read_text())


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def data_digest(data, cr_only=False):
    digest = hashlib.sha256(); entries = 0; events = 0.
    for k in range(data.numEntries()):
        row = data.get(k); channel = row.find('CMS_channel').getLabel()
        if cr_only and not channel.startswith('CR_'):
            continue
        xname = 'mass_' + channel.split('_')[0]
        digest.update(channel.encode())
        digest.update(struct.pack('dd', row.getRealValue(xname), data.weight()))
        entries += 1; events += data.weight()
    return dict(sha256=digest.hexdigest(), entries=entries, events=events)


def main():
    truth = load(R / 'gof/generation/truth.json')
    generated = ROOT.TFile.Open(truth['toy_file'])
    source_file = ROOT.TFile.Open(truth['source']); source = source_file.Get('w')
    generator_file = ROOT.TFile.Open(str(R / 'gof/generation/generator.root')); generator = generator_file.Get('w')
    mc = generator.obj('ModelConfig')
    globals_names = sorted(x.GetName() for x in mc.GetGlobalObservables())
    assert len(globals_names) == 10
    assert truth['ntoys'] == 100 and truth['rtrue'] == 0
    assert generator.var('r').getVal() == 0
    assert sha(truth['toy_file']) == truth['toy_file_sha256']
    assert sha(truth['source']) == truth['source_sha256']
    assert sha(R / 'gof/generation/generator.root') == truth['generator_workspace_sha256']
    nuisance_names = [n[:-3] for n in globals_names]
    truth_nuisances = {name:generator.var(name).getVal() for name in nuisance_names}
    assert all(not generator.var(name).isConstant() for name in nuisance_names)
    assert all(generator.var(name).getVal() == source.var(name).getVal() for name in nuisance_names)
    true_means = {v.GetName():v.getVal() for v in generator.allVars() if v.GetName().startswith('Ymean_')}
    indices = {c.GetName():c.getIndex() for c in generator.allCats() if c.GetName().startswith('pdfindex_')}
    assert len(indices) == 4 and all(generator.cat(name).isConstant() for name in indices)
    assert generator.genobj('discreteParams').getSize() == 0
    observed_cr = load(P / 'runs/20261030_upper_sr_cr_mu_fixed/cr_only/results.json')
    assert all(true_means[p['name']] == p['value'] for p in observed_cr['means'])
    draws = []; toy_data = {}; snapshots = {}
    for toy_id in range(1,101):
        toy = generated.Get(f'toys/toy_{toy_id}')
        snapshot = generated.Get(f'toys/toy_{toy_id}_snapshot')
        assert toy and snapshot and snapshot.getSize() == 10
        toy_data[toy_id] = data_digest(toy)
        assert toy_data[toy_id]['entries'] == 1800
        snapshots[toy_id] = {name:snapshot.find(name).getVal() for name in globals_names}
        draws.append([snapshots[toy_id][name] for name in globals_names])
    assert len({q['sha256'] for q in toy_data.values()}) == 100
    assert len({tuple(snapshots[i].values()) for i in snapshots}) == 100
    draw_array = np.asarray(draws)
    globals_summary = {name:dict(sample_mean=float(draw_array[:,i].mean()),sample_sd=float(draw_array[:,i].std(ddof=1)),truth=truth_nuisances[name[:-3]]) for i,name in enumerate(globals_names)}
    samples = []
    for toy_id in [1,7,15]:
        base = R / f'gof/toy{toy_id:03d}'
        if not (base / 'cr_result.json').exists():
            continue
        cr = load(base / 'cr_result.json')
        selected_means = {p['name']:p['value'] for p in cr['means']}
        cr_file = ROOT.TFile.Open(str(Path(cr['selected']['directory']) / 'workspace.root')); crw = cr_file.Get('w')
        cr_digest = data_digest(crw.data('toy_data'),True)
        assert cr_digest == data_digest(generated.Get(f'toys/toy_{toy_id}'),True)
        assert cr_digest['entries'] == 424
        assert all(not crw.var(name).isConstant() for name in selected_means)
        cr_globals = {g.GetName():g.getVal() for g in crw.obj('ModelConfig').GetGlobalObservables()}
        assert cr_globals == snapshots[toy_id]
        cr_file.Close()
        for mass in ['6.0','7.5','9.5','11.5','12.0']:
            directory = base / f'm{mass}'
            path = directory / 'workspace.root'
            if not (directory/'result.json').exists() or not path.exists():
                continue
            q = load(directory/'result.json')
            before = sha(path)
            f = ROOT.TFile.Open(str(path)); w=f.Get('w')
            actual_digest = data_digest(w.data('toy_data'))
            assert actual_digest == toy_data[toy_id]
            joint_globals = {g.GetName():g.getVal() for g in w.obj('ModelConfig').GetGlobalObservables()}
            assert joint_globals == cr_globals
            assert all(w.var(name).isConstant() and w.var(name).getVal() == value for name,value in selected_means.items())
            assert not w.var('r').isConstant()
            cats = [c for c in w.allCats() if c.GetName().startswith('pdfindex_')]
            assert len(cats) == 4 and w.genobj('discreteParams').getSize() == 4
            corrections = [p.getCorrection() for p in w.obj('ModelConfig').GetPdf().getComponents() if p.InheritsFrom('RooMultiPdf')]
            assert len(corrections) == 8 and abs(2*sum(corrections)-28) < 1e-12
            diag = q.get('nominal_diagnostics',{})
            samples.append(dict(toy=toy_id,mass=float(mass),data=actual_digest,CR_data=cr_digest,
                globals_equal_snapshot_and_CR=True,means_reestimated_and_fixed=True,
                max_abs_mean_shift_from_generator=max(abs(value-true_means[name]) for name,value in selected_means.items()),
                r_floating=True,r_range=[w.var('r').getMin(),w.var('r').getMax()],
                CR_status=cr['selected']['status'],CR_EDM=cr['selected']['edm'],CR_valid=cr['valid'],
                discrete_indices=4,twice_NLL_discrete_penalty=28.,
                nominal_status=diag.get('status'),nominal_EDM=diag.get('edm'),nominal_QC_valid=diag.get('valid'),
                nominal_score_reconstruction_difference=diag.get('score_state_difference')))
            f.Close();assert sha(path) == before
    observed=[]
    for mass in ['6.0','7.5','9.5','11.5','12.0']:
        directory=R/f'observed/m{mass}'
        command=directory/'joint/gof_observed/combine.command.json'
        if not command.exists():continue
        cmd=load(command);workspace=Path(cmd[cmd.index('-d')+1])
        f=ROOT.TFile.Open(str(workspace));w=f.Get('w')
        assert all(w.var(name).isConstant() and w.var(name).getVal()==value for name,value in true_means.items())
        assert not w.var('r').isConstant() and '--fixedSignalStrength' not in cmd
        obs_globals={g.GetName():g.getVal() for g in w.obj('ModelConfig').GetGlobalObservables()}
        assert all(value==0 for value in obs_globals.values())
        observed.append(dict(mass=float(mass),CR_means_equal_reference=True,r_floating=True,global_observables=obs_globals,
                             r_range=[float(cmd[cmd.index('--rMin')+1]),float(cmd[cmd.index('--rMax')+1])]))
        f.Close()
    results=[load(path) for path in (R/'gof').glob('toy*/m*/result.json')]
    unique_tasks={(float(q['mass']),q['toy']) for q in results}
    assert len(unique_tasks)==len(results)
    all_masses={k/10 for k in range(60,121)}
    assert all(mass in all_masses and 1<=toy<=100 for mass,toy in unique_tasks)
    status=dict(completed_outcomes=len(results),target=6100,unique_task_keys=len(unique_tasks),
                finite_scores=sum(q.get('score') is not None and math.isfinite(q['score']) for q in results),
                CR_invalid=sum(q.get('CR_valid') is not True for q in results),
                r_boundaries=sum(q.get('r_boundary',False) for q in results),
                nominal_diagnostic_flags=dict(collections.Counter(flag for q in results for flag in q.get('nominal_diagnostics',{}).get('flags',[]))))
    sources=['generator_toys.py','toys_worker.py','toys_gof_worker.py','toys_helpers.py','gof_diagnostics.py','toys_qc_outcomes.py','toys_run.py','report.py','audit_observed_gof.py','toys_profiles_run.py']
    output=dict(audit_date=datetime.datetime.now().astimezone().isoformat(),read_only=True,new_fits=0,
                scope='GoF statistical protocol and observed/toy comparability; model geometry/integrals audited separately',
                passed_structural_checks=True,generation=dict(r=0,source=truth['source'],source_sha256=truth['source_sha256'],
                generator_family_indices=indices,unique_datasets=100,unique_auxiliary_snapshots=100,
                paired_across_61_masses=True,truth_nuisances=truth_nuisances,global_observable_sample_summary=globals_summary),
                observed_CR=dict(status=observed_cr['selected']['status'],EDM=observed_cr['selected']['edm'],
                attempts=len(observed_cr['attempts']),deviances=[q['deviance'] for q in observed_cr['attempts']]),
                observed_samples=observed,toy_samples=samples,current_status=status,
                required_model_or_generation_corrections=[],
                interpretation_limits=[
                    'Conditional parametric bootstrap at the fitted B-only nuisance values and selected generating family; not an exact uniformly calibrated composite-null p-value.',
                    'CR-only mean estimation is repeated on each toy, then means are fixed in the joint fit. This reproduces the plug-in two-stage procedure but does not profile the mean uncertainty in observed limits or significance.',
                    'The statistic profiles r, but the reference ensemble is generated at r=0 only. This is the requested background-reference GoF, not a calibration over arbitrary nonzero signal truths.',
                    'A shared 100-toy ensemble gives 6100 correlated mass-dependent scores, not 6100 independent toys. No look-elsewhere calibration is produced.',
                    'Raw p uses all finite scores. Numerical/positivity failures and missing scores must remain visible; a converged subset is not an unconditional ensemble.',
                    'Observed CR reference used three starts; toy CR uses two fixed starts plus one retry on numerical failure. The mathematical estimator/model is the same, but equality of numerical global minimization is not guaranteed by status/EDM alone.',
                    'Observed joint GoF starts from an observed B-only joint fit; toy joint GoF starts from observed SR seeds plus toy CR seeds. This affects numerical exploration, not the likelihood; retain numerical flags and the known-better-FitDiagnostics check.',
                    '100 toys imply coarse Monte Carlo precision; zero exceedances are not proof that the underlying p-value is zero.'
                ],
                reporting_followups=['Explicitly label GoF as a conditional parametric bootstrap; avoid calling it an exact composite-null calibration.',
                    'Keep failed observed nominal fits separate from valid ones; diagnostic completion alone is not scientific validity.',
                    'For any future r-boundary outcome, also exclude it from resolved numerical-QC counts or clearly define the GoF as a bounded-r test. None is present in this audit snapshot.'],
                source_sha256={name:sha(S/name) for name in sources})
    for f in [generated,source_file,generator_file]:f.Close()
    (R/'gof_protocol_independent_audit.json').write_text(json.dumps(output,indent=2,allow_nan=False)+'\n')
    print(json.dumps(dict(passed=output['passed_structural_checks'],samples=len(samples),current_status=status),indent=2))


if __name__=='__main__':main()
