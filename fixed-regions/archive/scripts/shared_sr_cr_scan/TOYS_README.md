# Shared-SR toy pipeline

Run directory: `/Users/imb/Documents/Codex/Analisi-statistica-Xee/runs/20261031_shared_sr_cr_scan`.

The new model has 12 channels (three disjoint mass regions × four eras). A global Chebyshev6/Bernstein6 envelope has shared coefficients, family choice and continuum intensity across SRL/SRH; CR uses an independent Chebyshev4. Signal DCB and both Υ components appear in all three regions. Υ means are estimated in a CR-only B fit for each toy, then fixed in the simultaneous fit; yields remain shared and floating.

## Requested ensembles

- Bias: 61 masses from6 to12 GeV, two injections (0 and the new native expected median), 50 toys per cell =6100 point estimates.
- The r=0 bias ensemble is one common set of50 Chebyshev7 B-only datasets, paired across masses. CR-only fits are cached once per dataset. Injected ensembles have50 independent datasets per mass. There are3100 independent bias datasets, not6100.
- Bias truth is a shared global Chebyshev7 plus CR4, fitted B-only to observed data with the observed CR-only means fixed. Three starts agree and pass positivity/EDM audits.
- GoF:100 separate B-only datasets from the fitted envelope model, reused across all61 masses =6100 r-free native saturated scores. These are not the alternate Chebyshev7 bias datasets.
- Profiles: preselected toy1 per mass and injection,122 outcomes. These run after all point fits and GoF scores. No ensemble coverage is measured.

A preselected toy1 with an invalid primary point or invalid CR receives the terminal outcome `skipped_invalid_point_or_CR`, remains one of the 122 planned controls, and is never replaced by another toy. This is a numerical prerequisite, not a selection on the bias or interval size. Covariance quality 2 alone does not cause skipping. Valid controls receive the bounded full profile protocol.

## Numerical protocol

Point estimates use the same two starts, r=0 and the mass-specific expected median, for both injection ensembles. At most one numerical retry per start is allowed. Selection uses unshifted extended deviance plus auxiliary penalties, not offset `minNll`. Validity requires status0, finite EDM<0.01, positive continuum and total expectations, no r boundary, and a valid CR fit. Hesse quality is a separate diagnostic; covQual2 does not by itself discard a converged point estimate.

The extended deviance uses `2*(expectedEvents - Nobs + sum(n*log(n/mu_at_bin_center)))`. Do not replace the analytic expectedEvents term by the sum of bin-center approximations.

Failed toys are never replaced. All finite results, numerical failures, boundaries and CR mean shifts are retained. Point bias statistics are conditional on numerical fit validity, with an all-finite sensitivity shown separately.

## Outputs

- `bias/fits/mM/{r0,rExpected}/toyNNN/result.json`: primary two-anchor point result.
- `bias/profiles/mM/{r0,rExpected}/toy001/result.json`: separate control interval result; primary point estimate remains unchanged. A lower profile minimum is flagged.
- `bias/generation/common_r0`: common50 B-only bias toys.
- `bias/generation/mM/rExpected`: injected toy ensembles.
- `bias/common_r0_cr/toyNNN`: shared CR-only fits for the common B-only bias toys.
- `gof/generation`:100 common envelope B-only GoF datasets.
- `gof/toyNNN/mM/result.json`: native score, fitted r and quality flags.
- `toy_templates/mM/workspace.root`: immutable private copies of the observed B-only seed workspaces.
- `bias/anchors/mM.json`: immutable expected-median anchors; unrelated later additions to observed results do not change injections.
- `toys_progress_main.json`, `bias/profiles/progress.json`, `toys_audit.json`: progress and audits.

Standard GoodnessOfFit does not expose nominal-fit EDM and performs a saturated minimization after the nominal fit. The isolated immutable runtime in `work/shared_gof_diagnostic_runtime/v2` adds diagnostic-only double-precision branches captured before the saturated fit. Exact native-score parity has been verified. `gof_diagnostics.read_nominal_diagnostics` checks status, nonnegative EDM<0.01, saved-state agreement, full floating-parameter completeness, exact polynomial and total-bin positivity, and reconstructed deviance plus the constant envelope penalty. No extra Hesse is run, so covQual3 is not required. Earlier raw scores are preserved and audited in separate `nominal_diagnostic` subdirectories; pending diagnostics must remain explicit.

Every completed GoF job has a terminal `nominal_diagnostics.outcome_complete` flag, including jobs with no usable score, CR failures, timeouts and parse errors. These failed outcomes preserve the missing score rather than inventing a number. A scored legacy record awaiting diagnostic backfill is different: it remains pending until audited. The lightweight `toys_qc_outcomes.py --watch` finalizer covers crash fallbacks from the already-running scheduler without restarting fits. Completion requires 6100 point outcomes, 6100 GoF outcomes with 6100 completed QC outcomes, and 122 control-profile outcomes.

## Running and monitoring

Use `bin/python` for ROOT scripts. `toys_run.py` itself uses system Python and starts ROOT subprocesses. The active scheduler PID/session is recorded in `scheduler_process.json`; check it before launching another feeder.

```sh
python3 scripts/shared_sr_cr_scan/toys_run.py --workers 3 --worker-budget-file /Users/imb/Documents/Codex/Analisi-statistica-Xee/runs/20261031_shared_sr_cr_scan/toys_worker_budget.json --name main
```

The scheduler reads `toys_worker_budget.json` every3 seconds. Change its `workers` field to adjust concurrency without restarting. It waits for each observed median, generates required ensembles, runs point-only fits, and automatically starts the final profile queue. Existing results and locked jobs are skipped.

Do not terminate active fit process groups to change scheduling. The earlier feeder processes were stopped separately, leaving their active control profiles to finish; all four initial control outcomes were retained or are migrated by `toys_refresh_edm.py`. This migration is read-only with respect to ROOT fits and archives earlier JSON.

`caffeinate -i -w <scheduler PID>` prevents idle system sleep only for the scheduler lifetime. It does not change screen locking or display settings.

`toys_audit.py` verifies immutable inputs and samples exact toy cloning, auxiliary observations and fixed CR means. `--all` expands semantic checks to every point result and may take considerable time. The lightweight EDM checks cover all point records regardless of the sampling option.
