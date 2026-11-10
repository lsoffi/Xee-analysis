# Xee sliding-window scan

Sliding windows and selected control regions by mass · [PDF](docs/tabella_SR_CR_5MeV.pdf) · [CSV](docs/tabella_SR_CR_5MeV.csv)

Snapshot of the four-era CERN production scan for 69 mass hypotheses. The production worker sets Chebyshev degree 6, enables C0 at shared boundaries and leaves the slope free. The per-mass signal windows and selected control regions are in scan_261106/records.json.

The scripts use the CERN production paths and CMSSW_14_1_0_pre4 with Combine. Observed event counts (bins.csv), signal cards and ROOT workspaces (signal/), and job outputs are external. Run from the configured production environment with those inputs available. The repository also includes the limit retry submission and the toy study workflows described below.

## Reproduce the nominal scan on CERN infrastructure

Use Bash on lxplus, an AFS work area, access to the EOS inputs, and a CMSSW_14_1_0_pre4 area with CombinedLimit v10.6.1. Run the commands below in a separate checkout. Several production scripts contain absolute AFS and EOS paths. Before using Condor or publishing plots from a new checkout, point those paths at your own checkout, CMSSW installation and output directory. The original jobs write into the production results directory.

1. Activate the environment with `source /cvmfs/cms.cern.ch/cmsset_default.sh`, enter the CMSSW release's `src` directory and run `eval "$(scram runtime -sh)"`. Verify that `combine --version` reports v10.6.1.
2. From the checkout's `scan_261106` directory, copy `bins.csv` and the `signal/` directory from the CERN production scan. Run `sha256sum -c ../provenance/scan-inputs.sha256`. The manifest covers the observed counts and the individual signal workspaces. The committed `records.json` and `cr/` supply the actual windows, control regions and fitted shapes.
3. To regenerate those inputs, inspect `prepare.py`. It reads the 5 MeV window table from EOS, the four data workspaces under `/eos/user/n/npalmeri/shared/datasets_261106/261106`, the signal inputs under `/eos/user/n/npalmeri/shared/datacards_261021/261021`, and the era-specific CR fit JSON files from the sibling production directory `cr_261106`. Run it only in the reproduction checkout and recheck the hashes.

For a single nominal point, execute the following from `scan_261106` after restoring the inputs.

```bash
export XEE_BUILD_OUT="$PWD/reproduction" XEE_DEGREE=6 XEE_C0=1 XEE_PHYSICAL_SR=1
python3 build.py --mass 7.4
cd "$XEE_BUILD_OUT/m7.4"
common=(workspace.root -m 7.4 --cminDefaultMinimizerStrategy 0 --X-rtd MINIMIZER_freezeDisassociatedParams)
combine -M FitDiagnostics "${common[@]}" --saveNormalizations --saveShapes --keepFailures --rMin -30 --rMax 30 --robustFit 1
combine -M AsymptoticLimits "${common[@]}"
combine -M Significance "${common[@]}"
```

For the 69-point production configuration, create the logs directory and run `condor_submit jobs.sub` from the staged scan directory after adapting absolute paths. The worker stages each point in Condor scratch. Inspect native fit status and covariance. For limits, require a finite observed value and all five ordered expected quantiles. The plotting script `plot.py` reads the method receipts and ROOT trees, then writes the limit and local p-value plots. Its production entry point polls for up to 720 minutes.

## Toy studies

The bias workflow directories preserve the available generation, fit, submission and report scripts, including mass-specific job lists and injection manifests where present. For the quartic study between phi and J/psi, inspect `bias_phi_jpsi_quartic/study.dag` and the sibling `bias_phi_jpsi_quartic_r0` and `bias_phi_jpsi_quartic_expected` directories. The DAG orders background fits, snapshot preparation, toys and reports. The production design uses 200 toys per mass, generator and injection, with C4 and L4 fitted to each toy. The two branches carry the same score penalty of 16. Preserve the original seed lists and record exclusions when validating the bias.

This is a source snapshot of an ongoing study. Generated snapshots, ROOT results, job receipts and published plots live in the CERN production area. The L4 prefit at 2.6 GeV was undergoing recovery when this snapshot was made. A complete replay from a fresh checkout has not yet been verified.

## Fixed-trigger validation update, 9 November 2026

The [analysis note](docs/Xee_analysis_note.pdf), section 5.2, includes the C4+L4 median-injection bias with `CMS_eff_trigger_e` fixed to zero across 1.2–2.7 GeV. All 6400 toys pass the existing two-branch validity criterion. Thirty of 32 central estimates satisfy |bias| <= 0.2. The [summary CSV](docs/data/bias_trigger_fixed_all_expected.csv) includes the bootstrap intervals and branch counts.

Run `python3 bias_trigger_fixed_all_expected/plot.py` in the CERN production environment to aggregate the two fixed-trigger result directories and publish the plot. The script uses the recorded mass, generator, injection and seed manifests. Run `python3 docs/update_trigger_note.py` from this directory with pypdf, reportlab and Pillow installed to regenerate the final two PDF pages from the committed figures, preserving the first seven pages.

The coverage workspace builder registers the discrete envelope category as a generic `discreteParams` RooArgSet, as required by Combine, and excludes that category from the continuous nuisance constraint set. The trigger remains profiled. Coverage results are pending validation.


## Status — 10 November 2026

The middle-region C4+L4 GoF calculation is complete: 8,000 toys, one Condor job per toy. See [method, commands and reproduction](scan_261106/gof_sr_single_toy/README.md), [results](docs/data/gof_middle_sr.csv), and section 5.4 of the [analysis note](docs/Xee_analysis_note.pdf). Thirteen of sixteen masses have local p ≥ 0.05; the qualifications at 2.2, 2.6 and 2.7 GeV remain documented. No additional middle-region production is scheduled by this update.

The active work is the **high region (4.1–9.0 GeV)**: `high_c4l4_diagnostics` (fits), `high_c4l4_limits_bias` (limits and 40,000 planned bias toys), and `upsilon_limit_diagnostics` (profiles and CLs checks at 8.6/9.0). These workflows remain provisional, with missing results and fit-quality flags. Native envelope FitDiagnostics uses `--cminRunAllDiscreteCombinations` in addition to the original settings to avoid the reproduced iterative-minimizer crash. `partial_reports/report.py` publishes available results without replacing missing points with zeros. `task_monitor/monitor.py` is a run-specific read-only dashboard, not a general scheduler.

Persistent execution rule: one Condor job per toy, unique reproducible seeds, immediate result and exit receipts. Reuse completed results and never duplicate active seeds. ROOT run products are not committed. The historical `high_c4l4_diagnostics/worker.sh` records the failed first attempt; `retry.sh` is the recovery command. Site paths and cluster IDs in monitoring scripts describe this production and must be updated for a new campaign.
