# Xee sliding-window scan

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
