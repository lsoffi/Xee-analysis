# Middle-region SR-only goodness of fit

Completed production: 16 masses (1.2–2.7 GeV), 500 toys per mass, 8,000 unique seeds. The summary and figure are committed under `../../docs/data/gof_middle_sr.csv` and `../../docs/figures/gof_middle_sr.png`; the analysis note describes the method in section 5.4, pages 12–13.

## Statistical procedure

Fit the simultaneous SR + selected CR likelihood in four eras, profiling `r`, the shared C4/L4 index and the nuisance parameters including trigger efficiency. Use the existing total-C0 model and table geometry. All channel masks are zero for minimisation. CR masks become one only when evaluating the native saturated GoF statistic; SR masks remain zero. Combine applies these settings in both nominal and saturated stages. This retains the CR constraint on the fit but removes its direct channel contribution at statistic evaluation. Shared auxiliary constraints retain Combine's native treatment; this is not a separately computed Pearson chi-square.

Each independent toy command is:

```bash
combine -M GoodnessOfFit "$workspace" --algo saturated \
  --cminDefaultMinimizerStrategy 0 \
  --setParametersForFit "$fit_masks" \
  --setParametersForEval "$evaluation_masks" \
  -t 1 --toysFreq -s "$seed"
```

Observed data use the same command without `-t`, `--toysFreq`, and `-s`. No additional r-range option is imposed. The old 500-toy sequential batch campaign was stopped; its incomplete outputs are not pooled into this sample. `rows.txt` and `manifest.json` specify the 8,000 unique seeds: `91000000 + 1000*i + j`, with `i=12..27` and `j=0..499`.

The p-value is the fraction of toys whose statistic is at least the observed statistic. Bars are exact two-sided Clopper–Pearson 95% binomial intervals. These are pointwise tests of correlated, overlapping windows, not a globally calibrated scan test.

## Reproduction on CERN

The scripts preserve the actual production paths, rooted at `/afs/cern.ch/work/s/soffi/DIELE-WORK/Analisi-statistica-Xee`. Deploy the versioned `scan_261106` code into that working area; activate `CMSSW_14_1_0_pre4` with Combine v10.6.1. ROOT workspaces, datasets and generated toy ROOT files are external run products, not Git inputs. See the parent README and coverage builders for input provenance.

1. Rebuild the unmasked envelope workspaces with `coverage_all_env/build.py` (the 1.4 GeV source is `coverage_m14_env/envelope.root`). These depend on the quartic branch snapshots described in the earlier bias/coverage workflow.
2. Run `gof_sr_only/build.py MASS` for each mass. It writes the masked workspace, ordered channel masks, and an audit JSON. The `gof_sr_only/study.dag` is historical: it launches sequential 500-toy jobs and must **not** be resubmitted for new toy production.
3. Reuse or calculate the one observed statistic per mass with `gof_sr_only/worker.sh MASS observed` under Condor (the worker requires `_CONDOR_SCRATCH_DIR`).
4. Submit `gof_sr_single_toy/jobs.sub` from its directory using `condor_submit -name bigbird18.cern.ch jobs.sub`, only in a fresh run or after filtering the manifest to missing seeds without active jobs. Every row is one toy. Do not submit an entire existing manifest again.
5. Run `python3 report.py --output OUTPUT_DIR`. Optional `--results` and `--observed` arguments locate external result folders. It validates all 500 unique seeds per mass and regenerates CSV/PNG/PDF.

## Outcome and scope

13/16 masses have p ≥ 0.05. Local low p-values remain at 2.2 (0.024), 2.6 (0.034) and 2.7 GeV (0.028). C4+L4 is the adopted working configuration, not a uniform GoF pass. The middle-region computation is complete; the current production and model-development work is in the high region. These qualifications are retained, not declared resolved. GoF does not replace bias/coverage studies or validate CR shapes.
