# Analysis-note update — 23 November 2026

The main note now has 46 pages. Pages 41–42 document the current full-spectrum configuration; pages 43–46 add the completed 9 GeV CR LogPol2/3 comparison, diagnostics and numerical recovery. Earlier pages remain historical, and the updated comparison explicitly supersedes pages 34/40.

## Rebuild

Requires Python with pypdf and reportlab (pypdfium2 only for optional rendering).
From repository root:

```bash
git show 932e121:sliding-window/docs/Xee_analysis_note.pdf > /tmp/xee-note-base.pdf
python sliding-window/docs/updates_20261123/current_configuration_source.py
python sliding-window/docs/updates_20261123/build.py --base /tmp/xee-note-base.pdf --output sliding-window/docs/Xee_analysis_note.pdf
```

## Evidence and scope

`evidence/` contains the four native limit summaries, five-configuration CSV, FitDiagnostics status/EDM/covariance and r estimates extracted from the saved ROOT files, and graph-equivalence audits. `figures/` contains the updated limit comparison and B-only/S+B prefit/postfit diagnostics.

Final remote results: `scan_261106/m90_cr_logpol{2,3}_c1/limits_fdstart_completion/{expected,observed}`. Fit outputs: `fit_diagnostics_fast`. All four final limit receipts report exit zero on 19 November 2026. Earlier timeout attempts are not the source of the reported limits.

Historical note audit: low/middle/high bias for both injections, GoF, coverage, fixed-trigger diagnostics, initialization checks and HybridNew were already described in the baseline note. Preserve their original completeness and validity qualifications. No new toy validation is claimed for the 23 November configuration. The remaining historical partial-result statements have not been promoted to completed results without output verification.
