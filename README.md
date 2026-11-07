# Xee sliding-window scan

Snapshot of the four-era CERN production scan for 69 mass hypotheses. The production worker sets Chebyshev degree 6, enables C0 at shared boundaries and leaves the slope free. The per-mass signal windows and selected control regions are in scan_261106/records.json.

The scripts use the CERN production paths and CMSSW_14_1_0_pre4 with Combine. Observed event counts (bins.csv), signal cards and ROOT workspaces (signal/), and job outputs are external. Run from the configured production environment with those inputs available. This snapshot covers the nominal scan and the limit retry submission. The separate toy bias campaigns are outside its scope.
