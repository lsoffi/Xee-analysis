#!/bin/bash
set -eu
main="$PWD/main"
mkdir -p "$main" "$PWD/retry"
start="$PWD"
export _CONDOR_SCRATCH_DIR="$main"
bash /afs/cern.ch/work/s/soffi/DIELE-WORK/Analisi-statistica-Xee/scan_261106/bias_phi_jpsi_r0/m1.3/worker.sh "$@" || true
if ! python3 -c 'import json,sys; assert all(v["ok"] for v in json.load(open(sys.argv[1]))["branches"].values())' "/afs/cern.ch/work/s/soffi/DIELE-WORK/Analisi-statistica-Xee/scan_261106/bias_phi_jpsi_r0/m1.3/results/$4/result.json"; then
 export _CONDOR_SCRATCH_DIR="$start/retry"
 bash /afs/cern.ch/work/s/soffi/DIELE-WORK/Analisi-statistica-Xee/scan_261106/bias_phi_jpsi_r0/m1.3/retry.sh "$@"
fi
