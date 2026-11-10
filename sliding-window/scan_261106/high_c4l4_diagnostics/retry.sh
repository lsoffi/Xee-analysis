#!/bin/bash
set -eo pipefail
base=/afs/cern.ch/work/s/soffi/DIELE-WORK/Analisi-statistica-Xee
source /cvmfs/cms.cern.ch/cmsset_default.sh
cd "$base/CMSSW_14_1_0_pre4/src"
eval "$(scram runtime -sh)"
run="$base/scan_261106/high_c4l4_diagnostics"
m="$1";out="$run/m$m"
cd "${_CONDOR_SCRATCH_DIR:?}"
trap 'rc=$?; printf "%s\n" "$rc" > "$out/retry_exit_code.txt"; cp -f combine_logger.out "$out/retry_internal.log" 2>/dev/null || true' EXIT
combine -M FitDiagnostics "$out/envelope.root" --saveNormalizations --saveShapes --keepFailures --rMin -30 --rMax 30 --cminDefaultMinimizerStrategy 0 --robustFit 1 --X-rtd MINIMIZER_freezeDisassociatedParams --cminRunAllDiscreteCombinations > "$out/retry_combine.log" 2>&1
cp fitDiagnostics*.root "$out/"
python3 "$run/report.py" "$m"
