#!/bin/bash
set -eo pipefail
base=/afs/cern.ch/work/s/soffi/DIELE-WORK/Analisi-statistica-Xee
source /cvmfs/cms.cern.ch/cmsset_default.sh
cd "$base/CMSSW_14_1_0_pre4/src"
eval "$(scram runtime -sh)"
run="$base/scan_261106"
export XEE_BUILD_OUT="$PWD"
cd "${_CONDOR_SCRATCH_DIR:?}"
export XEE_BUILD_OUT="$PWD" XEE_DEGREE=2 XEE_C0=1 XEE_PHYSICAL_SR=1
m="$1"
out="$run/bias_phi_jpsi_r0/L2/results/m$m"
mkdir -p "$out"
python3 "$run/bias_phi_jpsi_r0/L2/build.py" --mass "$m" > "$out/build.log" 2>&1
cd "m$m"
common=(workspace.root -m "$m" --cminDefaultMinimizerStrategy 0 --X-rtd MINIMIZER_freezeDisassociatedParams)
rc=0
for method in FitDiagnostics; do
  extra=()
  if [[ $method == FitDiagnostics ]]; then extra=(--saveNormalizations --saveShapes --keepFailures --rMin -30 --rMax 30 --robustFit 1); fi
  if combine -M "$method" "${common[@]}" "${extra[@]}" > "$out/$method.log" 2>&1; then
    cp higgsCombine*."${method}".*.root "$out/"
    if [[ $method == FitDiagnostics ]]; then cp fitDiagnostics*.root "$out/"; fi
    touch "$out/$method.done"
  else rc=1; fi
done
cp model.json workspace.root "$out/"
exit "$rc"
