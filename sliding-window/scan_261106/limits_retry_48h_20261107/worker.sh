#!/bin/bash
set -eo pipefail
base=/afs/cern.ch/work/s/soffi/DIELE-WORK/Analisi-statistica-Xee
source /cvmfs/cms.cern.ch/cmsset_default.sh
cd "$base/CMSSW_14_1_0_pre4/src"
eval "$(scram runtime -sh)"
run="$base/scan_261106"
export XEE_BUILD_OUT="$PWD"
cd "${_CONDOR_SCRATCH_DIR:?}"
export XEE_BUILD_OUT="$PWD" XEE_DEGREE=6 XEE_C0=1 XEE_PHYSICAL_SR=1
m="$1"
out="$run/results/m$m"
test ! -f "$out/AsymptoticLimits.done" || exit 0
python3 "$run/build.py" --mass "$m" > "$run/limits_retry_48h_20261107/logs/$m.build.log" 2>&1
cd "m$m"
common=(workspace.root -m "$m" --cminDefaultMinimizerStrategy 0 --X-rtd MINIMIZER_freezeDisassociatedParams)
combine -M AsymptoticLimits "${common[@]}" > "$run/limits_retry_48h_20261107/logs/$m.combine.log" 2>&1
files=(higgsCombine*.AsymptoticLimits.*.root)
test "${#files[@]}" -eq 1 && test -s "${files[0]}" || exit 1
cp "${files[0]}" "$out/${files[0]}.retry.tmp"
mv "$out/${files[0]}.retry.tmp" "$out/${files[0]}"
touch "$out/AsymptoticLimits.done"
