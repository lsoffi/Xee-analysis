#!/bin/bash
set -eo pipefail
base=/afs/cern.ch/work/s/soffi/DIELE-WORK/Analisi-statistica-Xee
source /cvmfs/cms.cern.ch/cmsset_default.sh
cd "$base/CMSSW_14_1_0_pre4/src";eval "$(scram runtime -sh)"
scan="$base/scan_261106";run="$scan/high_c4l4_limits_bias";m=$1
cd "${_CONDOR_SCRATCH_DIR:?}"
export XEE_DEGREE=4 XEE_C0=1 XEE_PHYSICAL_SR=1
for model in C4 L4;do
 export XEE_BUILD_OUT="$PWD/$model"
 python3 "$scan/bias_phi_jpsi_quartic/$model/build.py" --mass "$m" > "$run/logs/build_${m}_${model}.log" 2>&1
 cd "$model/m$m"
 combine -M FitDiagnostics workspace.root --skipSBFit --minos none --rMin -80 --rMax 80 --cminDefaultMinimizerStrategy 0 --X-rtd MINIMIZER_freezeDisassociatedParams > "$run/logs/prefit_${m}_${model}.log" 2>&1
 python3 "$run/snapshot.py" "$m" "$model" "$run/m$m"
 cd ../..
done
