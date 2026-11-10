#!/bin/bash
set -eo pipefail
base=/afs/cern.ch/work/s/soffi/DIELE-WORK/Analisi-statistica-Xee
source /cvmfs/cms.cern.ch/cmsset_default.sh
cd "$base/CMSSW_14_1_0_pre4/src"; eval "$(scram runtime -sh)"
run="$base/scan_261106/high_c4l4_limits_bias";m=$1
cd "${_CONDOR_SCRATCH_DIR:?}"
combine -M AsymptoticLimits "$base/scan_261106/high_c4l4_diagnostics/m$m/envelope.root" --cminDefaultMinimizerStrategy 0 --X-rtd MINIMIZER_freezeDisassociatedParams --cminRunAllDiscreteCombinations > "$run/logs/limit_$m.log" 2>&1
python3 "$run/read_limit.py" "$m" "$run/limits"
