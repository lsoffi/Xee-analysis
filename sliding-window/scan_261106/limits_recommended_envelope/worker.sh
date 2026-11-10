#!/bin/bash
set -eo pipefail
base=/afs/cern.ch/work/s/soffi/DIELE-WORK/Analisi-statistica-Xee
source /cvmfs/cms.cern.ch/cmsset_default.sh
cd "$base/CMSSW_14_1_0_pre4/src"
eval "$(scram runtime -sh)"
m=$1
workspace=$2
run="$base/scan_261106/limits_recommended_envelope"
cd "${_CONDOR_SCRATCH_DIR:?}"
combine -M AsymptoticLimits "$workspace" --cminDefaultMinimizerStrategy 0 --X-rtd MINIMIZER_freezeDisassociatedParams > "$run/logs/m${m}.combine" 2>&1
python3 "$run/read.py" "$m" "$run/results"
