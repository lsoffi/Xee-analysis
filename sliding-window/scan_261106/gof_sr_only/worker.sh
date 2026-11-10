#!/bin/bash
set -eo pipefail
base=/afs/cern.ch/work/s/soffi/DIELE-WORK/Analisi-statistica-Xee
source /cvmfs/cms.cern.ch/cmsset_default.sh
cd "$base/CMSSW_14_1_0_pre4/src"
eval "$(scram runtime -sh)"
m=$1;kind=$2
run="$base/scan_261106/gof_sr_only"
fit=$(cat "$run/m$m/fit.txt")
evaluate=$(cat "$run/m$m/eval.txt")
cd "${_CONDOR_SCRATCH_DIR:?}"
opts=()
if [[ "$kind" == toys ]]; then opts=(-t 500 --toysFreq); fi
combine -M GoodnessOfFit "$run/m$m/workspace.root" --algo saturated --cminDefaultMinimizerStrategy 0 --setParametersForFit "$fit" --setParametersForEval "$evaluate" "${opts[@]}" > "$run/logs/${m}_${kind}.combine" 2>&1
python3 "$run/read.py" "$m" "$kind" "$run/results"
