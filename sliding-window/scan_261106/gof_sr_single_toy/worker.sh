#!/bin/bash
set -eo pipefail
base=/afs/cern.ch/work/s/soffi/DIELE-WORK/Analisi-statistica-Xee
source /cvmfs/cms.cern.ch/cmsset_default.sh
cd "$base/CMSSW_14_1_0_pre4/src"
eval "$(scram runtime -sh)"
m=$1;seed=$2
source_run="$base/scan_261106/gof_sr_only"
run="$base/scan_261106/gof_sr_single_toy"
fit=$(cat "$source_run/m$m/fit.txt")
evaluate=$(cat "$source_run/m$m/eval.txt")
cd "${_CONDOR_SCRATCH_DIR:?}"
trap 'rc=$?; printf "%s\n" "$rc" > "$run/logs/${m}_${seed}.exit"' EXIT
combine -M GoodnessOfFit "$source_run/m$m/workspace.root" --algo saturated --cminDefaultMinimizerStrategy 0 --setParametersForFit "$fit" --setParametersForEval "$evaluate" -t 1 --toysFreq -s "$seed" > "$run/logs/${m}_${seed}.combine" 2>&1
python3 "$run/read.py" "$m" "$seed" "$run/results"
