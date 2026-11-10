#!/bin/bash
set -eo pipefail
base=/afs/cern.ch/work/s/soffi/DIELE-WORK/Analisi-statistica-Xee
source /cvmfs/cms.cern.ch/cmsset_default.sh
cd "$base/CMSSW_14_1_0_pre4/src";eval "$(scram runtime -sh)"
run="$base/scan_261106/upsilon_limit_diagnostics"
m=$1;kind=$2;point=$3;start=$4;tag=$5
w="$base/scan_261106/high_c4l4_diagnostics/m$m/envelope.root"
out="$run/results/$tag";mkdir -p "$out"
cd "${_CONDOR_SCRATCH_DIR:?}"
trap 'rc=$?; printf "%s\n" "$rc" > "$out/exit_code.txt"; cp -f combine_logger.out "$out/internal.log" 2>/dev/null || true' EXIT
common=(--cminDefaultMinimizerStrategy 0 --X-rtd MINIMIZER_freezeDisassociatedParams --cminRunAllDiscreteCombinations --rMin 0 --rMax 200 --setParameters "r=$start")
if [[ "$kind" == profile ]];then
 nuis=$(cat "$run/nuis_$m.txt");funcs=$(cat "$run/funcs_$m.txt")
 combine -M MultiDimFit "$w" --algo fixed --fixedPointPOIs "r=$point" --saveNLL --saveFitResult --saveWorkspace --saveSpecifiedNuis "$nuis" --saveSpecifiedFunc "$funcs" --saveSpecifiedIndex pdfindex "${common[@]}" > "$out/combine.log" 2>&1
else
 combine -M AsymptoticLimits "$w" --singlePoint "$point" "${common[@]}" > "$out/combine.log" 2>&1
fi
cp ./*.root "$out/"
python3 "$run/extract.py" "$m" "$kind" "$point" "$start" "$tag"
