#!/bin/bash
set -eo pipefail
base=/afs/cern.ch/work/s/soffi/DIELE-WORK/Analisi-statistica-Xee
source /cvmfs/cms.cern.ch/cmsset_default.sh
cd "$base/CMSSW_14_1_0_pre4/src"
eval "$(scram runtime -sh)"
run="$base/scan_261106/bias_expected_06_200"
m=$1; gen=$2; inj=$3; seed=$4
cd "${_CONDOR_SCRATCH_DIR:?}"
out="$PWD"
combine -M GenerateOnly "$run/gen_${gen}_$m.root" -m "$m" -s "$seed" -t 1 --snapshotName injection --setParameters "r=$inj" --toysFrequentist --bypassFrequentistFit --expectSignal "$inj" --saveToys > "$out/generate.log" 2>&1
toy=$(realpath higgsCombine*.GenerateOnly.*.root)
model=L2
for strategy in 0 1; do
 mkdir -p "attempt$strategy/L2"
 cd "attempt$strategy/L2"
combine -M FitDiagnostics "$run/${model}_$m.root" -m "$m" -s "$seed" -t 1 --toysFile "$toy" --toysFrequentist --bypassFrequentistFit --snapshotName background --skipBOnlyFit --minos none --rMin -80 --rMax 80 --cminDefaultMinimizerStrategy "$strategy" --X-rtd MINIMIZER_freezeDisassociatedParams > "$out/$model.log" 2>&1 || true
 cd ..
 if [[ -s L2/fitDiagnosticsTest.root ]]; then
 python3 /afs/cern.ch/work/s/soffi/DIELE-WORK/Analisi-statistica-Xee/scan_261106/bias_low_C3L2_paired/expected_6/score.py "$run" "$toy" "$gen" "$seed" > score.log
 sed -n 's/^SCORE //p' score.log > result.json
 if python3 -c 'import json;assert json.load(open("result.json"))["branches"]["L2"]["ok"]'; then break; fi
 fi
 cd "$out"
done
cd "$out"
f=$(find attempt* -name result.json -size +0c | sort | tail -1)
test -n "$f"
cp "$f" "/afs/cern.ch/work/s/soffi/DIELE-WORK/Analisi-statistica-Xee/scan_261106/bias_low_C3L2_paired/expected_6/results/$seed.tmp"
mv "/afs/cern.ch/work/s/soffi/DIELE-WORK/Analisi-statistica-Xee/scan_261106/bias_low_C3L2_paired/expected_6/results/$seed.tmp" "/afs/cern.ch/work/s/soffi/DIELE-WORK/Analisi-statistica-Xee/scan_261106/bias_low_C3L2_paired/expected_6/results/$seed.json"
