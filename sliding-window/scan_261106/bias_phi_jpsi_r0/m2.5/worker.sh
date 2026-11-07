#!/bin/bash
set -eo pipefail
base=/afs/cern.ch/work/s/soffi/DIELE-WORK/Analisi-statistica-Xee
source /cvmfs/cms.cern.ch/cmsset_default.sh
cd "$base/CMSSW_14_1_0_pre4/src"
eval "$(scram runtime -sh)"
run="$base/scan_261106/bias_phi_jpsi_r0/m2.5"
m=$1; gen=$2; inj=$3; seed=$4
out="$run/results/$seed"; mkdir -p "$out"
cd "$_CONDOR_SCRATCH_DIR"
combine -M GenerateOnly "$run/${gen}_$m.root" -m "$m" -s "$seed" -t 1 --snapshotName background --toysFrequentist --bypassFrequentistFit --expectSignal "$inj" --saveToys > "$out/generate.log" 2>&1
toy=$(realpath higgsCombine*.GenerateOnly.*.root)
for model in C2 L3; do
 mkdir "$model"; cd "$model"
 combine -M FitDiagnostics "$run/${model}_$m.root" -m "$m" -s "$seed" -t 1 --toysFile "$toy" --toysFrequentist --bypassFrequentistFit --snapshotName background --skipBOnlyFit --minos none --rMin -80 --rMax 80 --cminDefaultMinimizerStrategy 0 --X-rtd MINIMIZER_freezeDisassociatedParams > "$out/$model.log" 2>&1
 cd ..
done
python3 "$run/score.py" "$run" "$toy" "$gen" "$seed" > "$out/score.log"
grep "^SCORE " "$out/score.log" | cut -c7- > "$out/result.tmp"
mv "$out/result.tmp" "$out/result.json"
