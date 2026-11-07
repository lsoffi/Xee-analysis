#!/bin/bash
set -eo pipefail
base=/afs/cern.ch/work/s/soffi/DIELE-WORK/Analisi-statistica-Xee
source /cvmfs/cms.cern.ch/cmsset_default.sh
cd "$base/CMSSW_14_1_0_pre4/src"
eval "$(scram runtime -sh)"
run="$base/scan_261106/bias_low20"
m=$1; gen=$2; inj=$3; seed=$4
out="$run/results/$seed"; mkdir -p "$out"
cd "$_CONDOR_SCRATCH_DIR"
combine -M GenerateOnly "$run/${gen}_$m.root" -m "$m" -s "$seed" -t 1 --snapshotName background --toysFrequentist --bypassFrequentistFit --expectSignal "$inj" --saveToys > "$out/generate.log" 2>&1
toy=$(realpath higgsCombine*.GenerateOnly.*.root)
for model in cheb log; do
 mkdir "$model"; cd "$model"
 combine -M FitDiagnostics "$run/${model}_$m.root" -m "$m" -s "$seed" -t 1 --toysFile "$toy" --toysFrequentist --bypassFrequentistFit --snapshotName background --skipBOnlyFit --minos none --rMin -80 --rMax 80 --cminDefaultMinimizerStrategy 0 --X-rtd MINIMIZER_freezeDisassociatedParams > "$out/$model.log" 2>&1
 cp fitDiagnostics*.root "$out/$model.root"; cd ..
done
touch "$out/done"
