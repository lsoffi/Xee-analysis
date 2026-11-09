#!/bin/bash
set -eo pipefail
base=/afs/cern.ch/work/s/soffi/DIELE-WORK/Analisi-statistica-Xee
source /cvmfs/cms.cern.ch/cmsset_default.sh
cd "$base/CMSSW_14_1_0_pre4/src"
eval "$(scram runtime -sh)"
run="$base/scan_261106/bias_phi_jpsi_quartic/m1.4"
m=$1; gen=$2; inj=$3; seed=$4
cd "${_CONDOR_SCRATCH_DIR:?}"
out="$PWD"
combine -M GenerateOnly "/afs/cern.ch/work/s/soffi/DIELE-WORK/Analisi-statistica-Xee/scan_261106/bias_phi_jpsi_quartic/m${m}/gen_${gen}.root" -m "$m" -s "$seed" -t 1 --snapshotName injection --setParameters "r=$inj" --toysFrequentist --bypassFrequentistFit --expectSignal "$inj" --saveToys > "$out/generate.log" 2>&1
toy=$(realpath higgsCombine*.GenerateOnly.*.root)
combine -M AsymptoticLimits "/afs/cern.ch/work/s/soffi/DIELE-WORK/Analisi-statistica-Xee/scan_261106/coverage_all_env/m${m}/envelope.root" -m "$m" -s "$seed" -t 1 --toysFile "$toy" --toysFrequentist --bypassFrequentistFit --singlePoint "$inj" --rMax 150 --cminDefaultMinimizerStrategy 0 --X-rtd MINIMIZER_freezeDisassociatedParams > "/afs/cern.ch/work/s/soffi/DIELE-WORK/Analisi-statistica-Xee/scan_261106/coverage_all_env/m${m}/logs/$seed.combine" 2>&1
python3 "/afs/cern.ch/work/s/soffi/DIELE-WORK/Analisi-statistica-Xee/scan_261106/coverage_all_env/read.py" "$seed" "$gen" "$inj" "/afs/cern.ch/work/s/soffi/DIELE-WORK/Analisi-statistica-Xee/scan_261106/coverage_all_env/m${m}/results"
