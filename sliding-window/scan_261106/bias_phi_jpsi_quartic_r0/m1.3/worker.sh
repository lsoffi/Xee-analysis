#!/bin/bash
set -eo pipefail
base=/afs/cern.ch/work/s/soffi/DIELE-WORK/Analisi-statistica-Xee
source /cvmfs/cms.cern.ch/cmsset_default.sh
cd "$base/CMSSW_14_1_0_pre4/src"
eval "$(scram runtime -sh)"
run="$base/scan_261106/bias_phi_jpsi_quartic/m1.3"
m=$1; gen=$2; inj=$3; seed=$4
cd "${_CONDOR_SCRATCH_DIR:?}"
out="$PWD"
combine -M GenerateOnly "/afs/cern.ch/work/s/soffi/DIELE-WORK/Analisi-statistica-Xee/scan_261106/bias_phi_jpsi_quartic/m1.3/${gen}_$m.root" -m "$m" -s "$seed" -t 1 --snapshotName background --setParameters "r=$inj" --toysFrequentist --bypassFrequentistFit --expectSignal "$inj" --saveToys > "$out/generate.log" 2>&1
toy=$(realpath higgsCombine*.GenerateOnly.*.root)
for model in C4 L4; do
 cd "$out"
 for strategy in 0 1; do
  mkdir -p "$out/attempt$strategy/$model"; cd "$out/attempt$strategy/$model"
  combine -M FitDiagnostics "$run/${model}_$m.root" -m "$m" -s "$seed" -t 1 --toysFile "$toy" --toysFrequentist --bypassFrequentistFit --snapshotName background --skipBOnlyFit --minos none --rMin -80 --rMax 80 --cminDefaultMinimizerStrategy "$strategy" --X-rtd MINIMIZER_freezeDisassociatedParams > "$out/$model.log" 2>&1 || true
  cd ..
  if [[ -s "$model/fitDiagnosticsTest.root" ]]; then
   python3 "/afs/cern.ch/work/s/soffi/DIELE-WORK/Analisi-statistica-Xee/scan_261106/bias_phi_jpsi_quartic_r0/m1.3/$model.py" "$run" "$toy" "$gen" "$seed" > "$model.score"
   sed -n 's/^SCORE //p' "$model.score" | python3 -c 'import json,sys; print(json.dumps(json.load(sys.stdin)["branches"][sys.argv[1]]))' "$model" > "$out/$model.json"
   if python3 -c 'import sys,json;assert json.load(open(sys.argv[1]))["ok"]' "$out/$model.json"; then break; fi
  fi
 done
done
cd "$out"
python3 -c 'import sys,json; print(json.dumps(dict(seed=int(sys.argv[1]),generator=sys.argv[2],branches={k:json.load(open(k+".json")) for k in ("C4","L4")})))' "$seed" "$gen" > "/afs/cern.ch/work/s/soffi/DIELE-WORK/Analisi-statistica-Xee/scan_261106/bias_phi_jpsi_quartic_r0/m1.3/results/$seed.tmp"
mv "/afs/cern.ch/work/s/soffi/DIELE-WORK/Analisi-statistica-Xee/scan_261106/bias_phi_jpsi_quartic_r0/m1.3/results/$seed.tmp" "/afs/cern.ch/work/s/soffi/DIELE-WORK/Analisi-statistica-Xee/scan_261106/bias_phi_jpsi_quartic_r0/m1.3/results/$seed.json"
