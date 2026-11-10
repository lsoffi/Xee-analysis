#!/bin/bash
set -eo pipefail
base=/afs/cern.ch/work/s/soffi/DIELE-WORK/Analisi-statistica-Xee
source /cvmfs/cms.cern.ch/cmsset_default.sh
cd "$base/CMSSW_14_1_0_pre4/src"
eval "$(scram runtime -sh)"
scan="$base/scan_261106"
run="$scan/high_c4l4_diagnostics"
m="$1"
out="$run/m$m"
mkdir -p "$out"
cd "${_CONDOR_SCRATCH_DIR:?}"
export XEE_DEGREE=4 XEE_C0=1 XEE_PHYSICAL_SR=1
for model in C4 L4; do
 export XEE_BUILD_OUT="$PWD/$model"
 python3 "$scan/bias_phi_jpsi_quartic/$model/build.py" --mass "$m" > "$out/build_$model.log" 2>&1
 cp "$model/m$m/model.json" "$model.json"
 python3 - "$model/m$m/workspace.root" "${model}_${m}.root" <<'PY'
import ROOT as R,sys
R.gSystem.Load('libHiggsAnalysisCombinedLimit')
f=R.TFile.Open(sys.argv[1]);w=f.Get('w');w.saveSnapshot('background',w.allVars());w.writeToFile(sys.argv[2])
PY
done
python3 "$run/envelope.py" "$m" > "$out/envelope.log" 2>&1
cp envelope.root C4.json L4.json "$out/"
combine -M FitDiagnostics envelope.root --saveNormalizations --saveShapes --keepFailures --rMin -30 --rMax 30 --cminDefaultMinimizerStrategy 0 --robustFit 1 --X-rtd MINIMIZER_freezeDisassociatedParams > "$out/combine.log" 2>&1
cp fitDiagnostics*.root "$out/"
python3 "$run/report.py" "$m"
