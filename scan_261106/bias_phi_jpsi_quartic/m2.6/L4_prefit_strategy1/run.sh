#!/bin/bash
set -eo pipefail
base=/afs/cern.ch/work/s/soffi/DIELE-WORK/Analisi-statistica-Xee
source /cvmfs/cms.cern.ch/cmsset_default.sh
cd "$base/CMSSW_14_1_0_pre4/src"
eval "$(scram runtime -sh)"
d="$base/scan_261106/bias_phi_jpsi_quartic"
cd "$d/m2.6/L4_prefit_strategy1"
combine -M FitDiagnostics "$d/L4/results/m2.6/workspace.root" -m 2.6 --cminDefaultMinimizerStrategy 1 --X-rtd MINIMIZER_freezeDisassociatedParams --saveNormalizations --saveShapes --keepFailures --rMin -30 --rMax 30 --robustFit 1
