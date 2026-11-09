#!/bin/bash
set -eo pipefail
source /cvmfs/cms.cern.ch/cmsset_default.sh
cd /afs/cern.ch/work/s/soffi/DIELE-WORK/Analisi-statistica-Xee/CMSSW_14_1_0_pre4/src
eval "$(scram runtime -sh)"
cd /afs/cern.ch/work/s/soffi/DIELE-WORK/Analisi-statistica-Xee/scan_261106/coverage_all_env/m${1}
python3 ../build.py "$1"
