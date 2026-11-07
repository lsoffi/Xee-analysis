#!/bin/bash
set -eo pipefail
source /cvmfs/cms.cern.ch/cmsset_default.sh
cd /afs/cern.ch/work/s/soffi/DIELE-WORK/Analisi-statistica-Xee/CMSSW_14_1_0_pre4/src
eval "$(scram runtime -sh)"
cd /afs/cern.ch/work/s/soffi/DIELE-WORK/Analisi-statistica-Xee/scan_261106/bias_env_06_200
python3 report_final.py
