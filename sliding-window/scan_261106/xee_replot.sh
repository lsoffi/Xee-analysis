#!/bin/bash
set -e
base=/afs/cern.ch/work/s/soffi/DIELE-WORK/Analisi-statistica-Xee
source /cvmfs/cms.cern.ch/cmsset_default.sh
cd "$base/CMSSW_14_1_0_pre4/src"
eval "$(scram runtime -sh)"
cd "$base/scan_261106"
sed -e "s/name+'.tmp.png'/name+'.refresh2.tmp.png'/" -e "s/'results.tmp'/'results.refresh2.tmp'/" plot.py > plot_refresh2.py
nohup python3 -u plot_refresh2.py > plot_refresh2.log 2>&1 < /dev/null &
echo "monitor $!"
