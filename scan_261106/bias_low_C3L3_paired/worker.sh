#!/bin/bash
set -eu
config=$1; shift
exec bash "/afs/cern.ch/work/s/soffi/DIELE-WORK/Analisi-statistica-Xee/scan_261106/bias_low_C3L3_paired/$config/worker.sh" "$@"
