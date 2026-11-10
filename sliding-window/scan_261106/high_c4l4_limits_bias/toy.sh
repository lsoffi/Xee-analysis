#!/bin/bash
set -eo pipefail
base=/afs/cern.ch/work/s/soffi/DIELE-WORK/Analisi-statistica-Xee
source /cvmfs/cms.cern.ch/cmsset_default.sh
cd "$base/CMSSW_14_1_0_pre4/src";eval "$(scram runtime -sh)"
root="$base/scan_261106/high_c4l4_limits_bias"
m=$1;gen=$2;tag=$3;seed=$4;run="$root/m$m"
inj=$(python3 - "$root/limits/m$m.json" "$tag" <<'PY'
import json,sys
v=next(x['limit'] for x in json.load(open(sys.argv[1]))['quantiles'] if abs(x['quantile']-.5)<1e-5)
assert 0<v<80, 'Expected injection outside validated toy range; requires range extension'
print(v if sys.argv[2]=='expected' else 0)
PY
)
cd "${_CONDOR_SCRATCH_DIR:?}";out="$PWD"
trap 'rc=$?; printf "%s\n" "$rc" > "$run/results/${seed}.exit"' EXIT
combine -M GenerateOnly "$run/gen_$gen.root" -m "$m" -s "$seed" -t 1 --snapshotName injection --setParameters "r=$inj" --toysFrequentist --bypassFrequentistFit --expectSignal "$inj" --saveToys > generate.log 2>&1
toy=$(realpath higgsCombine*.GenerateOnly.*.root)
for model in C4 L4;do
 for strategy in 0 1;do
  mkdir -p "$out/attempt$strategy/$model";cd "$out/attempt$strategy/$model"
  combine -M FitDiagnostics "$run/${model}_$m.root" -m "$m" -s "$seed" -t 1 --toysFile "$toy" --toysFrequentist --bypassFrequentistFit --snapshotName background --skipBOnlyFit --minos none --rMin -80 --rMax 80 --cminDefaultMinimizerStrategy "$strategy" --X-rtd MINIMIZER_freezeDisassociatedParams > "$out/$model.log" 2>&1 || true
  cd ..
  if [[ -s "$model/fitDiagnosticsTest.root" ]];then
   python3 "$run/$model.py" "$run" "$toy" "$gen" "$seed" > "$model.score"
   sed -n 's/^SCORE //p' "$model.score" | python3 -c 'import json,sys;print(json.dumps(json.load(sys.stdin)["branches"][sys.argv[1]]))' "$model" > "$out/$model.json"
   if python3 -c 'import json,sys;assert json.load(open(sys.argv[1]))["ok"]' "$out/$model.json";then break;fi
  fi
 done
done
cd "$out"
python3 - "$m" "$gen" "$tag" "$seed" "$inj" "$run" <<'PY'
import json,sys,pathlib
m,g,tag,seed,inj,run=sys.argv[1:];b={k:json.load(open(k+'.json')) for k in ['C4','L4']};valid=all(v['ok'] for v in b.values());winner=min(b,key=lambda k:b[k]['score']) if valid else None
p=pathlib.Path(run)/'results'/f'{seed}.json';tmp=p.with_suffix('.tmp');tmp.write_text(json.dumps(dict(mass=float(m),generator=g,tag=tag,seed=int(seed),injection=float(inj),branches=b,valid=valid,winner=winner)));tmp.replace(p)
PY
