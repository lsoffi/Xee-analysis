#!/usr/bin/env bash
set -euo pipefail
XEE_REPO="$(cd "$(dirname "$0")/.." && pwd)"
XEE_COMBINE_COMMIT=6c543b6c0a29df51d390652c2b0ff559c8a211f1
command -v pixi >/dev/null || { echo 'Install pixi first: https://pixi.sh'; exit 1; }
mkdir -p "$XEE_REPO/runtime"
if [[ ! -d "$XEE_REPO/runtime/CombinedLimit/.git" ]]; then
  git clone https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit.git "$XEE_REPO/runtime/CombinedLimit"
fi
git -C "$XEE_REPO/runtime/CombinedLimit" checkout --detach "$XEE_COMBINE_COMMIT"
cp "$XEE_REPO/environment/combine-pixi.toml" "$XEE_REPO/runtime/CombinedLimit/pixi.toml"
cp "$XEE_REPO/environment/combine-pixi.lock" "$XEE_REPO/runtime/CombinedLimit/pixi.lock"
(cd "$XEE_REPO/runtime/CombinedLimit" && pixi run --frozen install)
# Use the locked Python version; an arbitrary system Python can be too new for
# the pinned plotting wheels or absent on a fresh computer.
pixi run --frozen --manifest-path "$XEE_REPO/runtime/CombinedLimit/pixi.toml" \
  python -m venv "$XEE_REPO/.venv-plots"
"$XEE_REPO/.venv-plots/bin/python" -m pip install -r "$XEE_REPO/environment/plot-requirements.txt"
"$XEE_REPO/bin/python" "$XEE_REPO/scripts/workflow.py" doctor
