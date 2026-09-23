#!/usr/bin/env bash
# Unattended second half of the modelling pipeline. Stages run strictly one after another: on a
# 16 GB laptop, the rank-data builder (~3.5 GB) and GPU training (~8.7 GB commit incl. the CUDA
# driver) together hit the commit limit. Each stage is skipped if its output already exists.
cd "$(dirname "$0")/.." || exit 1
PY=.venv/Scripts/python
A=/c/Users/MONSTER/kaggle-data/casmi26/artifacts
log() { echo "[$(date '+%H:%M:%S')] $*"; }

log "waiting for ranker tables"
while [ ! -f runs/rank_chain.done ]; do sleep 60; done

for f in 0 1; do
  if ! grep -q "^saved" runs/fp_fold$f.log 2>/dev/null; then
    log "training D3 fold $f"
    $PY -u scripts/train_fpmodel.py --train-folds $f --epochs 8 > runs/fp_fold$f.log 2>&1 || { log "D3 fold $f FAILED"; exit 1; }
  fi
  log "D3 fold $f: $(grep '^epoch' runs/fp_fold$f.log | tail -n 1)"
done

log "joining out-of-fold D3 features"
$PY -u scripts/add_d3_features.py > runs/add_d3.log 2>&1 || { log "add_d3 FAILED"; exit 1; }

BASE=$(ls -t runs/*_d1-d2-analog_*.json | head -n 1)
for variant in "--no-d3" ""; do
  log "training ranker $variant"
  $PY -u scripts/train_ranker.py $variant > runs/ranker${variant// /}.log 2>&1 || { log "ranker $variant FAILED"; exit 1; }
done
for tag in ranker-d1d2 ranker-d1d2d3; do
  RUN=$(ls -t runs/*_${tag}_*.json | head -n 1)
  log "=== $tag vs interim rule"
  $PY scripts/compare_runs.py "$BASE" "$RUN"
done
log "pipeline finished"
date > runs/pipeline.done
