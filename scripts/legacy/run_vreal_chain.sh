#!/usr/bin/env bash
# V-real tables for every pool the leaderboard has scored, then the calibration check.
cd "$(dirname "$0")/../.." || exit 1
PY=.venv/Scripts/python
log() { echo "[$(date '+%H:%M:%S')] $*"; }
for tag in "" _pcl _pclx _pclp; do
  [ -f "C:/Users/MONSTER/kaggle-data/casmi26/artifacts/rank_vreal${tag}_0.parquet" ] && { log "vreal$tag exists"; continue; }
  log "V-real table, pool '$tag'"
  CASMI_POOL_TAG=$tag $PY -u scripts/build_rank_data.py --set vreal --suffix "$tag" > "runs/vreal${tag}.log" 2>&1 || { log "FAILED $tag"; tail -n 5 "runs/vreal${tag}.log"; exit 1; }
  tail -n 2 "runs/vreal${tag}.log"
done
log "evaluate"
$PY -u scripts/eval_vreal.py 2>&1 | grep -v -i warning
date > runs/vreal.done
