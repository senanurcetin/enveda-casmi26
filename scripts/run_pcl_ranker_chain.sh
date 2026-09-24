#!/usr/bin/env bash
# Retrain the v13-clean ranker on candidate lists from the PubChemLite pool, then score V-in (_pcl).
cd "$(dirname "$0")/.." || exit 1
PY=.venv/Scripts/python
export CASMI_POOL_TAG=_pcl
log() { echo "[$(date '+%H:%M:%S')] $*"; }
for s in 0 1; do
  grep -q "^wrote" runs/pcl_train_$s.log 2>/dev/null || { log "train shard $s"; $PY -u scripts/build_rank_data.py --set train --n 1500 --c1-mode own --suffix _pcl --shard $s --n-shards 2 > runs/pcl_train_$s.log 2>&1 || { log "train $s FAILED"; exit 1; }; }
  grep -q "^wrote" runs/pcl_e180_$s.log 2>/dev/null || { log "e180 shard $s"; $PY -u scripts/build_rank_data.py --set e180 --n 1200 --classes C1 --suffix _pcl --shard $s --n-shards 2 > runs/pcl_e180_$s.log 2>&1 || { log "e180 $s FAILED"; exit 1; }; }
done
log "merge + D3"; $PY -u scripts/merge_train_pcl.py > runs/pcl_merge.log 2>&1 || { log "merge FAILED"; tail -n 5 runs/pcl_merge.log; exit 1; }; tail -n 1 runs/pcl_merge.log
log "ranker"; $PY -u scripts/train_ranker.py --train-files=rank_train_v13_pcl_d3.parquet --no-membership --tag ranker-v13-clean-pcl > runs/pcl_ranker.log 2>&1 || { log "ranker FAILED"; exit 1; }
log "score on V-in (_pcl)"; $PY -u scripts/eval_pool_variant.py --ranker ranker-v13-clean-pcl > runs/pcl_eval2.log 2>&1 || { log "eval FAILED"; exit 1; }
grep -E "^V-in" runs/pcl_eval2.log
$PY scripts/compare_runs.py "$(ls -t runs/*_ranker-v13-clean_pcl_*.json | head -n 1)" "$(ls -t runs/*_ranker-v13-clean-pcl_pcl_*.json | head -n 1)"
date > runs/pcl_ranker.done
