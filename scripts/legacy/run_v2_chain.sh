#!/usr/bin/env bash
# v2: rebuild Class-1 training groups with cross-library knockout, retrain and re-judge rankers.
cd "$(dirname "$0")/../.." || exit 1
PY=.venv/Scripts/python
log() { echo "[$(date '+%H:%M:%S')] $*"; }
for s in 0 1; do
  if ! grep -q "^wrote" runs/rank_c1x_$s.log 2>/dev/null; then
    log "C1 cross-library groups, shard $s"
    $PY -u scripts/build_rank_data.py --set train --n 1500 --classes C1 --suffix _c1x --shard $s --n-shards 2 > runs/rank_c1x_$s.log 2>&1 || { log "c1x $s FAILED"; exit 1; }
  fi
  log "$(grep -E '^C1 skipped' runs/rank_c1x_$s.log)"
done
log "merge v2 + D3"
$PY -u scripts/merge_train_v2.py > runs/merge_v2.log 2>&1 || { log "merge FAILED"; exit 1; }
cat runs/merge_v2.log | grep -E "^v2|C1 truth"
T=--train-files=rank_train_v2_d3.parquet
$PY -u scripts/train_ranker.py $T --c2-weight 1 --tag ranker-v2-w1 > runs/ranker_v2_w1.log 2>&1 || log "w1 FAILED"
$PY -u scripts/train_ranker.py $T --c2-weight 2 --tag ranker-v2-w2 > runs/ranker_v2_w2.log 2>&1 || log "w2 FAILED"
$PY -u scripts/train_ranker.py $T --c2-weight 1 --drop-d3 --tag ranker-v2-w1-nod3 > runs/ranker_v2_w1_nod3.log 2>&1 || log "nod3 FAILED"
last() { ls -t runs/*_$1_*.json | head -n 1; }
for pair in "d1-d2-analog ranker-v2-w1" "ranker-d1d2d3 ranker-v2-w1" "ranker-v2-w1 ranker-v2-w2" "ranker-v2-w1-nod3 ranker-v2-w1"; do
  set -- $pair
  log "=== $1 -> $2"
  $PY scripts/compare_runs.py "$(last $1)" "$(last $2)"
done
date > runs/v2.done
