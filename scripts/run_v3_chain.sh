#!/usr/bin/env bash
cd "$(dirname "$0")/.." || exit 1
PY=.venv/Scripts/python
log() { echo "[$(date '+%H:%M:%S')] $*"; }
for s in 0 1; do
  if ! grep -q "^wrote" runs/rank_e180_$s.log 2>/dev/null; then
    log "enveda-180 Class-1 groups, shard $s"
    $PY -u scripts/build_rank_data.py --set e180 --n 1200 --classes C1 --shard $s --n-shards 2 > runs/rank_e180_$s.log 2>&1 || { log "e180 $s FAILED"; exit 1; }
  fi
  log "$(grep -E '^C1 skipped|^wrote' runs/rank_e180_$s.log | tr '\n' ' ')"
done
log "merge v3 + D3"
$PY -u scripts/merge_train_v3.py > runs/merge_v3.log 2>&1 || { log "merge FAILED"; tail -n 5 runs/merge_v3.log; exit 1; }
grep -E "^v3|C1 truth" runs/merge_v3.log
$PY -u scripts/train_ranker.py --train-files=rank_train_v3_d3.parquet --c2-weight 1 --tag ranker-v3 > runs/ranker_v3.log 2>&1 || { log "ranker FAILED"; exit 1; }
last() { ls -t runs/*_$1_*.json | head -n 1; }
for pair in "d1-d2-analog ranker-v3" "ranker-d1d2d3 ranker-v3" "ranker-v2-w1 ranker-v3"; do
  set -- $pair; log "=== $1 -> $2"; $PY scripts/compare_runs.py "$(last $1)" "$(last $2)"
done
date > runs/v3.done
