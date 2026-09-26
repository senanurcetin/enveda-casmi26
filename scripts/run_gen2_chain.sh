#!/usr/bin/env bash
# v16: extended transforms (42, 15 analog seeds) + 2x training queries (3000), tables suffix _gen2
cd "$(dirname "$0")/.." || exit 1
PY=.venv/Scripts/python
export CASMI_POOL_TAG=_pcl
A=C:/Users/MONSTER/kaggle-data/casmi26/artifacts
log() { echo "[$(date '+%H:%M:%S')] $*"; }
run() { name=$1; shift; [ -f "$A/$name" ] && { log "$name exists"; return 0; }; $PY -u scripts/build_rank_data.py "$@" > "runs/gen2_${name%.parquet}.log" 2>&1 || { log "FAILED $name"; tail -n 5 "runs/gen2_${name%.parquet}.log"; exit 1; }; log "$(tail -n 1 runs/gen2_${name%.parquet}.log)"; }
log "training tables (3000 queries, 2 shards)"
run rank_train_gen2_0.parquet --set train --n 3000 --classes C1,C2,C3 --c1-mode own --generate --suffix _gen2 --shard 0 --n-shards 2 &
run rank_train_gen2_1.parquet --set train --n 3000 --classes C1,C2,C3 --c1-mode own --generate --suffix _gen2 --shard 1 --n-shards 2 &
wait
[ -f "$A/rank_train_gen2_0.parquet" ] && [ -f "$A/rank_train_gen2_1.parquet" ] || { log "train tables missing"; exit 1; }
log "enveda-180 Class-1 and V-real"
run rank_e180_gen2_0.parquet --set e180 --n 1200 --classes C1 --generate --suffix _gen2 &
run rank_vreal_gen2_0.parquet --set vreal --classes C1,C2,C3 --generate --suffix _gen2 &
wait
run rank_vin_gen2_0.parquet --set vin --classes C1,C2,C3 --generate --suffix _gen2
log "merge + D3 + D4"
$PY -u scripts/merge_train_gen.py --suffix _gen2 --n-train 3000 2>&1 | grep -vi warn || exit 1
log "ranker"
$PY -u scripts/train_ranker.py --train-files rank_train_gen2_d3_d4.parquet --no-membership --tag ranker-v16-gen2 --skip-vin 2>&1 | grep -E "features|best iter|saved"
log "evaluate"
for set in vreal vin; do
  $PY -u scripts/eval_gen.py --set $set --suffix _gen2 --scorers ranker-v13-clean-pcl:nogen,ranker-v15-gen,ranker-v16-gen2 2>&1 | grep -vi warn
done
date > runs/gen2.done
