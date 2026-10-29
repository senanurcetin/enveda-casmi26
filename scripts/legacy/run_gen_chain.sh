#!/usr/bin/env bash
# Class-3 hedge: rebuild ranker training + validation tables with generated candidates (pool _pcl),
# join D3 + D4, train the ranker, evaluate.
cd "$(dirname "$0")/../.." || exit 1
PY=.venv/Scripts/python
export CASMI_POOL_TAG=_pcl
A=C:/Users/MONSTER/kaggle-data/casmi26/artifacts
log() { echo "[$(date '+%H:%M:%S')] $*"; }
run() { name=$1; shift; [ -f "$A/$name" ] && { log "$name exists"; return 0; }; $PY -u scripts/build_rank_data.py "$@" > "runs/gen_${name%.parquet}.log" 2>&1 || { log "FAILED $name"; tail -n 5 "runs/gen_${name%.parquet}.log"; exit 1; }; log "$(tail -n 1 runs/gen_${name%.parquet}.log)"; }
log "training tables (2 shards in parallel)"
run rank_train_gen_0.parquet --set train --n 1500 --classes C1,C2,C3 --c1-mode own --generate --suffix _gen --shard 0 --n-shards 2 &
run rank_train_gen_1.parquet --set train --n 1500 --classes C1,C2,C3 --c1-mode own --generate --suffix _gen --shard 1 --n-shards 2 &
wait
log "enveda-180 Class-1 and V-real tables"
run rank_e180_gen_0.parquet --set e180 --n 1200 --classes C1 --generate --suffix _gen &
run rank_vreal_gen_0.parquet --set vreal --classes C1,C2,C3 --generate --suffix _gen &
wait
run rank_vin_gen_0.parquet --set vin --classes C1,C2,C3 --generate --suffix _gen
[ -f "$A/rank_train_gen_0.parquet" ] && [ -f "$A/rank_vin_gen_0.parquet" ] && date > runs/gen_tables.done
