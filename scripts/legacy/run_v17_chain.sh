#!/usr/bin/env bash
# v17: consensus + transform-type + NP-likeness features, depth-3 D4 block, tuned ranker params.
# One heavy process at a time, each gated on free memory (the machine is shared).
cd "$(dirname "$0")/../.." || exit 1
PY=.venv/Scripts/python
export CASMI_POOL_TAG=_pcl
A=C:/Users/MONSTER/kaggle-data/casmi26/artifacts
FLAGS="--generate --gen-consensus --np-score --suffix _gen3"
log() { echo "[$(date '+%H:%M:%S')] $*"; }
free_gb10() { powershell -NoProfile -c "[int]((Get-CimInstance Win32_OperatingSystem).FreePhysicalMemory/1MB*10)" | tr -d '\r'; }
gate() { until [ "$(free_gb10)" -ge 30 ]; do sleep 60; done; }
step() { name=$1; shift; [ -f "$A/$name" ] && grep -q "wrote $name" "runs/v17_${name%.parquet}.log" 2>/dev/null && { log "$name exists"; return 0; }
         gate; $PY -u scripts/build_rank_data.py "$@" > "runs/v17_${name%.parquet}.log" 2>&1 || { log "FAILED $name"; exit 1; }
         log "$(grep wrote runs/v17_${name%.parquet}.log)"; }
step rank_train_gen3_0.parquet --set train --n 3000 --classes C1,C2,C3 --c1-mode own $FLAGS --shard 0 --n-shards 2
step rank_train_gen3_1.parquet --set train --n 3000 --classes C1,C2,C3 --c1-mode own $FLAGS --shard 1 --n-shards 2
step rank_e180_gen3_0.parquet --set e180 --n 1200 --classes C1 $FLAGS
step rank_vreal_gen3_0.parquet --set vreal --classes C1,C2,C3 $FLAGS
step rank_vin_gen3_0.parquet --set vin --classes C1,C2,C3 $FLAGS
gate; log "merge + D3 + D4 + D4x"
$PY -u scripts/merge_train_gen.py --suffix _gen3 --n-train 3000 --d4x 2>&1 | grep -vi "warn\|model in" || exit 1
P='{"lambda_l2": 10.0, "min_data_in_leaf": 100}'
gate; log "rankers"
$PY -u scripts/train_ranker.py --train-files rank_train_gen3_d3_d4.parquet --no-membership --tag ranker-v17 --params "$P" --skip-vin 2>&1 | grep -E "features|best iter|saved"
$PY -u scripts/train_ranker.py --train-files rank_train_gen3_d3_d4.parquet --no-membership --tag ranker-v17-abl --params "$P" \
    --drop gen_sim_sum,gen_n_seeds,gen_n_paths,gen_t1,gen_t2,np_score,d4x_ --skip-vin 2>&1 | grep -E "features|saved"
log "evaluate"
for set in vreal vin; do
  $PY -u scripts/eval_gen.py --set $set --suffix _gen3 --scorers ranker-v16-gen2,ranker-v17-abl,ranker-v17 2>&1 | grep -vi "warn\|model in"
done
date > runs/v17.done
