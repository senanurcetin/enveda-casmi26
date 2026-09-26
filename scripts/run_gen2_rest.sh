#!/usr/bin/env bash
# rest of the v16 chain, one process at a time (the machine is shared; parallel builds ran out of memory)
cd "$(dirname "$0")/.." || exit 1
PY=.venv/Scripts/python
export CASMI_POOL_TAG=_pcl
A=C:/Users/MONSTER/kaggle-data/casmi26/artifacts
log() { echo "[$(date '+%H:%M:%S')] $*"; }
step() { name=$1; shift; [ -f "$A/$name" ] && { log "$name exists"; return 0; }; $PY -u scripts/build_rank_data.py "$@" > "runs/gen2_${name%.parquet}.log" 2>&1 || { log "FAILED $name"; exit 1; }; log "$(tail -n 1 runs/gen2_${name%.parquet}.log)"; }
shard1_running() { powershell -NoProfile -c "@(Get-CimInstance Win32_Process -Filter \"name='python.exe'\" | Where-Object { \$_.CommandLine -like '*--shard 1 --n-shards 2*' }).Count" | tr -d ''; }
while [ "$(shard1_running)" != "0" ]; do sleep 30; done
grep -q "wrote rank_train_gen2_1" runs/gen2_rank_train_gen2_1.log || { log "shard 1 did not finish"; exit 1; }
log "shard 1 done; tail of shard 0 (queries 1400-1499)"
step rank_train_gen2_3.parquet --set train --n 3000 --classes C1,C2,C3 --c1-mode own --generate --suffix _gen2 --shard 0 --n-shards 2 --skip 1400 --out-shard 3
step rank_e180_gen2_0.parquet --set e180 --n 1200 --classes C1 --generate --suffix _gen2
step rank_vreal_gen2_0.parquet --set vreal --classes C1,C2,C3 --generate --suffix _gen2
step rank_vin_gen2_0.parquet --set vin --classes C1,C2,C3 --generate --suffix _gen2
log "merge + D3 + D4"
$PY -u scripts/merge_train_gen.py --suffix _gen2 --n-train 3000 2>&1 | grep -vi warn || exit 1
log "ranker"
$PY -u scripts/train_ranker.py --train-files rank_train_gen2_d3_d4.parquet --no-membership --tag ranker-v16-gen2 --skip-vin 2>&1 | grep -E "features|best iter|saved"
log "evaluate"
for set in vreal vin; do
  $PY -u scripts/eval_gen.py --set $set --suffix _gen2 --scorers ranker-v13-clean-pcl:nogen,ranker-v15-gen,ranker-v16-gen2 2>&1 | grep -vi warn
done
date > runs/gen2.done
