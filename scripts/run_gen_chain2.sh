#!/usr/bin/env bash
# after run_gen_chain.sh: D3 + D4, rankers with / without the hedge, evaluation on V-real and V-in
cd "$(dirname "$0")/.." || exit 1
PY=.venv/Scripts/python
export CASMI_POOL_TAG=_pcl
log() { echo "[$(date '+%H:%M:%S')] $*"; }
log "merge + D3 + D4"; $PY -u scripts/merge_train_gen.py 2>&1 | grep -vi warn || exit 1
log "rankers"
$PY -u scripts/train_ranker.py --train-files rank_train_gen_d3_d4.parquet --no-membership --tag ranker-v15-gen --skip-vin 2>&1 | grep -E "features|best iter|saved" 
$PY -u scripts/train_ranker.py --train-files rank_train_gen_d3_d4.parquet --no-membership --tag ranker-v15-gen-c3w2 --c3-weight 2 --skip-vin 2>&1 | grep -E "saved"
log "evaluate"
for set in vreal vin; do
  $PY -u scripts/eval_gen.py --set $set --scorers ranker-v13-clean-pcl:nogen,ranker-v14-d4:nogen,ranker-v15-gen:nogen,ranker-v15-gen,ranker-v15-gen-c3w2 2>&1 | grep -vi warn
done
date > runs/gen_eval.done
