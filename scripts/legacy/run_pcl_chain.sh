#!/usr/bin/env bash
# PubChemLite pool variant: fingerprints -> V-in table -> cost on V-in, paired vs the default pool.
cd "$(dirname "$0")/../.." || exit 1
PY=.venv/Scripts/python
export CASMI_POOL_TAG=_pcl
log() { echo "[$(date '+%H:%M:%S')] $*"; }
log "fingerprints"; $PY -u scripts/build_fingerprints.py > runs/pcl_fp.log 2>&1 || { log "fp FAILED"; exit 1; }; tail -n 1 runs/pcl_fp.log
log "V-in table on the _pcl pool"; $PY -u scripts/build_rank_data.py --set vin --suffix _pcl > runs/pcl_vin.log 2>&1 || { log "vin FAILED"; exit 1; }; grep "^wrote" runs/pcl_vin.log
log "score"; $PY -u scripts/eval_pool_variant.py --ranker ranker-v13-clean > runs/pcl_eval.log 2>&1 || { log "eval FAILED"; tail -n 5 runs/pcl_eval.log; exit 1; }
grep -E "^pool|^V-in" runs/pcl_eval.log
unset CASMI_POOL_TAG
$PY scripts/compare_runs.py "$(ls -t runs/*_ranker-v13-clean_*.json | grep -v _pcl_ | head -n 1)" "$(ls -t runs/*_ranker-v13-clean_pcl_*.json | head -n 1)"
date > runs/pcl.done
