#!/usr/bin/env bash
# pool _pclx = train + COCONUT + PubChemLite + ChEBI + LIPID MAPS; V-in cost vs _pcl with the same ranker
cd "$(dirname "$0")/../.." || exit 1
PY=.venv/Scripts/python
export CASMI_POOL_TAG=_pclx CASMI_PUBCHEMLITE=C:/Users/MONSTER/kaggle-data/external/pubchemlite/PubChemLite_CCSbase_20260828.csv
export CASMI_EXTRA_SOURCES=C:/Users/MONSTER/kaggle-data/casmi26/artifacts
log() { echo "[$(date '+%H:%M:%S')] $*"; }
log "pool"; $PY -u scripts/build_pool.py 2>&1 | grep -vi warning | tail -n 4
log "fingerprints"; $PY -u scripts/build_fingerprints.py 2>&1 | tail -n 1
log "V-in table"; $PY -u scripts/build_rank_data.py --set vin --suffix _pclx > runs/pclx_vin.log 2>&1 || { log "vin FAILED"; exit 1; }
log "score"; $PY -u scripts/eval_pool_variant.py --ranker ranker-v13-clean-pcl 2>&1 | grep -E "^pool|^V-in"
unset CASMI_POOL_TAG
$PY scripts/compare_runs.py "$(ls -t runs/*_ranker-v13-clean-pcl_pcl_*.json | head -n 1)" "$(ls -t runs/*_ranker-v13-clean-pcl_pclx_*.json | head -n 1)"
date > runs/pclx.done
