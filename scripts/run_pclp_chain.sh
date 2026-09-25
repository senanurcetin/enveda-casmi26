#!/usr/bin/env bash
# pool _pclp = train + COCONUT + PubChemLite + top-250k PubChem literature; V-in cost vs _pcl with the same ranker
cd "$(dirname "$0")/.." || exit 1
PY=.venv/Scripts/python
export CASMI_POOL_TAG=_pclp CASMI_PUBCHEMLITE=C:/Users/MONSTER/kaggle-data/external/pubchemlite/PubChemLite_CCSbase_20260828.csv
export CASMI_EXTRA_SOURCES=C:/Users/MONSTER/kaggle-data/kaggle_upload/casmi-pubchem-lit CASMI_EXTRA_NAMES= CASMI_PUBCHEM_MAX=250000
log() { echo "[$(date '+%H:%M:%S')] $*"; }
log "pool (built)"
log "fingerprints"; $PY -u scripts/build_fingerprints.py 2>&1 | tail -n 1
log "V-in table"; $PY -u scripts/build_rank_data.py --set vin --suffix _pclp > runs/pclp_vin.log 2>&1 || { log "vin FAILED"; exit 1; }
log "score"; $PY -u scripts/eval_pool_variant.py --ranker ranker-v13-clean-pcl 2>&1 | grep -E "^pool|^V-in"
unset CASMI_POOL_TAG
$PY scripts/compare_runs.py "$(ls -t runs/*_ranker-v13-clean-pcl_pcl_*.json | head -n 1)" "$(ls -t runs/*_ranker-v13-clean-pcl_pclp_*.json | head -n 1)"
date > runs/pclp.done
