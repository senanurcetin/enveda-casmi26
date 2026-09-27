#!/usr/bin/env bash
# which v17 addition hurt? one feature group removed at a time (tables unchanged, tuned params)
cd "$(dirname "$0")/.." || exit 1
PY=.venv/Scripts/python
free_gb10() { powershell -NoProfile -c "[int]((Get-CimInstance Win32_OperatingSystem).FreePhysicalMemory/1MB*10)" | tr -d '\r'; }
gate() { until [ "$(free_gb10)" -ge 35 ]; do sleep 60; done; }
log() { echo "[$(date '+%H:%M:%S')] $*"; }
P='{"lambda_l2": 10.0, "min_data_in_leaf": 100}'
T=rank_train_gen3_d3_d4.parquet
declare -A DROP=( [ranker-v17-nonp]=np_score [ranker-v17-nod4x]=d4x_ [ranker-v17-nocons]=gen_sim_sum,gen_n_seeds,gen_n_paths,gen_t1,gen_t2 )
for tag in ranker-v17-nonp ranker-v17-nod4x ranker-v17-nocons; do
  gate; log "train $tag (drop ${DROP[$tag]})"
  $PY -u scripts/train_ranker.py --train-files $T --no-membership --tag $tag --params "$P" --drop "${DROP[$tag]}" --skip-vin 2>&1 | grep -E "features|saved"
done
gate; log "evaluate"
for set in vreal vin; do
  $PY -u scripts/eval_gen.py --set $set --suffix _gen3 --scorers ranker-v17-nonp,ranker-v17-nod4x,ranker-v17-nocons 2>&1 | grep -vi "warn\|model in"
done
date > runs/v17_abl.done
