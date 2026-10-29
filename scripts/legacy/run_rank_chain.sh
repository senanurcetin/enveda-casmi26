#!/usr/bin/env bash
# Builds all ranker tables one after another. Each process needs ~3.5 GB RAM; two in parallel
# pushed a 16 GB laptop into swap. Resumable: shards whose log says "wrote" are skipped.
cd "$(dirname "$0")/../.." || exit 1
PY=.venv/Scripts/python

still_running() {  # $1 = shard number of a train build that may already be running
  powershell -NoProfile -Command "if (Get-CimInstance Win32_Process -Filter \"Name='python.exe'\" | Where-Object { \$_.CommandLine -match ('build_rank_data.py --set train --n 1500 --shard ' + '$1') }) { exit 0 } else { exit 1 }"
}

while still_running 0; do sleep 30; done
grep -q "^wrote" runs/rank_train_1.log 2>/dev/null || $PY -u scripts/build_rank_data.py --set train --n 1500 --shard 1 --n-shards 2 > runs/rank_train_1.log 2>&1
grep -q "^wrote" runs/rank_vin.log 2>/dev/null || $PY -u scripts/build_rank_data.py --set vin > runs/rank_vin.log 2>&1
date > runs/rank_chain.done
