#!/usr/bin/env bash
# start the rest of the v16 chain once shard 1 has finished and >= 5 GB RAM is free (shared machine)
cd "$(dirname "$0")/../.." || exit 1
free_gb() { powershell -NoProfile -c "[int]((Get-CimInstance Win32_OperatingSystem).FreePhysicalMemory/1MB*10)" | tr -d '\r'; }
until grep -q "wrote rank_train_gen2_1" runs/gen2_rank_train_gen2_1.log && [ "$(free_gb)" -ge 50 ]; do sleep 60; done
echo "[$(date '+%H:%M:%S')] memory free ($(free_gb)/10 GB), starting" >> runs/gen2_rest.log
bash scripts/legacy/run_gen2_rest.sh >> runs/gen2_rest.log 2>&1
