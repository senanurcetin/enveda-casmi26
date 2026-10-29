#!/usr/bin/env bash
# memory-gated: ranker hyper-parameter sweep (cross-fitted on training queries)
cd "$(dirname "$0")/../.." || exit 1
free_gb10() { powershell -NoProfile -c "[int]((Get-CimInstance Win32_OperatingSystem).FreePhysicalMemory/1MB*10)" | tr -d '\r'; }
until [ "$(free_gb10)" -ge 45 ]; do sleep 60; done
echo "[$(date '+%H:%M:%S')] start ($(free_gb10)/10 GB free)" > runs/tune.log
.venv/Scripts/python -u scripts/tune_ranker.py >> runs/tune.log 2>&1
date > runs/tune.done
