#!/usr/bin/env bash
cd "$(dirname "$0")/.." || exit 1
until [ -f runs/tune.done ]; do sleep 60; done
free_gb10() { powershell -NoProfile -c "[int]((Get-CimInstance Win32_OperatingSystem).FreePhysicalMemory/1MB*10)" | tr -d '\r'; }
until [ "$(free_gb10)" -ge 30 ]; do sleep 60; done
.venv/Scripts/python -u scripts/probe_d4.py > runs/probe_d4.log 2>&1
date > runs/probe_d4.done
