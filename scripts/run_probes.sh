#!/usr/bin/env bash
# Wait for each calibration-probe notebook, submit each successful one exactly once (never retried:
# the CLI is silent on code-competition success, so success is read from the submissions list),
# then wait until every new submission is scored.
cd "$(dirname "$0")/.." || exit 1
C=enveda-CASMI26-molecule-id-mass-spectra
log() { echo "[$(date '+%H:%M:%S')] $*"; }
count() { kaggle competitions submissions -c $C 2>&1 | grep -c "SubmissionStatus"; }
declare -A MSG=(
  [rule]="CALIBRATION: interim rule max(D1, D2), no learned ranker"
  [priors]="CALIBRATION: v1 LambdaRank D1+D2+D3 with membership priors (V-in C1 0.680 / C2 0.676)"
)
for name in rule priors; do
  K=senanuretin/casmi26-probe-$name
  while true; do
    s=$(kaggle kernels status $K 2>&1 | grep -v Warning)
    case "$s" in *COMPLETE*) break;; *ERROR*|*CANCEL*) log "$name: notebook failed ($s) — not submitted"; continue 2;; esac
    sleep 60
  done
  before=$(count)
  kaggle competitions submit -c $C -k $K -v 1 -f submission.csv -m "${MSG[$name]}" > /dev/null 2>&1
  sleep 15
  if [ "$(count)" -gt "$before" ]; then log "$name: submitted"; else log "$name: submission NOT registered"; fi
done
for i in $(seq 1 60); do
  rows=$(kaggle competitions submissions -c $C 2>&1 | grep -v Warning | grep CALIBRATION)
  echo "$rows" | grep -q PENDING || break
  sleep 180
done
log "results:"; echo "$rows" | awk '{print $(NF), $0}' | cut -c1-170
date > runs/probes.done
