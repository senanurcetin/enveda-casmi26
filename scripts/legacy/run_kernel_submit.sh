#!/usr/bin/env bash
# run_probe.sh NAME "MESSAGE" [VERSION]
# Waits for notebook senanuretin/casmi26-probe-NAME, submits it exactly once if it completed
# (success is read from the submissions list — the CLI is silent on code-competition success),
# then waits for the score. Never retries a submission.
cd "$(dirname "$0")/../.." || exit 1
NAME=$1; MSG=$2; VERSION=${3:-1}
C=enveda-CASMI26-molecule-id-mass-spectra
K=senanuretin/${KNAME:-casmi26-probe-$NAME}
log() { echo "[$(date '+%H:%M:%S')] $*"; }
count() { timeout 120 kaggle competitions submissions -c $C 2>&1 | grep -c "SubmissionStatus"; }

while true; do
  s=$(timeout 120 kaggle kernels status $K 2>&1 | grep -v Warning)
  case "$s" in *COMPLETE*) break;; *ERROR*|*CANCEL*) log "$NAME: notebook failed ($s) - not submitted"; exit 1;; esac
  sleep 60
done
log "$NAME: notebook complete"
# daily submission limits reset at 00:00 UTC: optionally hold the submission until NOT_BEFORE_UTC
if [ -n "$NOT_BEFORE_UTC" ]; then
  while [ "$(date -u '+%Y-%m-%dT%H:%M')" \< "$NOT_BEFORE_UTC" ]; do sleep 120; done
  log "$NAME: past $NOT_BEFORE_UTC UTC, submitting"
fi
before=$(count)
until [ -n "$before" ] && [ "$before" -gt 0 ]; do sleep 30; before=$(count); done   # a timed-out list reads as 0
timeout 300 kaggle competitions submit -c $C -k $K -v "$VERSION" -f submission.csv -m "$MSG" > /dev/null 2>&1
# the submissions list can lag or fail transiently: look a few times before giving up (never resubmit)
ok=0
for i in 1 2 3 4 5 6; do sleep 20; n=$(count); [ -n "$n" ] && [ "$n" -gt "$before" ] && { ok=1; break; }; done
[ "$ok" = 1 ] || { log "$NAME: submission not visible after 2 min - check the list by hand, do NOT resubmit"; exit 1; }
log "$NAME: submitted"
for i in $(seq 1 60); do
  row=$(timeout 120 kaggle competitions submissions -c $C 2>&1 | grep -v Warning | grep -F "$MSG" | head -n 1)
  echo "$row" | grep -q PENDING || break
  sleep 180
done
log "$NAME: $(echo "$row" | awk '{print $NF}')  <- public LB"
