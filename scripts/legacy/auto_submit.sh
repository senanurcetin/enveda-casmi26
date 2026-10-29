#!/usr/bin/env bash
# Submits casmi26-submission (version given as $1) as soon as Kaggle's identity verification is
# approved, then waits for the score. An unverified attempt returns 403 and uses no daily submission.
# Exactly one submission is made; the loop gives up after 3 days.
# NOTE: for code competitions the CLI prints nothing on success, so success is detected from the
# submissions list, never from the CLI message (relying on the message caused a duplicate submission).
cd "$(dirname "$0")/../.." || exit 1
VERSION=${1:-2}
C=enveda-CASMI26-molecule-id-mass-spectra
MSG="v13-clean LambdaRank over D1 direct + D2 analog + D3 spectrum->FP; V-in C1 0.839 / C2 0.662"
log() { echo "[$(date '+%Y-%m-%d %H:%M')] $*"; }

count() { kaggle competitions submissions -c $C 2>&1 | grep -c "SubmissionStatus"; }
before=$(count)
for attempt in $(seq 1 144); do
  out=$(kaggle competitions submit -c $C -k senanuretin/casmi26-submission -v "$VERSION" -f submission.csv -m "$MSG" 2>&1 | grep -v Warning)
  if [ "$(count)" -gt "$before" ]; then log "SUBMITTED (attempt $attempt)"; break; fi
  if echo "$out" | grep -q "403"; then log "not verified yet (attempt $attempt)"
  else log "unexpected response, stopping: $out"; exit 1; fi
  [ "$attempt" = 144 ] && { log "gave up after 3 days"; exit 1; }
  sleep 1800
done

# the notebook re-runs on the hidden test (~20-40 min), then it is scored
for i in $(seq 1 60); do
  row=$(kaggle competitions submissions -c $C 2>&1 | grep -v Warning | sed -n '3p')
  if echo "$row" | grep -qE "COMPLETE|ERROR"; then log "RESULT: $row"; echo "$row" > runs/first_submission.result; exit 0; fi
  sleep 300
done
log "no result after 5 h: $row"
