#!/usr/bin/env bash
# After the modelling chain: end-to-end submission run on the public test file with the ranker
# that won validation — checks format, fallbacks and per-stage runtime before anything goes to Kaggle.
cd "$(dirname "$0")/.." || exit 1
PY=.venv/Scripts/python
while [ ! -f runs/pipeline.done ]; do sleep 60; done
# pick the ranker variant with the higher V-in C2 MRR (C1 is ~saturated; C2 is where they differ)
BEST=$($PY - <<'EOF'
import json, glob
best = max((json.load(open(f)) for f in glob.glob("runs/*_ranker-*.json")), key=lambda r: r["classes"]["C2"]["mrr"])
print(best["tag"])
EOF
)
echo "best ranker by V-in C2 MRR: $BEST"
$PY -u scripts/make_submission.py --ranker "$BEST" --sanity
date > runs/submission_check.done
