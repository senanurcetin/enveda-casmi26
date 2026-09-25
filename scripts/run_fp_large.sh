#!/usr/bin/env bash
# D3 capacity experiment: d_model 384, 6 layers, 12 epochs, both folds (models/fp_fold{0,1}_L.pt)
cd "$(dirname "$0")/.." || exit 1
for f in 0 1; do
  .venv/Scripts/python -u scripts/train_fpmodel.py --train-folds $f --epochs 12 --d-model 384 --n-layers 6 --d-ff 1536 --tag-suffix _L 2>&1 | grep -E "^train|^epoch|Error|Traceback" > runs/fp_large_$f.log
done
date > runs/fp_large.done
