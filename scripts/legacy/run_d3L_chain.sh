#!/usr/bin/env bash
# Large D3 end to end: train both folds -> recompute D3 on the _gen2 tables -> ranker (v16 recipe)
# -> V-real / V-in vs v16 -> if V-real ACCEPTs and V-in Class 1 is not significantly worse:
# upload, dry-run on Kaggle (no fallbacks), submit once.
cd "$(dirname "$0")/../.." || exit 1
PY=.venv/Scripts/python
export CASMI_POOL_TAG=_pcl
A=C:/Users/MONSTER/kaggle-data/casmi26/artifacts
U=C:/Users/MONSTER/kaggle-data/kaggle_upload
log() { echo "[$(date '+%H:%M:%S')] $*"; }
for f in 1 0; do
  [ -f "$A/models/fp_fold${f}_L.pt" ] && grep -q "epoch 6" runs/fp_L_$f.log 2>/dev/null && continue
  log "D3-large fold $f"
  $PY -u scripts/train_fpmodel.py --train-folds $f --epochs 6 --bs 384 --d-model 384 --n-layers 6 --d-ff 1536 --tag-suffix _L 2>&1 \
    | grep --line-buffered -E "^train|^epoch|Error|Traceback" > runs/fp_L_$f.log
  grep -q "epoch 6" runs/fp_L_$f.log || { log "fold $f failed"; exit 1; }
  tail -n 1 runs/fp_L_$f.log
done
log "D3 recompute"; $PY -u scripts/redo_d3.py --suffix _L 2>&1 | grep -vi warn || exit 1
log "ranker"; $PY -u scripts/train_ranker.py --train-files rank_train_gen2_L_d3_d4.parquet --no-membership --tag ranker-v18-d3L --skip-vin 2>&1 | grep -E "features|saved"
log "evaluate"
for set in vreal vin; do $PY -u scripts/eval_gen.py --set $set --suffix _gen2_L --scorers ranker-v18-d3L 2>&1 | grep -vi warn; done
ok=1
for set in vreal vin; do
  old=$(ls runs/*_${set}_gen2_ranker-v16-gen2_*.json | tail -1); new=$(ls runs/*_${set}_gen2_L_ranker-v18-d3L_*.json | tail -1)
  r=$($PY scripts/compare_runs.py $old $new); echo "== $set"; echo "$r" | grep -E "^C[123]|->"
  [ $set = vreal ] && { echo "$r" | grep -q "ACCEPT" || ok=0; }
  [ $set = vin ] && { echo "$r" | grep "^C1" | grep -q "worse" && ok=0; }
done
[ $ok = 1 ] || { log "not better than v16: no submission"; date > runs/d3L.done; exit 0; }
log "upload + dry run"
export KAGGLE_API_TOKEN=$(cat C:/Users/MONSTER/kaggle-data/.kt)
cp $A/models/fp_fold0_L.pt $A/models/fp_fold1_L.pt $A/models/ranker-v18-d3L_seed*.txt $U/casmi-models/
(cd $U/casmi-models && timeout 1800 kaggle datasets version -p . -m "D3-large + ranker-v18-d3L" > /dev/null 2>&1)
mkdir -p $U/probe_d3L
sed -e 's/probe-gen/probe-d3l/g' $U/probe_gen/kernel-metadata.json > $U/probe_d3L/kernel-metadata.json
$PY - <<'PYEOF'
import json
nb=json.load(open('C:/Users/MONSTER/kaggle-data/kaggle_upload/probe_gen/casmi26-probe-gen.ipynb',encoding='utf-8'))
nb['cells'][0]['source']=["# CASMI26 probe: v16 recipe with the large D3 model (d_model 384, 6 layers)\n"]
c=''.join(nb['cells'][1]['source']).replace("os.environ['CASMI_RANKER'] = 'ranker-v15-gen'","os.environ['CASMI_RANKER'] = 'ranker-v18-d3L'\nos.environ['CASMI_FP_SUFFIX'] = '_L'")
nb['cells'][1]['source']=c.splitlines(True)
json.dump(nb,open('C:/Users/MONSTER/kaggle-data/kaggle_upload/probe_d3L/casmi26-probe-d3l.ipynb','w',encoding='utf-8'),indent=1)
PYEOF
until timeout 120 kaggle datasets status senanuretin/casmi-models 2>&1 | grep -q ready; do sleep 30; done
(cd $U/probe_d3L && timeout 300 kaggle kernels push -p . 2>&1 | grep -v Warning)
K=senanuretin/casmi26-probe-d3l
while true; do s=$(timeout 120 kaggle kernels status $K 2>&1); case "$s" in *COMPLETE*) break;; *ERROR*|*CANCEL*) log "dry run failed"; exit 1;; esac; sleep 60; done
mkdir -p $U/probe_d3L_out; timeout 600 kaggle kernels output $K -p $U/probe_d3L_out --file-pattern ".*\.log$" > /dev/null 2>&1
$PY -c "
import json,glob
t=''.join(l['data'] for l in json.load(open(glob.glob('C:/Users/MONSTER/kaggle-data/kaggle_upload/probe_d3L_out/*.log')[0],encoding='utf-8')))
import re; m=re.search(r'(\d+) molecules, (\d+) fell back', t); print(m.group(0) if m else 'NO RESULT'); raise SystemExit(0 if m and m.group(2)=='0' else 1)" || { log "dry run had fallbacks: no submission"; exit 1; }
bash scripts/legacy/run_probe.sh d3l "v18-d3L: v16 recipe + large D3 model" 1
date > runs/d3L.done
