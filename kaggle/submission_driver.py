"""Enveda CASMI 2026 — submission driver (runs with internet off).

1. installs the pinned RDKit (2026.3.3, the metric's version) from the bundled wheel,
2. rebuilds every train-derived asset inside the notebook — competition data may not be
   redistributed, so the spectrum store, structure table, candidate pool and fingerprints are
   built here from train.parquet (+ the CC0 COCONUT file),
3. predicts test.parquet with the D1/D2/D3 channels and the LambdaRank ranker -> submission.csv.

Inputs (Kaggle datasets): casmi-code (src, scripts, wheels), casmi-models, casmi-coconut.
Locally, pass paths through CASMI_CODE / CASMI_MODELS / CASMI_COCONUT / CASMI_ARTIFACTS to
rehearse the exact same sequence.
"""
import glob
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

# which ranker to use; "none" = the interim rule (max of direct and propagated evidence)
RANKER = os.environ.get("CASMI_RANKER", "ranker-v13-clean")
T0 = time.time()


def log(msg):
    print(f"[{time.time() - T0:7.0f}s] {msg}", flush=True)


def find(pattern):
    hits = sorted(glob.glob(pattern, recursive=True))
    if not hits:
        raise FileNotFoundError(pattern)
    return Path(hits[0])


ON_KAGGLE = Path("/kaggle/input").exists()
if ON_KAGGLE:
    code_src = find("/kaggle/input/**/casmi_code_root.txt").parent
    models_src = find("/kaggle/input/**/fp_fold0.pt").parent
    coconut = find("/kaggle/input/**/coconut_csv_lite-09-2026.*")   # Kaggle extracts the zip to .csv
    work, artifacts = Path("/kaggle/working/casmi"), Path("/tmp/casmi_artifacts")
else:
    code_src, models_src = Path(os.environ["CASMI_CODE"]), Path(os.environ["CASMI_MODELS"])
    coconut, artifacts = Path(os.environ["CASMI_COCONUT"]), Path(os.environ["CASMI_ARTIFACTS"])
    work = artifacts / "code"
log(f"code {code_src} | models {models_src} | coconut {coconut.name}")

# numba caches compiled functions next to their source: the input mount is read-only
shutil.copytree(code_src, work, ignore=shutil.ignore_patterns("wheels", "*.pyc", "__pycache__"), dirs_exist_ok=True)

try:
    import rdkit  # noqa: F401
except ImportError:
    subprocess.run([sys.executable, "-m", "pip", "install", "--no-index", "--no-deps",
                    f"--find-links={code_src / 'wheels'}", "rdkit==2026.3.3"], check=True)
import rdkit  # noqa: E402
assert rdkit.__version__ == "2026.03.3", f"RDKit {rdkit.__version__}: metric keys would not match"
log(f"RDKit {rdkit.__version__}")

(artifacts / "models").mkdir(parents=True, exist_ok=True)
for f in models_src.iterdir():
    shutil.copy(f, artifacts / "models" / f.name)

env = {**os.environ, "CASMI_ARTIFACTS": str(artifacts), "CASMI_COCONUT": str(coconut), "PYTHONUNBUFFERED": "1"}
# optional pool expansion with PubChemLite for Exposomics (CC BY 4.0, Schymanski et al., Zenodo 22143470)
if os.environ.get("CASMI_USE_PUBCHEMLITE") == "1":
    pcl = find("/kaggle/input/**/PubChemLite_CCSbase_*.csv") if ON_KAGGLE else Path(os.environ["CASMI_PUBCHEMLITE"])
    env.update({"CASMI_PUBCHEMLITE": str(pcl), "CASMI_POOL_TAG": "_pcl"})
    log(f"pool expansion: {pcl.name}")
for step in ["build_structures.py", "build_store.py", "build_pool.py", "build_fingerprints.py"]:
    log(f"-> {step}")
    subprocess.run([sys.executable, str(work / "scripts" / step)], env=env, check=True)

log(f"-> make_submission.py (ranker {RANKER})")
args = [] if RANKER == "none" else ["--ranker", RANKER]
subprocess.run([sys.executable, str(work / "scripts" / "make_submission.py"), *args], env=env, check=True)
log("done")
