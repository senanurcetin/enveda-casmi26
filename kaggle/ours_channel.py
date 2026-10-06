"""Our retrieve-and-rank pipeline as a library call, for use inside another Kaggle notebook (internet off).

`run_ours(test_parquet, out_csv)` does what kaggle/submission_driver.py does, with two differences:
the test file is a parameter (the host notebook may score a smoke subset) and the finished list is moved to
`out_csv` so it cannot collide with the host notebook's own submission.csv.

It never raises into the host: any failure returns None and the host keeps its own ranking.
Inputs (Kaggle datasets): casmi-code, casmi-models, casmi-coconut, casmi-pubchemlite.
"""
import glob
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path


def _find(pattern):
    hits = sorted(glob.glob(pattern, recursive=True))
    if not hits:
        raise FileNotFoundError(pattern)
    return Path(hits[0])


def run_ours(test_parquet, out_csv, ranker="ranker-v16-gen2", log=print):
    t0 = time.time()
    say = lambda m: log(f"[ours {time.time() - t0:6.0f}s] {m}")
    out_csv = Path(out_csv)
    kaggle_sub = Path("/kaggle/working/submission.csv")   # where make_submission writes on Kaggle
    parked = None
    try:
        code_src = _find("/kaggle/input/**/casmi_code_root.txt").parent
        models_src = _find("/kaggle/input/**/fp_fold0.pt").parent
        coconut = _find("/kaggle/input/**/coconut_csv_lite-09-2026.*")
        pcl = _find("/kaggle/input/**/PubChemLite_CCSbase_*.csv")
        work, artifacts = Path("/kaggle/working/casmi_ours"), Path("/tmp/casmi_artifacts")
        shutil.copytree(code_src, work, ignore=shutil.ignore_patterns("wheels", "*.pyc", "__pycache__"), dirs_exist_ok=True)
        try:
            import rdkit  # noqa: F401
        except ImportError:
            subprocess.run([sys.executable, "-m", "pip", "install", "--no-index", "--no-deps",
                            f"--find-links={code_src / 'wheels'}", "rdkit==2026.3.3"], check=True)
        import rdkit
        assert rdkit.__version__ == "2026.03.3", f"RDKit {rdkit.__version__}: metric keys would not match"
        (artifacts / "models").mkdir(parents=True, exist_ok=True)
        for f in models_src.iterdir():
            shutil.copy(f, artifacts / "models" / f.name)
        env = {**os.environ, "CASMI_ARTIFACTS": str(artifacts), "CASMI_COCONUT": str(coconut),
               "CASMI_PUBCHEMLITE": str(pcl), "CASMI_POOL_TAG": "_pcl", "PYTHONUNBUFFERED": "1"}
        for step in ["build_structures.py", "build_store.py", "build_pool.py", "build_fingerprints.py"]:
            say(f"-> {step}")
            subprocess.run([sys.executable, str(work / "scripts" / step)], env=env, check=True)
        if kaggle_sub.exists():                       # never clobber the host's file
            parked = kaggle_sub.with_name("submission.host_parked.csv")
            kaggle_sub.rename(parked)
        say(f"-> make_submission.py (ranker {ranker}) on {test_parquet}")
        subprocess.run([sys.executable, str(work / "scripts" / "make_submission.py"),
                        "--ranker", ranker, "--test", str(test_parquet)], env=env, check=True)
        out_csv.parent.mkdir(parents=True, exist_ok=True)
        shutil.move(str(kaggle_sub), str(out_csv))
        say(f"done -> {out_csv}")
        return out_csv
    except Exception as e:  # noqa: BLE001 - the host must survive anything
        say(f"FAILED -> host ranking kept: {e!r}")
        return None
    finally:
        if parked is not None and parked.exists():
            if kaggle_sub.exists():
                kaggle_sub.unlink()
            parked.rename(kaggle_sub)
