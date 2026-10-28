"""Paths and runtime environment. Resolves the same code for local runs and Kaggle notebooks."""
import os
from pathlib import Path

ON_KAGGLE = Path("/kaggle/input").exists()
PROJECT_ROOT = Path(__file__).resolve().parents[2]


def _find_competition_data() -> Path:
    # Explicit research snapshots must take precedence over mounted competition
    # discovery, especially after creating a family-purged training corpus.
    if os.environ.get("CASMI_DATA"):
        return Path(os.environ["CASMI_DATA"])
    if ON_KAGGLE:
        hits = sorted(Path("/kaggle/input").rglob("test.parquet"))
        if hits:
            return hits[0].parent
    return Path(os.environ.get("CASMI_DATA", "C:/Users/MONSTER/kaggle-data/casmi26"))


DATA_DIR = _find_competition_data()
# Built artifacts (spectrum store, pool, models). Large, never committed.
# /tmp on Kaggle: scratch disk, not saved as notebook output (which would slow every commit)
ARTIFACT_DIR = Path(os.environ.get("CASMI_ARTIFACTS",
                                   "/tmp/casmi_artifacts" if ON_KAGGLE else str(DATA_DIR / "artifacts")))
OUTPUT_DIR = Path("/kaggle/working") if ON_KAGGLE else PROJECT_ROOT / "runs"
SEED = 42

# Candidate-pool variant: "" = train ∪ COCONUT (default), "_pcl" = + PubChemLite. All pool-derived
# artifacts carry the tag, so variants can live side by side and be compared on the same queries.
POOL_TAG = os.environ.get("CASMI_POOL_TAG", "")
POOL_FILE = f"pool{POOL_TAG}.parquet"
POOL_FP_FILE = f"pool_fp2048{POOL_TAG}.npy"
POOL_FP_KEYS_FILE = f"pool_fp_keys{POOL_TAG}.parquet"
