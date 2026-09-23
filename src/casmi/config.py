"""Paths and runtime environment. Resolves the same code for local runs and Kaggle notebooks."""
import os
from pathlib import Path

ON_KAGGLE = Path("/kaggle/input").exists()
PROJECT_ROOT = Path(__file__).resolve().parents[2]


def _find_competition_data() -> Path:
    if ON_KAGGLE:
        hits = sorted(Path("/kaggle/input").rglob("test.parquet"))
        if hits:
            return hits[0].parent
    return Path(os.environ.get("CASMI_DATA", "C:/Users/MONSTER/kaggle-data/casmi26"))


DATA_DIR = _find_competition_data()
# Built artifacts (spectrum store, pool, models). Large, never committed.
ARTIFACT_DIR = Path(os.environ.get("CASMI_ARTIFACTS",
                                   "/kaggle/working/artifacts" if ON_KAGGLE else str(DATA_DIR / "artifacts")))
OUTPUT_DIR = Path("/kaggle/working") if ON_KAGGLE else PROJECT_ROOT / "runs"
SEED = 42
