"""PubChem literature subset for the candidate pool -> artifacts/extra_pubchem.parquet.

A compound that appears in at least one PubMed article (CID-PMID) is a plausible "known but
spectrum-less" answer (Class 2); the full 120 M-compound PubChem is not (mostly combinatorial
and vendor records). Rows carry n_pmid so the pool build can cap the subset by literature depth.

  CID-PMID.gz, CID-SMILES.gz  from https://ftp.ncbi.nlm.nih.gov/pubchem/Compound/Extras/

Normalisation matches build_extra_sources.py: single component, net charge 0, mass 100–1300,
RDKit standard InChIKey block 1.
"""
import gzip
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from casmi.config import ARTIFACT_DIR  # noqa: E402

SRC = Path("C:/Users/MONSTER/kaggle-data/external/pubchem")
MASS_RANGE = (100.0, 1300.0)
ORGANIC = {"C", "H", "N", "O", "P", "S", "F", "Cl", "Br", "I", "Si", "Se", "As", "B"}


def _rows(chunk):
    from rdkit import Chem, RDLogger
    from rdkit.Chem.Descriptors import ExactMolWt
    RDLogger.DisableLog("rdApp.*")
    out = []
    for cid, n, s in chunk:
        m = Chem.MolFromSmiles(s)
        if m is None or Chem.GetFormalCharge(m) != 0:
            continue
        if any(a.GetSymbol() not in ORGANIC for a in m.GetAtoms()):
            continue
        mass = ExactMolWt(m)
        if not (MASS_RANGE[0] <= mass <= MASS_RANGE[1]):
            continue
        key = Chem.MolToInchiKey(m)
        if key:
            out.append((key[:14], Chem.MolToSmiles(m), mass, cid, n))
    return out


def pmid_counts():
    t0 = time.time()
    links = pd.read_csv(SRC / "CID-PMID.gz", sep="\t", header=None, usecols=[0, 1], names=["cid", "pmid"],
                        dtype=np.int64, engine="c")
    counts = links.drop_duplicates().cid.value_counts()   # distinct articles per compound
    print(f"CID-PMID: {len(counts):,} compounds with >=1 PubMed link ({time.time() - t0:.0f}s)", flush=True)
    return counts


def main():
    t0 = time.time()
    counts = pmid_counts()
    keep = set(counts.index)
    todo = []
    with gzip.open(SRC / "CID-SMILES.gz", "rt") as f:
        for line in f:
            cid, _, smi = line.partition("\t")
            cid = int(cid)
            if cid in keep:
                smi = smi.rstrip()
                if "." not in smi:
                    todo.append((cid, int(counts[cid]), smi))
    print(f"SMILES for {len(todo):,} single-component literature compounds ({time.time() - t0:.0f}s)", flush=True)
    chunks = [todo[i:i + 5000] for i in range(0, len(todo), 5000)]
    with ProcessPoolExecutor(max_workers=6) as ex:
        rows = [r for c in ex.map(_rows, chunks, chunksize=4) for r in c]
    df = pd.DataFrame(rows, columns=["pool_key", "smiles", "mass", "cid", "n_pmid"])
    # one row per key: the best-documented record
    df = df.sort_values(["n_pmid", "cid"], ascending=[False, True]).drop_duplicates("pool_key")
    df.to_parquet(ARTIFACT_DIR / "extra_pubchem.parquet", index=False)
    q = df.n_pmid.quantile([.5, .9, .99]).tolist()
    print(f"extra_pubchem: {len(df):,} structures, n_pmid median/p90/p99 {q} ({time.time() - t0:.0f}s)")


if __name__ == "__main__":
    main()
