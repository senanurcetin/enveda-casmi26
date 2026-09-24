"""Build ranker training data (knockout-simulated queries) and the V-in evaluation table.

Training queries: NP-like structures (COCONUT members) from non-Enveda libraries, each turned
into a molecule-like query (≤ 5 spectra from its best-covered library), scored under C1 and C2
knockouts. Any structure sharing a heavy-atom skeleton with a V-in structure is excluded, which
is stricter than metric-key equivalence. Structures get a fold (0/1) from a hash of their skeleton
key, so tautomers always share a fold; the fold drives out-of-fold channel models later.

Usage:
  python scripts/build_rank_data.py --set train --n 1500 --shard 0 --n-shards 2
  python scripts/build_rank_data.py --set vin
"""
import argparse
import hashlib
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from casmi.config import ARTIFACT_DIR, POOL_FILE, POOL_FP_FILE, POOL_FP_KEYS_FILE  # noqa: E402
from casmi.channels.analog import AnalogIndex  # noqa: E402
from casmi.channels.library import LibraryIndex  # noqa: E402
from casmi.chem.canon import skeleton_key  # noqa: E402
from casmi.chem.keys import MetricKeys  # noqa: E402
from casmi.engine import BaselineEngine, EngineConfig, molecule_mass  # noqa: E402
from casmi.pool.pool import SOURCE_BITS, CandidatePool  # noqa: E402
from casmi.rank.features import add_context, add_priors  # noqa: E402
from casmi.spectra.store import SpectrumStore  # noqa: E402
from casmi.validate.knockout import Query, build_queries, equivalent_keys, knockout_for  # noqa: E402

TRAIN_LIBS = ["riken", "gnps", "massbank", "mona", "msdial", "spectraverse", "pluskal_ms2"]
TEST_MASS = (150.0, 1170.0)


def fold_of(skel: str) -> int:
    return int(hashlib.md5(skel.encode()).hexdigest(), 16) % 2


def training_queries(store, key_of_row, pool, vin_skeletons, n, seed=0):
    meta = store.meta
    ok = (meta.ingest_lib.isin(TRAIN_LIBS).to_numpy() & meta.test_adduct.to_numpy()
          & (meta.n_peaks.to_numpy() >= 5) & (np.abs(meta.ppm.to_numpy()) <= 20))
    rows = np.flatnonzero(ok)
    t = pd.DataFrame({"row": rows, "key": key_of_row[rows], "lib": meta.ingest_lib.astype(str).to_numpy()[rows]})
    src = pool.table.set_index("pool_key").sources
    t = t[t.key.map(src).fillna(0).astype(np.int64).to_numpy() & SOURCE_BITS["coconut"] > 0]
    mass = pool.table.set_index("pool_key").mass
    t = t[t.key.map(mass).between(*TEST_MASS)]
    best_lib = t.groupby(["key", "lib"]).size().reset_index(name="n").sort_values("n", ascending=False)
    best_lib = best_lib.drop_duplicates("key")
    t = t.merge(best_lib[["key", "lib"]], on=["key", "lib"])          # rows from each key's best library
    rows_of = {k: g.to_numpy() for k, g in t.groupby("key").row}
    lib_of = dict(zip(best_lib.key, best_lib.lib))
    rng = np.random.default_rng(seed)
    keys = best_lib.key.to_numpy()
    rng.shuffle(keys)
    out = []
    for k in keys:
        skel = skeleton_key(pool.smiles_of(k))
        if skel is None or skel in vin_skeletons:
            continue
        r = rows_of[k]
        if len(r) > 5:
            r = np.sort(rng.choice(r, 5, replace=False))
        spectra = []
        for i in r:
            mz, it = store.peaks(i)
            m = meta.iloc[i]
            spectra.append(dict(mz=mz, it=it, precursor_mz=float(m.precursor_mz), adduct=str(m.adduct),
                                ce_ev=float(m.ce_ev), mode=str(m.ionization_mode)))
        out.append((Query(qid=f"{lib_of[k]}:{k}", pool_key=k, spectra=spectra, spectrum_rows=r), fold_of(skel)))
        if len(out) == n:
            break
    return out


def e180_class1_queries(store, key_of_row, pool, vin_skeletons, n, seed=0):
    """Instrument-matched Class-1 queries: enveda-180 (timsTOF, the test instrument) structures that
    also have reference spectra in another library — exactly V-in's situation, where a timsTOF query
    must be matched to references acquired on other instruments."""
    meta = store.meta
    lib = meta.ingest_lib.astype(str).to_numpy()
    ok = (meta.test_adduct.to_numpy() & (meta.n_peaks.to_numpy() >= 5) & (np.abs(meta.ppm.to_numpy()) <= 20))
    keys = key_of_row
    e180 = pd.DataFrame({"row": np.flatnonzero(ok & (lib == "enveda-180"))})
    e180["key"] = keys[e180.row]
    other = set(keys[np.flatnonzero(ok & (lib != "enveda-180") & (lib != "enveda-np-examples"))])
    e180 = e180[e180.key.isin(other)]
    mass = pool.table.set_index("pool_key").mass
    e180 = e180[e180.key.map(mass).between(*TEST_MASS)]
    rows_of = {k: g.to_numpy() for k, g in e180.groupby("key").row}
    rng = np.random.default_rng(seed)
    order = np.array(sorted(rows_of))
    rng.shuffle(order)
    out = []
    for k in order:
        skel = skeleton_key(pool.smiles_of(k))
        if skel is None or skel in vin_skeletons:
            continue
        r = rows_of[k]
        if len(r) > 5:
            r = np.sort(rng.choice(r, 5, replace=False))
        spectra = []
        for i in r:
            mz, it = store.peaks(i)
            m = meta.iloc[i]
            spectra.append(dict(mz=mz, it=it, precursor_mz=float(m.precursor_mz), adduct=str(m.adduct),
                                ce_ev=float(m.ce_ev), mode=str(m.ionization_mode)))
        out.append((Query(qid=f"enveda-180:{k}", pool_key=k, spectra=spectra, spectrum_rows=r), fold_of(skel)))
        if len(out) == n:
            break
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--set", choices=["train", "vin", "e180"], required=True)
    ap.add_argument("--n", type=int, default=1500)
    ap.add_argument("--shard", type=int, default=0)
    ap.add_argument("--n-shards", type=int, default=1)
    ap.add_argument("--classes", default="C1,C2")
    ap.add_argument("--suffix", default="", help="output name suffix, e.g. _c1x")
    ap.add_argument("--c1-mode", choices=["crosslib", "own"], default="crosslib",
                    help="Class-1 knockout for train/e180 sets: remove the query library's spectra of the "
                         "structure (crosslib) or only the query's own spectra (own, as in v1/v13)")
    args = ap.parse_args()
    classes = args.classes.split(",")

    t0 = time.time()
    cfg = EngineConfig()
    store = SpectrumStore(ARTIFACT_DIR / "store")
    key_of_row = store.meta.train_ik14.to_numpy(dtype=object)
    pool = CandidatePool(pd.read_parquet(ARTIFACT_DIR / POOL_FILE))
    fp_keys = pd.read_parquet(ARTIFACT_DIR / POOL_FP_KEYS_FILE).pool_key.to_numpy()
    assert (fp_keys == pool.key).all()
    fps = np.load(ARTIFACT_DIR / POOL_FP_FILE)
    library = LibraryIndex(store, key_of_row, cfg.library)
    analog = AnalogIndex(store, key_of_row, cfg.analog, cache=ARTIFACT_DIR / "analog_index.npz")
    keys = MetricKeys(ARTIFACT_DIR / "metric_keys.parquet")                    # start from the shared cache
    keys.path = ARTIFACT_DIR / f"metric_keys_{args.set}_{args.shard}.parquet"   # but never write it concurrently
    engine = BaselineEngine(pool, library, cfg, analog, fps)

    vin = build_queries(store, key_of_row)
    if args.set == "vin":
        items = [(q, -1) for q in vin]
    else:
        vin_skel = {skeleton_key(pool.smiles_of(q.pool_key)) for q in vin}
        make = e180_class1_queries if args.set == "e180" else training_queries
        items = make(store, key_of_row, pool, vin_skel, args.n)
    items = items[args.shard::args.n_shards]
    print(f"{args.set} shard {args.shard}: {len(items)} queries, loaded in {time.time() - t0:.0f}s", flush=True)

    out_path = ARTIFACT_DIR / f"rank_{args.set}{args.suffix}_{args.shard}.parquet"
    lib_of_row = store.meta.ingest_lib.astype(str).to_numpy()
    skipped_c1 = 0
    tables = []
    for qi, (q, fold) in enumerate(items):
        mass = molecule_mass(q.spectra, cfg.mass_offset_ppm)
        truth, same = equivalent_keys(q.pool_key, mass, pool, keys, cfg.tol_ppm * 3)
        for cls in classes:
            ko = knockout_for(q, cls, library.rows_of, same)
            if cls == "C1" and args.set in ("train", "e180") and args.c1_mode == "crosslib":
                # Test Class-1 matches are always cross-library (timsTOF query vs another instrument's
                # reference). So remove every spectrum of the structure from the query's own library,
                # and skip the query when no other library has one — then it is not Class 1 at all.
                qlib = q.qid.split(":")[0]
                own = [r for k in same for r in library.rows_of.get(k, []) if lib_of_row[r] == qlib]
                ko.excluded_rows = np.union1d(ko.excluded_rows, np.asarray(own, np.int64))
                refs = [r for k in same for r in library.rows_of.get(k, []) if lib_of_row[r] != qlib]
                if not refs:
                    skipped_c1 += 1
                    continue
            f = engine.features([dict(s) for s in q.spectra], ko.excluded_rows, ko.excluded_keys)
            if f.empty:
                continue
            f = add_context(add_priors(f, pool, same if cls == "C2" else frozenset()))
            f["label"] = f.pool_key.isin(same).astype(np.int8)
            f["qid"], f["cls"], f["fold"] = q.qid, cls, fold
            tables.append(f)
        if (qi + 1) % 50 == 0:
            pd.concat(tables, ignore_index=True).to_parquet(out_path, index=False)
            keys.save()
            print(f"  {qi + 1}/{len(items)} ({time.time() - t0:.0f}s)", flush=True)
    pd.concat(tables, ignore_index=True).to_parquet(out_path, index=False)
    keys.save()
    print(f"C1 skipped (no cross-library reference): {skipped_c1}")
    print(f"wrote {out_path.name}: {sum(len(t) for t in tables):,} candidate rows ({time.time() - t0:.0f}s)")


if __name__ == "__main__":
    main()
