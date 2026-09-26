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


def training_query_skeletons():
    """Skeletons of every structure any ranker was trained on (train and enveda-180 tables)."""
    smiles = structure_smiles()
    keys = set()
    for f in ARTIFACT_DIR.glob("rank_train*_[0-9].parquet"):
        keys |= {q.split(":", 1)[1] for q in pd.read_parquet(f, columns=["qid"]).qid.unique()}
    for f in ARTIFACT_DIR.glob("rank_e180*_[0-9].parquet"):
        keys |= {q.split(":", 1)[1] for q in pd.read_parquet(f, columns=["qid"]).qid.unique()}
    return {skeleton_key(smiles[k]) for k in keys if k in smiles}


_SMILES = {}


def structure_smiles():
    if not _SMILES:
        s = pd.read_parquet(ARTIFACT_DIR / "structures.parquet", columns=["pool_key", "smiles"])
        _SMILES.update(zip(s.pool_key, s.smiles))
    return _SMILES


def _vreal_build(store, key_of_row, vin_skeletons, n_per_stratum=300, seed=1):
    """V-real: Class-2-style queries that look like the hidden test rather than like V-in.

    Two strata of structures with spectra in the non-Enveda libraries, both never used as a ranker
    training query (skeleton-level exclusion, like V-in):
      pcl   in PubChemLite but not in COCONUT — answers reachable only through PubChem, the kind the
            leaderboard showed exists (pool + PubChemLite: +0.009) and V-in cannot contain
      coco  COCONUT members, like the training queries, but held out of every ranker
    Strata are defined on the PubChemLite pool so the query set is identical for every pool variant;
    in a pool without PubChemLite a `pcl` answer is simply absent (it scores 0, as on the test)."""
    meta = store.meta
    pcl_src = pd.read_parquet(ARTIFACT_DIR / "pool_pcl.parquet", columns=["pool_key", "sources", "mass"])
    src = dict(zip(pcl_src.pool_key, pcl_src.sources))
    mass = dict(zip(pcl_src.pool_key, pcl_src.mass))
    smiles = structure_smiles()
    exclude = vin_skeletons | training_query_skeletons()
    ok = (meta.ingest_lib.isin(TRAIN_LIBS).to_numpy() & meta.test_adduct.to_numpy()
          & (meta.n_peaks.to_numpy() >= 5) & (np.abs(meta.ppm.to_numpy()) <= 20))
    rows = np.flatnonzero(ok)
    t = pd.DataFrame({"row": rows, "key": key_of_row[rows], "lib": meta.ingest_lib.astype(str).to_numpy()[rows]})
    t = t[t.key.map(mass).between(*TEST_MASS)]
    s_ = t.key.map(src).fillna(0).astype(np.int64)
    t["stratum"] = np.where(s_ & SOURCE_BITS["coconut"], "coco", np.where(s_ & SOURCE_BITS["pubchem"], "pcl", ""))
    t = t[t.stratum != ""]
    best_lib = t.groupby(["key", "lib"]).size().reset_index(name="n").sort_values(["n", "key"], ascending=[False, True])
    best_lib = best_lib.drop_duplicates("key")
    t = t.merge(best_lib[["key", "lib"]], on=["key", "lib"])
    rows_of = {k: g.to_numpy() for k, g in t.groupby("key").row}
    lib_of = dict(zip(best_lib.key, best_lib.lib))
    stratum_of = dict(zip(t.key, t.stratum))
    rng = np.random.default_rng(seed)
    out = []
    for stratum in ("pcl", "coco"):
        keys = np.array(sorted(k for k in rows_of if stratum_of[k] == stratum))
        rng.shuffle(keys)
        n = 0
        for k in keys:
            skel = skeleton_key(smiles[k])
            if skel is None or skel in exclude:
                continue
            exclude = exclude | {skel}          # one query per skeleton
            r = rows_of[k]
            if len(r) > 5:
                r = np.sort(rng.choice(r, 5, replace=False))
            spectra = []
            for i in r:
                mz, it = store.peaks(i)
                m = meta.iloc[i]
                spectra.append(dict(mz=mz, it=it, precursor_mz=float(m.precursor_mz), adduct=str(m.adduct),
                                    ce_ev=float(m.ce_ev), mode=str(m.ionization_mode)))
            # the test's precursor masses are right; a query whose answer falls outside the mass
            # window (mis-annotated adduct or precursor) would only measure annotation noise
            if not abs(molecule_mass(spectra) - mass[k]) / mass[k] * 1e6 <= 10.0:
                continue
            out.append((Query(qid=f"{lib_of[k]}:{k}", pool_key=k, spectra=spectra, spectrum_rows=r), fold_of(skel)))
            n += 1
            if n == n_per_stratum:
                break
    return out, stratum_of


VREAL_FILE = "vreal_queries.parquet"


def vreal_queries(store, key_of_row, vin_skeletons, n_per_stratum=300, seed=1):
    """V-real, frozen on first use (artifacts/vreal_queries.parquet): later training tables must not
    change which structures V-real holds out, and must themselves exclude V-real's skeletons."""
    path = ARTIFACT_DIR / VREAL_FILE
    if not path.exists():
        items, stratum_of = _vreal_build(store, key_of_row, vin_skeletons, n_per_stratum, seed)
        pd.DataFrame({"qid": [q.qid for q, _ in items], "pool_key": [q.pool_key for q, _ in items],
                      "rows": [list(map(int, q.spectrum_rows)) for q, _ in items], "fold": [f for _, f in items],
                      "stratum": [stratum_of[q.pool_key] for q, _ in items]}).to_parquet(path, index=False)
    t = pd.read_parquet(path)
    meta, items = store.meta, []
    for qid, key, rows, fold in zip(t.qid, t.pool_key, t.rows, t.fold):
        r = np.asarray(rows, np.int64)
        spectra = []
        for i in r:
            mz, it = store.peaks(i)
            m = meta.iloc[i]
            spectra.append(dict(mz=mz, it=it, precursor_mz=float(m.precursor_mz), adduct=str(m.adduct),
                                ce_ev=float(m.ce_ev), mode=str(m.ionization_mode)))
        items.append((Query(qid=qid, pool_key=key, spectra=spectra, spectrum_rows=r), int(fold)))
    return items, dict(zip(t.pool_key, t.stratum))


def vreal_skeletons():
    path = ARTIFACT_DIR / VREAL_FILE
    if not path.exists():
        return set()
    smiles = structure_smiles()
    return {skeleton_key(smiles[k]) for k in pd.read_parquet(path, columns=["pool_key"]).pool_key}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--set", choices=["train", "vin", "e180", "vreal"], required=True)
    ap.add_argument("--n", type=int, default=1500)
    ap.add_argument("--shard", type=int, default=0)
    ap.add_argument("--n-shards", type=int, default=1)
    ap.add_argument("--classes", default="C1,C2")
    ap.add_argument("--skip", type=int, default=0, help="skip the shard's first K queries (resume after a crash)")
    ap.add_argument("--out-shard", type=int, default=None, help="output file index (default: --shard)")
    ap.add_argument("--generate", action="store_true", help="add Class-3 hedge candidates (biotransformed analogs)")
    ap.add_argument("--suffix", default="", help="output name suffix, e.g. _c1x")
    ap.add_argument("--c1-mode", choices=["crosslib", "own"], default="crosslib",
                    help="Class-1 knockout for train/e180 sets: remove the query library's spectra of the "
                         "structure (crosslib) or only the query's own spectra (own, as in v1/v13)")
    args = ap.parse_args()
    classes = args.classes.split(",")

    t0 = time.time()
    cfg = EngineConfig(generate=args.generate)
    store = SpectrumStore(ARTIFACT_DIR / "store")
    key_of_row = store.meta.train_ik14.to_numpy(dtype=object)
    pool = CandidatePool(pd.read_parquet(ARTIFACT_DIR / POOL_FILE))
    fp_keys = pd.read_parquet(ARTIFACT_DIR / POOL_FP_KEYS_FILE).pool_key.to_numpy()
    assert (fp_keys == pool.key).all()
    fps = np.load(ARTIFACT_DIR / POOL_FP_FILE)
    library = LibraryIndex(store, key_of_row, cfg.library)
    analog = AnalogIndex(store, key_of_row, cfg.analog, cache=ARTIFACT_DIR / "analog_index.npz")
    keys = MetricKeys(ARTIFACT_DIR / "metric_keys.parquet")                    # start from the shared cache
    keys.path = ARTIFACT_DIR / f"metric_keys_{args.set}_{args.shard if args.out_shard is None else args.out_shard}.parquet"   # but never write it concurrently
    engine = BaselineEngine(pool, library, cfg, analog, fps)

    vin = build_queries(store, key_of_row)
    stratum_of = {}
    if args.set == "vin":
        items = [(q, -1) for q in vin]
    elif args.set == "vreal":
        vin_skel = {skeleton_key(pool.smiles_of(q.pool_key)) for q in vin}
        items, stratum_of = vreal_queries(store, key_of_row, vin_skel)
    else:
        vin_skel = {skeleton_key(pool.smiles_of(q.pool_key)) for q in vin}
        make = e180_class1_queries if args.set == "e180" else training_queries
        items = make(store, key_of_row, pool, vin_skel | vreal_skeletons(), args.n)
    items = items[args.shard::args.n_shards][args.skip:]
    print(f"{args.set} shard {args.shard}: {len(items)} queries, loaded in {time.time() - t0:.0f}s", flush=True)

    out_path = ARTIFACT_DIR / f"rank_{args.set}{args.suffix}_{args.shard if args.out_shard is None else args.out_shard}.parquet"
    lib_of_row = store.meta.ingest_lib.astype(str).to_numpy()
    skipped_c1 = 0
    tables = []
    for qi, (q, fold) in enumerate(items):
        mass = molecule_mass(q.spectra, cfg.mass_offset_ppm)
        truth, same = equivalent_keys(q.pool_key, mass, pool, keys, cfg.tol_ppm * 3,
                                      truth_smiles=structure_smiles().get(q.pool_key))
        for cls in classes:
            ko = knockout_for(q, cls, library.rows_of, same | {q.pool_key})   # the query's own key, even if not in the pool
            if cls == "C2" and args.set == "vreal":
                # a spectrum-less test answer is a candidate only through a non-library source (COCONUT,
                # PubChemLite, ...): a V-real truth known only from the training libraries leaves the pool
                src = pool.table.sources.to_numpy()
                ko.excluded_keys = {k for k in same if src[pool.row_of[k]] & ~(SOURCE_BITS["train"] | SOURCE_BITS["enveda180_only"]) == 0}
            if cls == "C1" and args.set in ("train", "e180", "vreal") and args.c1_mode == "crosslib":
                # Test Class-1 matches are always cross-library (timsTOF query vs another instrument's
                # reference). So remove every spectrum of the structure from the query's own library,
                # and skip the query when no other library has one — then it is not Class 1 at all.
                qlib = q.qid.split(":")[0]
                own = [r for k in same | {q.pool_key} for r in library.rows_of.get(k, []) if lib_of_row[r] == qlib]
                ko.excluded_rows = np.union1d(ko.excluded_rows, np.asarray(own, np.int64))
                refs = [r for k in same | {q.pool_key} for r in library.rows_of.get(k, []) if lib_of_row[r] != qlib]
                if not refs:
                    skipped_c1 += 1
                    continue
            f = engine.features([dict(s) for s in q.spectra], ko.excluded_rows, ko.excluded_keys)
            if f.empty:
                continue
            f = add_context(add_priors(f, pool, same if cls == "C2" else frozenset()))
            f["label"] = f.pool_key.isin(same).astype(np.int8)
            if args.generate:
                g = (f.gen == 1).to_numpy()
                f["gen_smiles"] = ""
                f.loc[g, "gen_smiles"] = [pool.smiles_of(k) for k in f.pool_key[g]]
                tskel = skeleton_key(structure_smiles()[q.pool_key])
                f.loc[g, "label"] = [int(skeleton_key(s_) == tskel and keys.get(k, s_) == truth)
                                     for k, s_ in zip(f.pool_key[g], f.gen_smiles[g])]
            f["qid"], f["cls"], f["fold"] = q.qid, cls, fold
            if stratum_of:
                f["stratum"] = stratum_of[q.pool_key]
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
