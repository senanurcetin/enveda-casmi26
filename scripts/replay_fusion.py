"""Replay the final fusion of the v4n fusion notebook offline, on the stage dumps of a validation run.

One validation commit run on Kaggle (notebooks/build_fusion_notebook.py --validation) writes stage_cache/*.json:
the base lists, PubChem lists, engine-2 lists, ICEBERG/GLACIER/frag scores, our lists and val_labels.csv. This
script re-runs the notebook's own fusion_core.build_submission on them for a grid of CFG overrides and scores
each with the metric key (InChIKey14 of the RDKit 2026.03.3 canonical tautomer), so fusion weights can be chosen
without spending submissions.

usage: python scripts/replay_fusion.py STAGE_DIR --notebook NB.ipynb --fuse-dir DIR [--pop-dir DIR]
           --grid OURS_BETA=0,0.2,0.4,0.8 [--grid OURS_MAX_LIB=0.5,0.7,0.9]

--notebook is the built fusion notebook (fusion_core is read from its embedded sources, as data);
--fuse-dir holds fuse.py and gl_fuse.py from the public ICEBERG / GLACIER datasets (only their pure re-rank
functions are used). The first grid point is the reference for the paired bootstrap.
"""
import argparse
import ast
import importlib.util
import itertools
import json
import os
import sys
import types

import numpy as np
import pandas as pd


def load_fusion_core(nb_path):
    nb = json.load(open(nb_path))
    src = "".join(nb["cells"][2]["source"])
    node = next(n for n in ast.parse(src).body
                if isinstance(n, ast.Assign) and getattr(n.targets[0], "id", "") == "EMBED")
    code = ast.literal_eval(node.value)["fusion_core.py"]
    mod = types.ModuleType("fusion_core")
    exec(compile(code, "fusion_core.py", "exec"), mod.__dict__)
    return mod


def load_module(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def make_keys():
    from rdkit import Chem, RDLogger
    from rdkit.Chem.MolStandardize import rdMolStandardize
    from rdkit.Chem.rdMolDescriptors import CalcMolFormula
    RDLogger.DisableLog("rdApp.*")
    taut = rdMolStandardize.TautomerEnumerator()
    cache_k, cache_f = {}, {}

    def score_key(s):
        if s not in cache_k:
            try:
                cache_k[s] = Chem.MolToInchiKey(taut.Canonicalize(Chem.MolFromSmiles(s)))[:14]
            except Exception:
                cache_k[s] = None
        return cache_k[s]

    def formula(s):
        if s not in cache_f:
            try:
                cache_f[s] = CalcMolFormula(Chem.MolFromSmiles(s))
            except Exception:
                cache_f[s] = s
        return cache_f[s]

    return score_key, formula


def reciprocal_ranks(sub, labels, score_key):
    rr = {}
    for mid, s in zip(sub.molecule_id, sub.smiles):
        if mid not in labels:                 # commit runs score a subset; the sample file lists every molecule
            continue
        truth = score_key(labels[mid])
        rr[mid] = 0.0
        for i, g in enumerate(s.split(";"), 1):
            if truth is not None and score_key(g) == truth:
                rr[mid] = 1.0 / i
                break
    return rr


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("stage_dir")
    ap.add_argument("--notebook", required=True)
    ap.add_argument("--fuse-dir", required=True)
    ap.add_argument("--pop-dir", help="folder with pool_lsid.npy and pool_lpmid.npy (POP_MU > 0)")
    ap.add_argument("--grid", action="append", default=[], metavar="KEY=v1,v2,...")
    ap.add_argument("--boot", type=int, default=2000)
    ap.add_argument("--head", help="head_scores.json from the DreaMS head; grid HEAD_LAM=... re-ranks same-formula groups")
    args = ap.parse_args()

    st = lambda n: json.load(open(os.path.join(args.stage_dir, f"{n}.json"))) \
        if os.path.exists(os.path.join(args.stage_dir, f"{n}.json")) else {}
    base, pc, eng, ice, gl, frag, meta, ours = (st(n) for n in ("base", "pc", "eng", "ice", "gl", "frag", "meta", "ours"))
    labels = pd.read_csv(os.path.join(args.stage_dir, "val_labels.csv")).set_index("molecule_id").smiles.to_dict()
    fc = load_fusion_core(args.notebook)
    sys.path.insert(0, args.fuse_dir)
    ice_fuse, gl_fuse = load_module(os.path.join(args.fuse_dir, "fuse.py"), "fuse"), \
        load_module(os.path.join(args.fuse_dir, "gl_fuse.py"), "gl_fuse")
    pool_pop = None
    if args.pop_dir:
        pool_pop = np.load(os.path.join(args.pop_dir, "pool_lsid.npy")).astype(np.float64) + \
            np.load(os.path.join(args.pop_dir, "pool_lpmid.npy")).astype(np.float64)
    score_key, formula = make_keys()
    head = json.load(open(args.head)) if args.head else {}

    grid = []
    for g in args.grid:
        k, vs = g.split("=", 1)
        grid.append([(k, ast.literal_eval(v)) for v in vs.split(",")])
    combos = [dict(c) for c in itertools.product(*grid)] or [{}]
    lib = {m: v[2] for m, v in base.items()}
    bins = {"lib<0.5": lambda x: x < 0.5, "0.5-0.9": lambda x: 0.5 <= x < 0.9, "lib>=0.9": lambda x: x >= 0.9}
    print(f"molecules {len(labels)} | our lists {len(ours)} | grid {len(combos)}")
    results = []
    for over in combos:
        over = dict(over); head_lam = float(over.pop("HEAD_LAM", 0.0)); head_max_lib = float(over.pop("HEAD_MAX_LIB", 9.0))
        cfg = dict(meta["cfg"]); cfg.update(over)
        sub, stats = fc.build_submission(base, pc, eng, ice, gl, cfg, meta["mol_order"], meta["sample_ids"],
                                         ice_fuse=ice_fuse, gl_fuse=gl_fuse, score_key=score_key, formula=formula,
                                         log=lambda *a: None, pool_pop=pool_pop, frag_scores=frag, ours=ours)
        if head_lam > 0:                      # same-formula re-rank of the final list, as ICEBERG/GLACIER do
            fused = []
            for mid, s in zip(sub.molecule_id, sub.smiles):
                smis = [x for x in s.split(";") if x]
                h = head.get(str(mid), {})
                if len(smis) > 1 and h and lib.get(mid, 0.0) < head_max_lib:   # library hits can be left alone
                    o = gl_fuse.rerank_multi(smis, [score_key(x) or x for x in smis], [float(len(smis) - i) for i in range(len(smis))],
                                             [formula(x) for x in smis], [h], [head_lam], top_n=len(smis))
                    smis = [smis[i] for i in o]
                fused.append(";".join(smis))
            sub = sub.assign(smiles=fused)
        over["HEAD_LAM"] = head_lam; over["HEAD_MAX_LIB"] = head_max_lib
        rr = reciprocal_ranks(sub, labels, score_key)
        mids = sorted(rr)
        row = {"cfg": over, "mrr": np.mean([rr[m] for m in mids]), "top1": np.mean([rr[m] == 1 for m in mids]),
               "hit25": np.mean([rr[m] > 0 for m in mids])}
        for b, f in bins.items():
            sel = [rr[m] for m in mids if f(lib.get(m, 0.0))]
            row[b] = (np.mean(sel) if sel else float("nan"), len(sel))
        row["ours_changed_top1"] = stats.get("ours", {}).get("top1_changed", 0)
        results.append((row, np.array([rr[m] for m in mids])))
    ref = results[0][1]
    rng = np.random.default_rng(0)
    idx = rng.integers(0, len(ref), size=(args.boot, len(ref)))
    for row, v in results:
        d = (v - ref)[idx].mean(1)
        lo, hi = np.percentile(d, [2.5, 97.5])
        print(f"{json.dumps(row['cfg']):40s} MRR {row['mrr']:.4f} top1 {row['top1']:.3f} hit25 {row['hit25']:.3f} | "
              + " ".join(f"{b} {row[b][0]:.3f}(n={row[b][1]})" for b in bins)
              + f" | top1 changed by ours {row['ours_changed_top1']} | vs ref {v.mean() - ref.mean():+.4f} [{lo:+.4f}, {hi:+.4f}]")


if __name__ == "__main__":
    main()
