"""Patch a copy of the public v4n fusion notebook so our own pipeline becomes a third list in its final fusion.

Usage: python build_fusion_notebook.py SRC.ipynb OUT_DIR [--beta 0.4] [--max-lib 0.7] [--id senanuretin/NAME]

Only our additions live in this repository; SRC.ipynb (Apache-2.0, by huseyinemreaksoy, built on seyitkaangunes'
and others' public notebooks) is read at build time and is not committed.

What changes in the copy
  1. fusion_core.py (embedded): `ours` list per molecule, weighted reciprocal-rank fused into the final list
     of molecules whose library match is weak (lib_max < OURS_MAX_LIB). Library hits are never touched.
  2. CFG: OURS_BETA / OURS_MAX_LIB / OURS_KRR / OURS_DEPTH (OURS_BETA = 0 switches the channel off).
  3. A new cell before the fusion cell runs kaggle/ours_channel.py (fail-open) and loads its lists.
  4. kernel-metadata: our four datasets are attached.
"""
import argparse
import ast
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
OURS_SRC = (HERE.parent / "kaggle" / "ours_channel.py").read_text()

FUSE_FN = '''

def fuse_ours(vs, ours, score_key, beta, krr, depth, n=25):
    """Fused list vs (host ranking) with ours: sum 1/(krr + r) over vs plus beta/(krr + r) over our list, metric-key dedup."""
    sc, smi_of = {}, {}
    for r, s in enumerate(vs, 1):
        k = score_key(s) or s
        sc[k] = sc.get(k, 0.0) + 1.0 / (krr + r); smi_of.setdefault(k, s)
    for r, s in enumerate(ours[:depth], 1):
        k = score_key(s) or s
        sc[k] = sc.get(k, 0.0) + beta / (krr + r); smi_of.setdefault(k, s)
    order = sorted(sc, key=lambda k: -sc[k])[:n]
    return [smi_of[k] for k in order]
'''

OURS_CELL = '''# OURS: our retrieve-and-rank pipeline (casmi26 repo, kaggle/ours_channel.py) as a third list in the final fusion.
# Runs after ICEBERG / GLACIER have released the GPU. Fail-open: any error -> OURS = {} -> host ranking unchanged.
OURS = {}
if CFG['OURS_BETA'] > 0:
    _t0 = time.time()
    import importlib.util as _ilu
    _spec = _ilu.spec_from_file_location('ours_channel', os.path.join(OURS_DIR, 'ours_channel.py'))
    _mod = _ilu.module_from_spec(_spec); _spec.loader.exec_module(_mod)
    _csv = _mod.run_ours(os.path.join(COMP, 'test.parquet'), os.path.join(WORK, 'ours_submission.csv'), log=print,
                         data_dir=COMP if CFG['VALIDATION'] else None,   # validation: library from the purged train
                         fp_keep=('fp_fold0.pt',) if CFG['VALIDATION'] and CFG['VAL_SET'] == 'test_like' else None)
    if _csv is not None and os.path.exists(_csv):
        _o = pd.read_csv(_csv)
        OURS = {str(m): [x for x in str(s).split(';') if x and x != 'CCO'] for m, s in zip(_o.molecule_id, _o.smiles)}
    log(f'OURS: {len(OURS)} molecules, {sum(len(v) for v in OURS.values())} candidates, {time.time()-_t0:.0f}s')
    assert IS_RERUN or OURS, 'OURS channel produced nothing on the commit run: fix before submitting'
    dump('ours', OURS)                        # stage cache: lets fusion weights be replayed offline
'''


TEST_LIKE_FN = '''
def test_like_keys(n_np, n_pc, seed, max_mass=1170.0, min_mass=150.0):
    """Test-like hold-out (added by casmi26 build_fusion_notebook.py --val-set test_like): structures that NEITHER
    pipeline trained on, with answers reachable the way hidden-test answers are.
      * in the fold-safe FPNet's 'hold' fold (host FP_BANK 'fold0' never saw them),
      * skeleton fold 1 of our 2-fold split (md5(skeleton) % 2), so our fp_fold0 model is out-of-fold for them,
      * in the host pool through an external source (src bit 2), not only because they are training structures,
      * one eligible spectrum at least (test adduct, base peak >= 1000), mass in the test window.
    Two strata: 'np' = in a natural-product catalogue (pool_np > 0), 'pc' = external but not an NP catalogue.
    With the v2 pool the external source is COCONUT, so 'pc' is empty: a PubChem-only answer cannot be simulated
    fairly here (every training structure stays in the host pool through src bit 1). Returns {inchikey14: stratum}."""
    import numpy as np, hashlib as _h
    from rdkit import Chem
    def _skel(smi):
        m = Chem.MolFromSmiles(smi) if isinstance(smi, str) else None
        if m is None:
            return None
        rw = Chem.RWMol(Chem.RemoveHs(m, sanitize=False))
        for a in rw.GetAtoms():
            a.SetFormalCharge(0); a.SetIsAromatic(False); a.SetNoImplicit(True); a.SetNumExplicitHs(0)
            a.SetChiralTag(Chem.ChiralType.CHI_UNSPECIFIED)
        for b in rw.GetBonds():
            b.SetBondType(Chem.BondType.SINGLE); b.SetIsAromatic(False); b.SetStereo(Chem.BondStereo.STEREONONE)
        return Chem.MolToSmiles(rw, canonical=True)
    sp = pd.read_parquet(find('casmi26-fold-safe-fpnet/**/split.parquet'))
    pm = pd.read_parquet(os.path.join(POOL_DIR, 'pool_meta.parquet'), columns=['key', 'pid', 'src'])
    npf = np.load(find('casmi26-pubchem-popularity-prior/**/pool_np.npy'))
    ext = pm[(pm.src.values & 2) > 0]
    nat = set(ext.key.values[npf[ext.pid.values] > 0]); ext = set(ext.key.values)
    st = pd.read_parquet(os.path.join(POOL_DIR, 'train_structs.parquet'), columns=['inchikey14', 'key', 'smiles', 'mass'])
    st = st[st.inchikey14.isin(set(sp.inchikey14[sp.fold == 'hold'])) & st.key.isin(ext) & st.mass.between(min_mass, max_mass)]
    import pyarrow.parquet as _pq
    m = _pq.read_table(os.path.join(COMP, 'train.parquet'), columns=['inchikey14', 'adduct', 'base_peak_intensity']).to_pandas()
    m = m[m.adduct.isin(TEST_ADDUCTS) & (m.base_peak_intensity.isna() | (m.base_peak_intensity >= 1000))]
    st = st[st.inchikey14.isin(set(m.inchikey14))]
    rng = np.random.RandomState(seed)
    st = st.iloc[rng.permutation(len(st))]
    out, seen, cnt, want = {}, set(), {'np': 0, 'pc': 0}, {'np': n_np, 'pc': n_pc}
    for ik, key, smi in zip(st.inchikey14, st.key, st.smiles):
        s_ = 'np' if key in nat else 'pc'
        if cnt[s_] >= want[s_] or key in seen:
            continue
        sk = _skel(smi)
        if sk is None or int(_h.md5(sk.encode()).hexdigest(), 16) % 2 != 1:
            continue
        out[ik] = s_; seen.add(key); cnt[s_] += 1
        if cnt == want:
            break
    log(f'test_like: {cnt} structures from {len(st)} eligible')
    return out

'''


def patch_validation_cell(src):
    a = "def build_validation("
    assert a in src
    src = src.replace(a, TEST_LIKE_FN.lstrip("\n") + "\n" + a, 1)
    b = "    if CFG['VAL_SET'] == 'fold0_np':"
    assert b in src
    src = src.replace(b, (
        "    if CFG['VAL_SET'] == 'test_like':\n"
        "        _tl = test_like_keys(CFG['VAL_N_NP'], CFG['VAL_N'] - CFG['VAL_N_NP'], CFG['VAL_SEED'])\n"
        "        _nm, _ns = build_validation(COMP, os.path.join(WORK, 'val_comp'), None, len(_tl), CFG['VAL_SEED'],\n"
        "                                    CFG['VAL_MAX_SPEC'], 'all', only=set(_tl))\n"
        "        _lab = pd.read_csv(os.path.join(STAGE, 'val_labels.csv'))\n"
        "        _lab['stratum'] = _lab.inchikey14.map(_tl)\n"
        "        _lab.to_csv(os.path.join(STAGE, 'val_labels.csv'), index=False)\n"
        "    el" + b.lstrip()), 1)
    return src


def patch_fusion_core(src):
    a = "                PC_FWD_T1=None, PC_FWD_N=25)"
    assert a in src
    src = src.replace(a, "                PC_FWD_T1=None, PC_FWD_N=25,\n"
                         "                OURS_BETA=0.0, OURS_MAX_LIB=0.7, OURS_KRR=3.0, OURS_DEPTH=25)")
    b = "def fuse_rrf("
    assert b in src
    src = src.replace(b, FUSE_FN.lstrip("\n") + "\n\n" + b, 1)
    c = "                     frag_scores=None):"
    assert c in src
    src = src.replace(c, "                     frag_scores=None, ours=None):")
    d = "        vs = fill(mid, vs[:25])\n        out.append((mid, ';'.join(vs[:25]) if vs else c['FALLBACK']))"
    assert d in src
    new = ("        o = (ours or {}).get(str(mid))\n"
           "        if o and c['OURS_BETA'] > 0 and vs and lib_of.get(mid, 0.0) < c['OURS_MAX_LIB']:\n"
           "            f = fuse_ours(vs[:25], o, score_key, c['OURS_BETA'], c['OURS_KRR'], c['OURS_DEPTH'])\n"
           "            stats.setdefault('ours', {'molecules': 0, 'top1_changed': 0})\n"
           "            stats['ours']['molecules'] += 1; stats['ours']['top1_changed'] += int(f[:1] != vs[:1])\n"
           "            vs = f\n") + d
    src = src.replace(d, new)
    return src


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("src"); ap.add_argument("out_dir")
    ap.add_argument("--beta", type=float, default=0.4)
    ap.add_argument("--max-lib", type=float, default=0.7)
    ap.add_argument("--id", default="senanuretin/casmi26-fusion-ours")
    ap.add_argument("--meta", help="kernel-metadata.json of the source notebook")
    ap.add_argument("--validation", action="store_true", help="commit-run validation (class-2-like), fold-safe FPNet")
    ap.add_argument("--val-set", choices=["np_examples", "test_like"], default="np_examples",
                    help="np_examples: famous NPs (all in COCONUT, optimistic for our channel); test_like: structures "
                         "neither pipeline trained on, NP-catalogue and PubChem-type strata")
    ap.add_argument("--val-n", type=int, default=250); ap.add_argument("--val-n-np", type=int, default=250)
    ap.add_argument("--add", action="append", default=[], help="extra dataset source (e.g. a re-uploaded copy of one of ours); repeatable")
    ap.add_argument("--drop", action="append", default=[], help="dataset source to remove (validation-only inputs); repeatable")
    args = ap.parse_args()

    nb = json.load(open(args.src))
    cells = nb["cells"]
    src2 = "".join(cells[2]["source"])
    tree = ast.parse(src2)
    node = next(n for n in tree.body if isinstance(n, ast.Assign) and getattr(n.targets[0], "id", "") == "EMBED")
    embed = ast.literal_eval(node.value)               # data only: nothing from the source notebook is executed
    embed["fusion_core.py"] = patch_fusion_core(embed["fusion_core.py"])
    header = src2.split("EMBED = ")[0]
    cells[2]["source"] = (header + "EMBED = " + repr(embed) + "\n").splitlines(keepends=True)

    cfg = "".join(cells[4]["source"])
    cfg += (f"\nCFG.update({{'OURS_BETA': {args.beta}, 'OURS_MAX_LIB': {args.max_lib}, 'OURS_KRR': 3.0, 'OURS_DEPTH': 25}})"
            "   # our own pipeline as a third list; OURS_BETA 0 = off\n")
    if args.validation and args.val_set == "np_examples":
        cfg += ("CFG.update({'VALIDATION': True, 'VAL_SET': 'np_examples', 'VAL_PURGE': 'all', 'FP_BANK': 'fold0'})"
                "   # validation run: never submit this version\n")
    elif args.validation:
        cfg += (f"CFG.update({{'VALIDATION': True, 'VAL_SET': 'test_like', 'VAL_PURGE': 'all', 'FP_BANK': 'fold0', "
                f"'VAL_N': {args.val_n}, 'VAL_N_NP': {args.val_n_np}, 'VAL_MAX_SPEC': 6}})"
                "   # validation run: never submit this version\n")
        vi = next(i for i, c in enumerate(cells) if "def build_validation(" in "".join(c["source"]))
        cells[vi]["source"] = patch_validation_cell("".join(cells[vi]["source"])).splitlines(keepends=True)
    cells[4]["source"] = cfg.splitlines(keepends=True)

    idx = next(i for i, c in enumerate(cells) if c["cell_type"] == "code" and "fusion_core.build_submission(" in "".join(c["source"]))
    fcell = "".join(cells[idx]["source"])
    assert "frag_scores=FRAG)" in fcell
    cells[idx]["source"] = fcell.replace("frag_scores=FRAG)", "frag_scores=FRAG, ours=OURS)").splitlines(keepends=True)

    ours_cell = (f"OURS_DIR = '/kaggle/working/ours_code'\n"
                 "os.makedirs(OURS_DIR, exist_ok=True)\n"
                 f"open(os.path.join(OURS_DIR, 'ours_channel.py'), 'w').write({OURS_SRC!r})\n" + OURS_CELL)
    cells.insert(idx, {"cell_type": "code", "id": "ours-channel", "metadata": {}, "execution_count": None, "outputs": [],
                       "source": ours_cell.splitlines(keepends=True)})   # nbformat 4.5: every cell needs an id

    out = Path(args.out_dir); out.mkdir(parents=True, exist_ok=True)
    json.dump(nb, open(out / "notebook.ipynb", "w"))
    meta = json.load(open(args.meta)) if args.meta else {}
    meta.pop("id_no", None)
    meta.update(id=args.id, title=args.id.split("/")[1], code_file="notebook.ipynb", is_private=True,
                enable_gpu=True, enable_internet=False)
    ds = list(meta.get("dataset_sources", []))
    for d in ["senanuretin/casmi-code", "senanuretin/casmi-models", "senanuretin/casmi-coconut", "senanuretin/casmi-pubchemlite"]:
        if d not in ds:
            ds.append(d)
    ds += [d for d in args.add if d not in ds]
    meta["dataset_sources"] = [d for d in ds if d not in args.drop]
    json.dump(meta, open(out / "kernel-metadata.json", "w"), indent=1)
    print("wrote", out, "| fusion cell at index", idx, "| datasets", len(meta["dataset_sources"]))


if __name__ == "__main__":
    sys.exit(main())
