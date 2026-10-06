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
    _csv = _mod.run_ours(os.path.join(COMP, 'test.parquet'), os.path.join(WORK, 'ours_submission.csv'), log=print)
    if _csv is not None and os.path.exists(_csv):
        _o = pd.read_csv(_csv)
        OURS = {str(m): [x for x in str(s).split(';') if x and x != 'CCO'] for m, s in zip(_o.molecule_id, _o.smiles)}
    log(f'OURS: {len(OURS)} molecules, {sum(len(v) for v in OURS.values())} candidates, {time.time()-_t0:.0f}s')
    assert IS_RERUN or OURS, 'OURS channel produced nothing on the commit run: fix before submitting'
'''



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
