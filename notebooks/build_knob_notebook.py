"""Copy of the public v4n fusion notebook with only CFG overrides changed (no code, no extra datasets).

Usage: python build_knob_notebook.py SRC.ipynb OUT_DIR --meta KERNEL-METADATA.json --id senanuretin/NAME --set GL_MH_ONLY=True [--set POP_MU=0.3 ...]

The base notebook is the one reproducible, platform-stable submission we have (public LB 0.420); every
variant that attached our own datasets failed on Kaggle's rerun ('A system error'). This builder keeps the
15 public datasets and changes one knob per submission, so the public LB can judge each one.
Source notebook: huseyinemreaksoy/casmi26-v4n-fusion-pubchem-on-public-0-421 (Apache-2.0); not committed here.
"""
import argparse
import ast
import json
from pathlib import Path


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("src"); ap.add_argument("out_dir")
    ap.add_argument("--meta", required=True)
    ap.add_argument("--id", required=True)
    ap.add_argument("--set", action="append", default=[], metavar="KEY=VALUE", help="CFG override, VALUE is a Python literal")
    args = ap.parse_args()

    over = {}
    for kv in args.set:
        k, v = kv.split("=", 1)
        over[k] = ast.literal_eval(v)
    nb = json.load(open(args.src))
    cfg = "".join(nb["cells"][4]["source"])
    assert "CFG = dict(" in cfg
    cfg += f"\nCFG.update({over!r})   # knob variant: {', '.join(args.set)}\n"
    nb["cells"][4]["source"] = cfg.splitlines(keepends=True)

    out = Path(args.out_dir); out.mkdir(parents=True, exist_ok=True)
    json.dump(nb, open(out / "notebook.ipynb", "w"))
    meta = json.load(open(args.meta)); meta.pop("id_no", None)
    meta.update(id=args.id, title=args.id.split("/")[1], code_file="notebook.ipynb", is_private=True,
                enable_gpu=True, enable_internet=False)
    json.dump(meta, open(out / "kernel-metadata.json", "w"), indent=1)
    print("wrote", out, "| overrides", over, "| datasets", len(meta["dataset_sources"]))


if __name__ == "__main__":
    main()
