"""Write a Kaggle script kernel that trains the DreaMS isomer-ranking head (scripts/dreams_head.py) on a T4.

usage: python kaggle/build_head_kernel.py OUT_DIR --docker-from BASE_KERNEL_METADATA.json [--epochs 3] [--id senanuretin/NAME]

The kernel installs the competition's pinned RDKit wheel, runs `prep` and `train`, and leaves head_ep*.pt plus the
log (hold-out MRR, blank-spectrum control, chance level) in /kaggle/working/head_art. Internet stays off.
The docker image is copied from the base notebook's metadata so the Python version matches the cp312 RDKit wheel.
"""
import argparse
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
HEAD_SRC = (HERE.parent / "scripts" / "dreams_head.py").read_text()

KERNEL = '''import glob, os, subprocess, sys, time
T0 = time.time()
log = lambda m: print(f"[{{time.time() - T0:7.0f}}s] {{m}}", flush=True)
whl = sorted(w for w in glob.glob("/kaggle/input/**/rdkit-2026.3.3-*.whl", recursive=True)
             if f"cp{{sys.version_info.major}}{{sys.version_info.minor}}" in w)
log(f"python {{sys.version.split()[0]}} | rdkit wheels {{whl}}")
subprocess.run([sys.executable, "-m", "pip", "install", "-q", "--no-index", "--no-deps", whl[0]], check=True)
DREAMS = os.path.dirname(sorted(glob.glob("/kaggle/input/**/Enveda180_DreaMS_index.parquet", recursive=True))[0])
ART = "/kaggle/working/head_art"
open("/kaggle/working/dreams_head.py", "w").write({src!r})
run = lambda *a: subprocess.run([sys.executable, "/kaggle/working/dreams_head.py", *a], check=True)
log("prep"); run("prep", "--dreams", DREAMS, "--out", ART)
log("train"); run("train", "--dreams", DREAMS, "--art", ART, "--epochs", "{epochs}", "--fp-in-ram", "--per-mol", "{per_mol}",
                  "--decoys", "{decoys}", "--batch", "{batch}", "--lr", "{lr}")
for f in ("mol_fp.npy", "isomers.json", "spectra.parquet", "mols.parquet"):      # keep the output small
    os.remove(os.path.join(ART, f))
log("done")
'''


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("out_dir")
    ap.add_argument("--docker-from", required=True)
    ap.add_argument("--id", default="senanuretin/casmi26-dreams-head")
    ap.add_argument("--epochs", type=int, default=3)
    ap.add_argument("--per-mol", type=int, default=2)
    ap.add_argument("--decoys", type=int, default=31)
    ap.add_argument("--batch", type=int, default=256)
    ap.add_argument("--lr", type=float, default=3e-4)
    a = ap.parse_args()
    out = Path(a.out_dir); out.mkdir(parents=True, exist_ok=True)
    (out / "train_head.py").write_text(KERNEL.format(src=HEAD_SRC, epochs=a.epochs, per_mol=a.per_mol,
                                                    decoys=a.decoys, batch=a.batch, lr=a.lr))
    base = json.load(open(a.docker_from))
    meta = dict(id=a.id, title=a.id.split("/")[1], code_file="train_head.py", language="python", kernel_type="script",
                is_private=True, enable_gpu=True, enable_tpu=False, enable_internet=False,
                dataset_sources=["hengck23/hengck23-dreams-enveda-casmi26", "metric/rdkit-2026-3-3-wheel"],
                competition_sources=[], kernel_sources=[], model_sources=[],
                docker_image=base.get("docker_image"), machine_shape=base.get("machine_shape", "NvidiaTeslaT4"))
    json.dump(meta, open(out / "kernel-metadata.json", "w"), indent=1)
    print("wrote", out)


if __name__ == "__main__":
    main()
