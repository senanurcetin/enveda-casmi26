"""Write a Kaggle script kernel that scores a validation run's candidate lists with the trained DreaMS head.

usage: python kaggle/build_score_kernel.py OUT_DIR --docker-from BASE_KERNEL_METADATA.json
           [--head-kernel senanuretin/casmi26-dreams-head] [--val-kernel senanuretin/casmi26-val-ours] [--weights head_ep2.pt]

Inputs are the outputs of two earlier kernels (kernel_sources): the head's weights and the validation run's
stage_cache + val_comp/test.parquet. The kernel
  1. checks that our DreaMS preprocessing reproduces the published Enveda-180 embeddings (cosine on 64 spectra),
  2. embeds the held-out test spectra on the GPU,
  3. writes head_scores.json ({molecule_id: {smiles: score}}, the ICE/GL dump format) for offline replay.
"""
import argparse
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
HEAD_SRC = (HERE.parent / "scripts" / "dreams_head.py").read_text()

KERNEL = '''import glob, json, os, subprocess, sys, time
import numpy as np, pandas as pd, pyarrow.parquet as pq
T0 = time.time()
log = lambda m: print(f"[{{time.time() - T0:7.0f}}s] {{m}}", flush=True)
whl = sorted(w for w in glob.glob("/kaggle/input/**/rdkit-2026.3.3-*.whl", recursive=True)
             if f"cp{{sys.version_info.major}}{{sys.version_info.minor}}" in w)
subprocess.run([sys.executable, "-m", "pip", "install", "-q", "--no-index", "--no-deps", whl[0]], check=True)
open("/kaggle/working/dreams_head.py", "w").write({src!r})
sys.path.insert(0, "/kaggle/working")
import dreams_head as H
first = lambda pat: sorted(glob.glob(pat, recursive=True), key=len)[0]
DREAMS = os.path.dirname(first("/kaggle/input/**/Enveda180_DreaMS_index.parquet"))
WEIGHTS = first("/kaggle/input/**/head_art/{weights}")
STAGE = os.path.dirname(first("/kaggle/input/**/stage_cache/val_labels.csv"))
TEST = first("/kaggle/input/**/val_comp/test.parquet")
TRAIN = first("/kaggle/input/**/enveda-CASMI26-molecule-id-mass-spectra/train.parquet")
log(f"dreams {{DREAMS}} | weights {{WEIGHTS}} | stage {{STAGE}} | test {{TEST}}")

# 1. parity of our preprocessing with the published embeddings
mod, model = H.load_dreams(DREAMS, device="cuda")
idx = pd.read_parquet(os.path.join(DREAMS, "Enveda180_DreaMS_index.parquet"), columns=["train_idx"])
pick = np.sort(np.random.default_rng(0).choice(len(idx), 64, replace=False))
want = set(idx.train_idx.values[pick].tolist())
pf = pq.ParquetFile(TRAIN); rows, start = {{}}, 0
for rg in range(pf.num_row_groups):
    n = pf.metadata.row_group(rg).num_rows
    hit = [r for r in want if start <= r < start + n]
    if hit:
        t = pf.read_row_group(rg, columns=["ms2_mzs", "ms2_normalized_intensities", "precursor_mz"]).to_pandas()
        for r in hit:
            rows[r] = (np.asarray(t.ms2_mzs.iat[r - start]), np.asarray(t.ms2_normalized_intensities.iat[r - start]),
                       float(t.precursor_mz.iat[r - start]))
    start += n
mine = H.embed_spectra(mod, model, [rows[r] for r in idx.train_idx.values[pick]], device="cuda")
ref = np.load(os.path.join(DREAMS, "Enveda180_DreaMS_embeddings.npy"), mmap_mode="r")[pick]
cos = (mine * ref).sum(1) / np.linalg.norm(mine, axis=1) / np.linalg.norm(ref, axis=1)
log(f"PARITY cosine to the published embeddings: min {{cos.min():.5f}} median {{np.median(cos):.5f}}")

# 2. embed the held-out spectra, 3. score every candidate list of the stage dump
subprocess.run([sys.executable, "/kaggle/working/dreams_head.py", "embed", "--dreams", DREAMS, "--parquet", TEST,
                "--out", "/kaggle/working/val_emb.npy", "--device", "cuda"], check=True)
subprocess.run([sys.executable, "/kaggle/working/dreams_head.py", "score", "--weights", WEIGHTS,
                "--emb", "/kaggle/working/val_emb.npy", "--stage", STAGE, "--out", "/kaggle/working/head_scores.json"], check=True)
log("done")
'''


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("out_dir")
    ap.add_argument("--docker-from", required=True)
    ap.add_argument("--id", default="senanuretin/casmi26-score-head")
    ap.add_argument("--head-kernel", default="senanuretin/casmi26-dreams-head")
    ap.add_argument("--val-kernel", default="senanuretin/casmi26-val-ours")
    ap.add_argument("--weights", default="head_ep2.pt")
    a = ap.parse_args()
    out = Path(a.out_dir); out.mkdir(parents=True, exist_ok=True)
    (out / "score_head.py").write_text(KERNEL.format(src=HEAD_SRC, weights=a.weights))
    base = json.load(open(a.docker_from))
    meta = dict(id=a.id, title=a.id.split("/")[1], code_file="score_head.py", language="python", kernel_type="script",
                is_private=True, enable_gpu=True, enable_tpu=False, enable_internet=False,
                dataset_sources=["hengck23/hengck23-dreams-enveda-casmi26", "metric/rdkit-2026-3-3-wheel"],
                competition_sources=["enveda-CASMI26-molecule-id-mass-spectra"],
                kernel_sources=[a.head_kernel, a.val_kernel], model_sources=[],
                docker_image=base.get("docker_image"), machine_shape=base.get("machine_shape", "NvidiaTeslaT4"))
    json.dump(meta, open(out / "kernel-metadata.json", "w"), indent=1)
    print("wrote", out)


if __name__ == "__main__":
    main()
