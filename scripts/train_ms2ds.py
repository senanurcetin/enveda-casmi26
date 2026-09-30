"""Learned spectral similarity (MS2DeepScore-style): embed a spectrum so that the cosine between two
embeddings predicts the Tanimoto similarity of their structures.

Same encoder as D3, 128-d output. Batches pair each anchor spectrum with either another spectrum of
the same structure or a spectrum of one of its nearest structures, and the loss is a Tanimoto-bin
balanced MSE over all pairs in the batch (most random pairs are dissimilar). Trained per skeleton
fold like D3, so a query is only ever embedded by a model that never saw its structure.

  python scripts/train_ms2ds.py --fold 0      # -> models/ms2ds_fold0.pt (trained on fold 0)
"""
import argparse
import importlib.util
import sys
import time
from pathlib import Path

import numpy as np
import torch
import torch.nn.functional as F

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from casmi.config import ARTIFACT_DIR, SEED  # noqa: E402
from casmi.channels.fpmodel import FPModel, FPModelConfig, save  # noqa: E402

spec = importlib.util.spec_from_file_location("tfp", ROOT / "scripts" / "train_fpmodel.py")
tfp = importlib.util.module_from_spec(spec)
spec.loader.exec_module(tfp)
EMB = 128


def tanimoto(a, b):
    inter = a @ b.T
    return inter / (a.sum(1, keepdim=True) + b.sum(1)[None, :] - inter + 1e-9)


def neighbours(fp_bits, k=8, device="cuda", chunk=4096):
    """Top-k most similar other structures by Tanimoto (GPU, chunked)."""
    X = torch.as_tensor(fp_bits, dtype=torch.float16, device=device)
    out = np.zeros((len(X), k), np.int64)
    for i in range(0, len(X), chunk):
        T = tanimoto(X[i:i + chunk].float(), X.float())
        T[torch.arange(T.shape[0]), torch.arange(i, i + T.shape[0])] = -1
        out[i:i + chunk] = T.topk(k, dim=1).indices.cpu().numpy()
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fold", type=int, required=True)
    ap.add_argument("--steps", type=int, default=6000)
    ap.add_argument("--bs", type=int, default=256)
    ap.add_argument("--lr", type=float, default=5e-4)
    args = ap.parse_args()
    torch.manual_seed(SEED)
    rng = np.random.default_rng(SEED)
    device = "cuda"
    t0 = time.time()
    data = tfp.Data(ARTIFACT_DIR / "fpdata")
    side = data.side
    rows = np.flatnonzero(side.fold.to_numpy() == args.fold)
    structs = side.struct_idx.to_numpy()
    su = np.unique(structs[rows])
    bits = np.unpackbits(data.fp_packed[su], axis=1)
    nb = su[neighbours(bits, device=device)]                        # neighbour struct ids per struct
    nb_of = dict(zip(su, nb))
    order = np.argsort(structs[rows], kind="stable")
    r_sorted, s_sorted = rows[order], structs[rows][order]
    starts = np.searchsorted(s_sorted, su)
    ends = np.searchsorted(s_sorted, su, side="right")
    span = dict(zip(su, zip(starts, ends)))
    pick = lambda s: r_sorted[rng.integers(*span[s])]
    print(f"fold {args.fold}: {len(rows):,} spectra, {len(su):,} structures ({time.time() - t0:.0f}s)", flush=True)

    model = FPModel(FPModelConfig(n_bits=EMB)).to(device)
    opt = torch.optim.AdamW(model.parameters(), lr=args.lr, weight_decay=0.01)
    sched = torch.optim.lr_scheduler.OneCycleLR(opt, args.lr, total_steps=args.steps, pct_start=0.05)
    scaler = torch.amp.GradScaler("cuda")
    edges = torch.linspace(0, 1.0001, 11, device=device)
    run = 0.0
    for step in range(args.steps):
        a = rng.choice(rows, args.bs)
        sa = structs[a]
        same = rng.random(args.bs) < 0.3
        partner_struct = np.where(same, sa, np.array([nb_of[s][rng.integers(8)] for s in sa]))
        b = np.array([pick(s) for s in partner_struct])
        idx = np.r_[a, b]
        x, y, _ = data.batch(idx, device)            # batch() sorts idx; y follows the sorted order
        with torch.autocast("cuda", dtype=torch.float16):
            e = F.normalize(model(**x).float(), dim=1)
        target = tanimoto(y, y)
        pred = e @ e.T
        iu = torch.triu_indices(len(idx), len(idx), 1, device=device)
        t, p = target[iu[0], iu[1]], pred[iu[0], iu[1]]
        bin_ = (torch.bucketize(t.contiguous(), edges) - 1).clamp(0, 10)
        cnt = torch.bincount(bin_, minlength=11).float()
        w = 1.0 / cnt[bin_].clamp(min=1)
        loss = (w * (p - t) ** 2).sum() / w.sum()
        opt.zero_grad(set_to_none=True)
        scaler.scale(loss).backward()
        scaler.unscale_(opt)
        torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
        scaler.step(opt)
        scaler.update()
        sched.step()
        run += loss.item()
        if (step + 1) % 500 == 0:
            print(f"  step {step + 1}/{args.steps} balanced MSE {run / 500:.4f} ({time.time() - t0:.0f}s)", flush=True)
            run = 0.0
    out = ARTIFACT_DIR / "models" / f"ms2ds_fold{args.fold}.pt"
    save(model, out, {"emb": EMB, "train_fold": args.fold})
    print(f"saved {out.name} ({time.time() - t0:.0f}s)")


if __name__ == "__main__":
    main()
