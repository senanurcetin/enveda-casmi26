"""Train the D3 spectrum → fingerprint model on one fold (or all folds for the final model).

  python scripts/train_fpmodel.py --train-folds 0      # predicts fold 1 + V-in out-of-fold
  python scripts/train_fpmodel.py --train-folds 1
  python scripts/train_fpmodel.py --train-folds 0,1,-1 # final model (includes V-in)

Validation is on structures the model never saw (the other fold), reported as BCE and as the
retrieval quantity that matters: how often the true structure's fingerprint out-scores a decoy of
similar mass under the predicted probabilities.
"""
import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import torch
import torch.nn.functional as F

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from casmi.config import ARTIFACT_DIR, SEED  # noqa: E402
from casmi.channels.fpmodel import FPModel, FPModelConfig, save  # noqa: E402


class Data:
    def __init__(self, root):
        self.mz = np.load(root / "mz.npy", mmap_mode="r")
        self.it = np.load(root / "it.npy", mmap_mode="r")
        self.npk = np.load(root / "npk.npy")
        self.side = pd.read_parquet(root / "side.parquet")
        self.fp_packed = np.load(root / "struct_fp.npy")                     # [n_struct, 256] packed bits

    def batch(self, idx, device):
        idx = np.sort(idx)                                                     # memmap locality
        s = self.side.iloc[idx]
        t = lambda a, dt=torch.float32: torch.as_tensor(np.ascontiguousarray(a), dtype=dt, device=device)
        x = dict(mz=t(self.mz[idx]), it=t(self.it[idx].astype(np.float32)), npk=t(self.npk[idx], torch.long),
                 prec=t(s.precursor_mz.to_numpy()), adduct=t(s.adduct_id.to_numpy(), torch.long),
                 ce=t(s.ce_ev.to_numpy()), positive=t(s.positive.to_numpy()))
        y = t(np.unpackbits(self.fp_packed[s.struct_idx.to_numpy()], axis=1))
        return x, y, s


@torch.no_grad()
def evaluate(model, data, idx, device, bs=2048, n_decoys=20, seed=0):
    model.eval()
    rng = np.random.default_rng(seed)
    losses, wins, n = [], 0, 0
    structs = data.side.struct_idx.to_numpy()
    for i in range(0, len(idx), bs):
        x, y, s = data.batch(idx[i:i + bs], device)
        with torch.autocast("cuda", dtype=torch.float16):
            logits = model(**x)
        logits = logits.float()
        losses.append(F.binary_cross_entropy_with_logits(logits, y).item())
        # decoys: fingerprints of random other structures in the same evaluation set
        lp_on, lp_off = F.logsigmoid(logits), F.logsigmoid(-logits)
        true_ll = (y * lp_on + (1 - y) * lp_off).mean(1)
        pick = rng.choice(structs[idx], size=(len(s), n_decoys))
        dec = torch.as_tensor(np.unpackbits(data.fp_packed[pick], axis=-1), dtype=torch.float32, device=device)
        dec_ll = (dec * lp_on[:, None] + (1 - dec) * lp_off[:, None]).mean(2)
        wins += (true_ll[:, None] > dec_ll).float().mean(1).sum().item()
        n += len(s)
    model.train()
    return float(np.mean(losses)), wins / n


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--train-folds", default="0")
    ap.add_argument("--epochs", type=int, default=8)
    ap.add_argument("--bs", type=int, default=768)
    ap.add_argument("--lr", type=float, default=3e-4)
    ap.add_argument("--timstof-weight", type=float, default=2.0, help="sampling weight of instrument-matched spectra")
    args = ap.parse_args()

    torch.manual_seed(SEED)
    device = "cuda"
    data = Data(ARTIFACT_DIR / "fpdata")
    folds = [int(f) for f in args.train_folds.split(",")]
    fold = data.side.fold.to_numpy()
    tr = np.flatnonzero(np.isin(fold, folds))
    held = [f for f in (0, 1) if f not in folds]
    rng = np.random.default_rng(SEED)
    va = rng.choice(np.flatnonzero(np.isin(fold, held)), 20000, replace=False) if held else \
        rng.choice(tr, 20000, replace=False)
    w = np.where(data.side.is_timstof.to_numpy()[tr] == 1, args.timstof_weight, 1.0)
    w = w / w.sum()

    model = FPModel(FPModelConfig()).to(device)
    n_params = sum(p.numel() for p in model.parameters())
    steps_per_epoch = len(tr) // args.bs
    opt = torch.optim.AdamW(model.parameters(), lr=args.lr, weight_decay=0.01)
    sched = torch.optim.lr_scheduler.OneCycleLR(opt, args.lr, total_steps=args.epochs * steps_per_epoch, pct_start=0.05)
    scaler = torch.amp.GradScaler("cuda")
    print(f"train folds {folds}: {len(tr):,} spectra | validate on {'fold ' + str(held) if held else 'train sample'} "
          f"| {n_params / 1e6:.1f}M params | {steps_per_epoch} steps/epoch", flush=True)

    tag = "fp_fold" + "_".join(str(f) for f in folds)
    out = ARTIFACT_DIR / "models"
    out.mkdir(exist_ok=True)
    history, t0 = [], time.time()
    for ep in range(args.epochs):
        order = rng.choice(tr, size=steps_per_epoch * args.bs, replace=True, p=w)
        run = 0.0
        for st in range(steps_per_epoch):
            x, y, _ = data.batch(order[st * args.bs:(st + 1) * args.bs], device)
            with torch.autocast("cuda", dtype=torch.float16):
                loss = F.binary_cross_entropy_with_logits(model(**x).float(), y)
            opt.zero_grad(set_to_none=True)
            scaler.scale(loss).backward()
            scaler.unscale_(opt)
            torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
            scaler.step(opt)
            scaler.update()
            sched.step()
            run += loss.item()
            if (st + 1) % 500 == 0:
                print(f"  ep {ep + 1} step {st + 1}/{steps_per_epoch} loss {run / 500:.4f} "
                      f"({(st + 1) * args.bs / (time.time() - t0 - sum(h['sec'] for h in history)):.0f} spectra/s)", flush=True)
                run = 0.0
        v_loss, v_win = evaluate(model, data, va, device)
        history.append({"epoch": ep + 1, "val_bce": v_loss, "val_win_vs_decoy": v_win,
                        "sec": time.time() - t0 - sum(h["sec"] for h in history)})
        print(f"epoch {ep + 1}: val BCE {v_loss:.4f} | truth beats random decoy {v_win:.3f} | "
              f"{history[-1]['sec'] / 60:.1f} min", flush=True)
        save(model, out / f"{tag}.pt", {"history": history, "train_folds": folds})
    (out / f"{tag}.json").write_text(json.dumps({"args": vars(args), "history": history}, indent=2))
    print(f"saved models/{tag}.pt ({(time.time() - t0) / 60:.0f} min)")


if __name__ == "__main__":
    main()
