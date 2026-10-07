"""Isomer-ranking head on frozen DreaMS spectrum embeddings (candidate-InfoNCE, same-formula decoys).

Why: the public pipeline loses most of its score ordering isomers (same formula) of Class-2 molecules, and its
fingerprint heads are trained with per-bit BCE, which is not a ranking objective. This head is trained directly to
put the true structure first among its isomers (MSAlign-style candidate-InfoNCE, arXiv 2605.19752).

Shortcut guard: every decoy is another Enveda-180 structure with the same molecular formula, i.e. drawn from the
same library as the truths, so "looks like the Enveda library" cannot separate truth from decoys. External
candidates (COCONUT / PubChem) are never used as decoys in training: with synthetic truths they would teach
"natural-product-like = wrong", the opposite of what the hidden test needs.

Data: hengck23/hengck23-dreams-enveda-casmi26 (DreaMS embeddings of the Enveda-180 spectra, row-aligned with
Enveda180_DreaMS_index.parquet). Only test adducts are used (no [2M+...] dimers).

  python scripts/dreams_head.py prep  --dreams DIR --out ART
  python scripts/dreams_head.py train --dreams DIR --art ART [--epochs 3]
"""
import argparse
import hashlib
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd

TEST_ADDUCTS = ['[M+H]+', '[M+NH4]+', '[M-H2O+H]+', '[M-2H2O+H]+', '[M+Na]+', '[M+K]+',
                '[M-H]-', '[M-H2O-H]-', '[M+CH2O2-H]-', '[M+Cl]-']
FP_BITS, FP_RADIUS = 2048, 2


def morgan_counts(smiles):
    """uint8 Morgan count fingerprints (radius 2, 2048 folded bits, counts clipped at 255); zeros if unparsable."""
    from rdkit import Chem, RDLogger
    from rdkit.Chem import rdFingerprintGenerator
    RDLogger.DisableLog("rdApp.*")
    gen = rdFingerprintGenerator.GetMorganGenerator(radius=FP_RADIUS, fpSize=FP_BITS)
    out = np.zeros((len(smiles), FP_BITS), np.uint8)
    for i, s in enumerate(smiles):
        m = Chem.MolFromSmiles(s)
        if m is not None:
            out[i] = np.minimum(gen.GetCountFingerprintAsNumPy(m), 255)
    return out


def holdout(key, frac=0.05):
    return int(hashlib.md5(key.encode()).hexdigest()[:8], 16) / 0xFFFFFFFF < frac


def prep(args):
    t0 = time.time()
    out = Path(args.out); out.mkdir(parents=True, exist_ok=True)
    idx = pd.read_parquet(Path(args.dreams) / "Enveda180_DreaMS_index.parquet",
                          columns=["inchikey14", "normalized_smiles", "molecular_formula", "adduct", "ionization_mode"])
    idx["row"] = np.arange(len(idx))                      # embedding row = index row
    spec = idx[idx.adduct.isin(TEST_ADDUCTS)].copy()
    mols = spec.drop_duplicates("inchikey14")[["inchikey14", "normalized_smiles", "molecular_formula"]].reset_index(drop=True)
    mols["mol_id"] = np.arange(len(mols))
    mols["holdout"] = mols.inchikey14.map(holdout)
    print(f"spectra with test adducts {len(spec):,} of {len(idx):,} | structures {len(mols):,} "
          f"(hold-out {mols.holdout.sum():,}) | {time.time() - t0:.0f}s", flush=True)
    fp = morgan_counts(mols.normalized_smiles.tolist())
    np.save(out / "mol_fp.npy", fp)
    print(f"fingerprints {fp.shape}, empty {int((fp.sum(1) == 0).sum())} | {time.time() - t0:.0f}s", flush=True)
    spec = spec.merge(mols[["inchikey14", "mol_id", "holdout"]], on="inchikey14")
    spec["neg"] = (spec.ionization_mode == "negative").astype(np.int8)
    spec[["row", "mol_id", "holdout", "neg"]].to_parquet(out / "spectra.parquet", index=False)
    mols.to_parquet(out / "mols.parquet", index=False)
    grp = mols.groupby("molecular_formula").mol_id.apply(list)
    iso = {int(m): [int(x) for x in g if x != m] for g in grp for m in g}
    n_iso = np.array([len(v) for v in iso.values()])
    json.dump({str(k): v for k, v in iso.items()}, open(out / "isomers.json", "w"))
    print(f"isomers per structure: median {np.median(n_iso):.0f}, share with none {np.mean(n_iso == 0):.3f} "
          f"| {time.time() - t0:.0f}s")


# ---------------------------------------------------------------- model / training (torch imported lazily)
def build_model(d_spec=1024, d_out=256, hidden=1024, dropout=0.1):
    import torch
    from torch import nn

    class Head(nn.Module):
        def __init__(self):
            super().__init__()
            self.spec = nn.Sequential(nn.Linear(d_spec + 1, hidden), nn.GELU(), nn.Dropout(dropout), nn.Linear(hidden, d_out))
            self.mol = nn.Sequential(nn.Linear(FP_BITS, hidden), nn.GELU(), nn.Dropout(dropout), nn.Linear(hidden, d_out))
            self.logit_scale = nn.Parameter(torch.tensor(2.66))          # 1 / 0.07

        def embed_spec(self, emb, neg):
            return nn.functional.normalize(self.spec(torch.cat([emb, neg[:, None]], 1)), dim=-1)

        def embed_mol(self, fp):
            return nn.functional.normalize(self.mol(torch.log1p(fp)), dim=-1)

        def forward(self, emb, neg, fp_cands):                           # fp_cands [B, C, BITS], truth at column 0
            s = self.embed_spec(emb, neg)
            m = self.embed_mol(fp_cands.flatten(0, 1)).view(fp_cands.shape[0], fp_cands.shape[1], -1)
            return torch.einsum("bd,bcd->bc", s, m) * self.logit_scale.clamp(max=4.6).exp()

    return Head()


def load_art(args):
    art = Path(args.art)
    spec = pd.read_parquet(art / "spectra.parquet")
    fp = np.load(art / "mol_fp.npy", mmap_mode="r")
    iso = {int(k): np.asarray(v, np.int64) for k, v in json.load(open(art / "isomers.json")).items()}
    emb = np.load(Path(args.dreams) / "Enveda180_DreaMS_embeddings.npy", mmap_mode="r")
    return spec, fp, iso, emb


def sample_decoys(rng, mol_id, iso, n_mols, k):
    pool = iso.get(int(mol_id), np.zeros(0, np.int64))
    if len(pool) == 0:                                   # no isomer in the library: random structures
        return rng.integers(0, n_mols, k)
    return pool[rng.integers(0, len(pool), k)] if len(pool) < k else rng.choice(pool, k, replace=False)


def train(args):
    import torch
    torch.set_num_threads(args.threads); torch.manual_seed(0)
    dev = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    rng = np.random.default_rng(0)
    spec, fp, iso, emb = load_art(args)
    tr = spec[spec.holdout == 0]
    fp_t = torch.from_numpy(np.asarray(fp, dtype=np.float32)) if args.fp_in_ram else None
    if fp_t is not None:
        fp_t = fp_t.to(dev)
    model = build_model().to(dev); opt = torch.optim.AdamW(model.parameters(), lr=args.lr, weight_decay=1e-4)
    print(f"device {dev} | training spectra {len(tr):,}", flush=True)
    n_mols, K, B = fp.shape[0], args.decoys, args.batch
    t0 = time.time()
    for ep in range(args.epochs):
        # at most --per-mol spectra of each structure per epoch, shuffled
        ep_rows = tr.groupby("mol_id", group_keys=False).apply(lambda g: g.sample(min(len(g), args.per_mol), random_state=ep))
        ep_rows = ep_rows.sample(frac=1.0, random_state=ep).reset_index(drop=True)
        model.train(); tot, n = 0.0, 0
        for st in range(0, len(ep_rows) - B + 1, B):
            b = ep_rows.iloc[st:st + B]
            cands = np.stack([np.concatenate([[m], sample_decoys(rng, m, iso, n_mols, K)]) for m in b.mol_id.values])
            x = torch.from_numpy(np.asarray(emb[np.sort(b.row.values)], dtype=np.float32))
            order = np.argsort(np.argsort(b.row.values)); x = x[torch.from_numpy(order)]       # back to batch order
            x = x.to(dev); neg = torch.from_numpy(b.neg.values.astype(np.float32)).to(dev)
            f = fp_t[torch.from_numpy(cands).to(dev)] if fp_t is not None else torch.from_numpy(np.asarray(fp[cands.ravel()], np.float32)).view(B, K + 1, -1).to(dev)
            logits = model(x, neg, f)
            loss = torch.nn.functional.cross_entropy(logits, torch.zeros(B, dtype=torch.long, device=dev))
            opt.zero_grad(); loss.backward(); opt.step()
            tot += float(loss) * B; n += B
            if (st // B) % 200 == 0:
                print(f"epoch {ep} step {st // B}/{len(ep_rows) // B} loss {tot / n:.4f} scale {model.logit_scale.exp().item():.1f} "
                      f"| {time.time() - t0:.0f}s", flush=True)
        torch.save({k: v.cpu() for k, v in model.state_dict().items()}, Path(args.art) / f"head_ep{ep}.pt")
        evaluate(args, model=model.cpu(), spec=spec, fp=fp, iso=iso, emb=emb, tag=f"epoch {ep}")
        model.to(dev)


def evaluate(args, model=None, spec=None, fp=None, iso=None, emb=None, tag="eval"):
    """Hold-out structures: rank the truth among ALL its same-formula Enveda-180 isomers (up to --max-cands).
    Also a blank-spectrum control: the same ranking with the spectrum embedding zeroed, which must fall to chance
    if the head reads the spectrum rather than the molecule."""
    import torch
    if model is None:
        spec, fp, iso, emb = load_art(args)
        model = build_model(); model.load_state_dict(torch.load(args.weights, map_location="cpu"))
    model.eval()
    ho = spec[spec.holdout == 1]
    ho = ho.groupby("mol_id", group_keys=False).apply(lambda g: g.head(2))         # <= 2 spectra per structure
    rr, rr_blank, rr_rand, ncand = [], [], [], []
    rng = np.random.default_rng(1)
    with torch.no_grad():
        for r in ho.itertuples():
            pool = iso.get(int(r.mol_id), np.zeros(0, np.int64))
            if len(pool) == 0:
                continue
            cands = np.concatenate([[r.mol_id], pool[:args.max_cands - 1]])
            f = torch.from_numpy(np.asarray(fp[cands], np.float32))[None]
            x = torch.from_numpy(np.asarray(emb[r.row], np.float32))[None]
            neg = torch.tensor([float(r.neg)])
            for xx, acc in ((x, rr), (torch.zeros_like(x), rr_blank)):
                s = model(xx, neg, f)[0].numpy()
                rank = 1 + int((s[1:] > s[0]).sum()) + 0.5 * int((s[1:] == s[0]).sum())
                acc.append(1.0 / rank)
            rr_rand.append(np.mean(1.0 / np.arange(1, len(cands) + 1)))
            ncand.append(len(cands))
    rr, rr_blank = np.array(rr), np.array(rr_blank)
    print(f"[{tag}] hold-out spectra {len(rr):,} | candidates median {int(np.median(ncand))} | "
          f"MRR {rr.mean():.4f} top1 {np.mean(rr == 1):.3f} | blank-spectrum MRR {rr_blank.mean():.4f} | "
          f"chance MRR {np.mean(rr_rand):.4f}", flush=True)
    return rr.mean()


# ---------------------------------------------------------------- inference: embed spectra, score candidate lists
def load_dreams(dreams_dir, device="cpu"):
    import importlib.util
    spec = importlib.util.spec_from_file_location("dreams_minimal_model", str(Path(dreams_dir) / "dreams_minimal_model.py"))
    mod = importlib.util.module_from_spec(spec); spec.loader.exec_module(mod)
    model = mod.load_model(str(Path(dreams_dir) / "embedding_model_state_dict.pt"), device=device).eval()
    return mod, model


def embed_spectra(mod, model, rows, batch=64, device="cpu"):
    """rows: iterable of (mzs, intensities, precursor_mz) -> [N, 1024] float32 DreaMS embeddings."""
    import torch
    rows = list(rows); out = []
    with torch.no_grad():
        for i in range(0, len(rows), batch):
            x = torch.from_numpy(mod.preprocess_batch(rows[i:i + batch])).to(device)
            out.append(model(x).float().cpu().numpy())
    return np.concatenate(out) if out else np.zeros((0, 1024), np.float32)


def embed(args):
    """Embed a test-format parquet (molecule_id, ms2_mzs, ms2_normalized_intensities, precursor_mz, ionization_mode)."""
    mod, model = load_dreams(args.dreams)
    df = pd.read_parquet(args.parquet)
    t0 = time.time()
    e = embed_spectra(mod, model, zip(df.ms2_mzs, df.ms2_normalized_intensities, df.precursor_mz))
    np.save(args.out, e)
    df[["molecule_id", "spectrum_id", "adduct", "ionization_mode"]].to_parquet(str(args.out).replace(".npy", "_meta.parquet"), index=False)
    print(f"embedded {len(df):,} spectra in {time.time() - t0:.0f}s -> {args.out}")


def score_lists(model, emb, neg, smiles):
    """Mean over a molecule's spectra of the head's cosine score for each candidate SMILES."""
    import torch
    fp = torch.from_numpy(morgan_counts(smiles).astype(np.float32))
    with torch.no_grad():
        s = model.embed_spec(torch.from_numpy(emb), torch.from_numpy(neg.astype(np.float32)))
        m = model.embed_mol(fp)
        return (s @ m.T).mean(0).numpy()


def score(args):
    """Head scores for every candidate in a stage dump: {mid: {smiles: score}} (the format of the ICE/GL dumps)."""
    import torch
    model = build_model(); model.load_state_dict(torch.load(args.weights, map_location="cpu")); model.eval()
    emb = np.load(args.emb); meta = pd.read_parquet(args.emb.replace(".npy", "_meta.parquet"))
    stage = Path(args.stage)
    lists = {}
    for name in ("base", "eng", "pc"):
        f = stage / f"{name}.json"
        if not f.exists():
            continue
        for m, v in json.load(open(f)).items():
            smis = v[0] if name == "base" else (v.get("smiles") if name == "eng" else v.get("pc")) or []
            lists.setdefault(str(m), []).extend(smis)
    ours = stage / "ours.json"
    if ours.exists():
        for m, v in json.load(open(ours)).items():
            lists.setdefault(str(m), []).extend(v)
    out, t0 = {}, time.time()
    for m, idx in meta.groupby("molecule_id").indices.items():
        c = list(dict.fromkeys(lists.get(str(m), [])))
        if not c:
            continue
        neg = (meta.ionization_mode.values[idx] == "negative")
        out[str(m)] = dict(zip(c, map(float, score_lists(model, emb[idx], neg, c))))
    json.dump(out, open(args.out, "w"))
    print(f"scored {sum(len(v) for v in out.values()):,} candidates of {len(out)} molecules in {time.time() - t0:.0f}s -> {args.out}")


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("prep"); p.add_argument("--dreams", required=True); p.add_argument("--out", required=True)
    for name in ("train", "eval"):
        p = sub.add_parser(name)
        p.add_argument("--dreams", required=True); p.add_argument("--art", required=True)
        p.add_argument("--epochs", type=int, default=3); p.add_argument("--batch", type=int, default=256)
        p.add_argument("--decoys", type=int, default=31); p.add_argument("--per-mol", type=int, default=2)
        p.add_argument("--lr", type=float, default=3e-4); p.add_argument("--threads", type=int, default=4)
        p.add_argument("--max-cands", type=int, default=200); p.add_argument("--weights")
        p.add_argument("--fp-in-ram", action="store_true", help="hold all fingerprints as float32 in RAM (~1.5 GB)")
    p = sub.add_parser("embed"); p.add_argument("--dreams", required=True); p.add_argument("--parquet", required=True)
    p.add_argument("--out", required=True)
    p = sub.add_parser("score"); p.add_argument("--weights", required=True); p.add_argument("--emb", required=True)
    p.add_argument("--stage", required=True); p.add_argument("--out", required=True)
    args = ap.parse_args()
    {"prep": prep, "train": train, "eval": evaluate, "embed": embed, "score": score}[args.cmd](args)


if __name__ == "__main__":
    main()
