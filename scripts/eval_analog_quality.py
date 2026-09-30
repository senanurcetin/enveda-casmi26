"""Does the learned similarity find better analogs than entropy similarity?

For V-in and V-real queries in the Class-2 setting (every library spectrum of the truth removed),
both methods retrieve analog structures within ±200 Da; the score is how close the best of the top-k
analogs is to the truth (Tanimoto). Better analogs are what D2 propagation and the Class-3 hedge
live on. Each query is embedded by the fold model that never saw its structure.
"""
import importlib.util
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import torch
import torch.nn.functional as F

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from casmi.config import ARTIFACT_DIR  # noqa: E402
from casmi.channels.analog import AnalogIndex  # noqa: E402
from casmi.channels.fpmodel import load  # noqa: E402
from casmi.channels.library import LibraryIndex  # noqa: E402
from casmi.chem.canon import skeleton_key  # noqa: E402
from casmi.engine import EngineConfig, molecule_mass  # noqa: E402
from casmi.pool.pool import CandidatePool  # noqa: E402
from casmi.spectra.clean import clean_peaks_fast  # noqa: E402
from casmi.spectra.store import SpectrumStore  # noqa: E402
from casmi.validate.knockout import build_queries  # noqa: E402


def _mod(name, file):
    s = importlib.util.spec_from_file_location(name, ROOT / "scripts" / file)
    m = importlib.util.module_from_spec(s)
    s.loader.exec_module(m)
    return m


brd, tfp = _mod("brd", "build_rank_data.py"), _mod("tfp", "train_fpmodel.py")


@torch.no_grad()
def embed_all(model, data, device, bs=4096):
    out = np.zeros((len(data.side), 128), np.float16)
    for i in range(0, len(data.side), bs):
        idx = np.arange(i, min(i + bs, len(data.side)))
        x, _, _ = data.batch(idx, device)
        with torch.autocast("cuda", dtype=torch.float16):
            out[idx] = F.normalize(model(**x).float(), dim=1).cpu().numpy()
    return out


def main():
    t0 = time.time()
    device = "cuda"
    cfg = EngineConfig()
    store = SpectrumStore(ARTIFACT_DIR / "store")
    key_of_row = store.meta.train_ik14.to_numpy(dtype=object)
    pool = CandidatePool(pd.read_parquet(ARTIFACT_DIR / "pool_pcl.parquet"))
    fps = np.load(ARTIFACT_DIR / "pool_fp2048_pcl.npy")
    library = LibraryIndex(store, key_of_row, cfg.library)
    analog = AnalogIndex(store, key_of_row, cfg.analog, cache=ARTIFACT_DIR / "analog_index.npz")
    data = tfp.Data(ARTIFACT_DIR / "fpdata")
    side = data.side
    skey = pd.read_parquet(ARTIFACT_DIR / "fpdata" / "structs.parquet").pool_key.to_numpy()[side.struct_idx.to_numpy()]
    fp_of_row = {r: i for i, r in enumerate(side.store_row.to_numpy())}
    lib_mass = side.precursor_mz.to_numpy().astype(np.float64)
    vin = build_queries(store, key_of_row)
    vin_skel = {skeleton_key(pool.smiles_of(q.pool_key)) for q in vin}
    vreal, stratum_of = brd.vreal_queries(store, key_of_row, vin_skel)
    items = [(q, -1, "vin") for q in vin] + [(q, f, stratum_of[q.pool_key]) for q, f in vreal]
    bits = lambda keys: np.unpackbits(fps[[pool.row_of[k] for k in keys]], axis=1).astype(np.float32)

    def tani(truth, keys):
        if not keys:
            return np.zeros(0)
        a, b = bits([truth]), bits(keys)
        inter = b @ a[0]
        return inter / (a.sum() + b.sum(1) - inter + 1e-9)

    rows = []
    for model_fold in (0, 1):
        model = load(ARTIFACT_DIR / "models" / f"ms2ds_fold{model_fold}.pt", device)
        E = torch.as_tensor(embed_all(model, data, device), device=device)
        print(f"embedded {len(E):,} spectra with ms2ds_fold{model_fold} ({time.time() - t0:.0f}s)", flush=True)
        for q, fold, st in items:
            use = 0 if fold == -1 else 1 - fold
            if use != model_fold or q.pool_key not in pool.row_of:
                continue
            qi = [fp_of_row[r] for r in q.spectrum_rows if r in fp_of_row]
            if not qi:
                continue
            mass = molecule_mass(q.spectra)
            ok = (np.abs(lib_mass - (q.spectra[0]["precursor_mz"])) <= 200) & (skey != q.pool_key)
            sims = (E[qi].float() @ E.float().T).max(0).values.cpu().numpy()
            sims[~ok] = -1
            top = np.argsort(-sims)[:400]
            ms_keys = list(dict.fromkeys(k for k in skey[top] if k in pool.row_of))[:10]
            # entropy analogs, same exclusion (every library spectrum of the truth)
            spectra = [dict(s) for s in q.spectra]
            for s in spectra:
                s["mz_c"], s["it_c"] = clean_peaks_fast(s["mz"], s["it"], s["precursor_mz"], library.cfg.clean)
            excl = np.asarray(library.rows_of.get(q.pool_key, []), np.int64)
            ent = analog.search(spectra, mass, excl)
            en_keys = [k for k, _ in sorted(ent.items(), key=lambda kv: -kv[1][0]) if k != q.pool_key and k in pool.row_of][:10]
            for name, keys in (("entropy", en_keys), ("learned", ms_keys)):
                t = tani(q.pool_key, keys)
                rows.append({"stratum": st, "method": name, **{f"best@{k}": float(t[:k].max()) if len(t) else 0.0 for k in (1, 5, 10)}})
        del E
        torch.cuda.empty_cache()
    r = pd.DataFrame(rows)
    print(r.groupby(["stratum", "method"]).mean().round(3).to_string())
    r2 = r.assign(close=r["best@5"] >= 0.7).groupby(["stratum", "method"]).close.mean().round(3)
    print("\nshare of queries with an analog of Tanimoto >= 0.7 in the top 5:\n" + r2.to_string())
    print(f"({time.time() - t0:.0f}s)")


if __name__ == "__main__":
    main()
