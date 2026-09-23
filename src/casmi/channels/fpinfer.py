"""D3 inference shared by validation (scripts/add_d3_features.py) and the submission."""
import numpy as np
import torch

from casmi.chem.masses import TEST_ADDUCTS
from casmi.spectra.clean import CleanConfig, clean_peaks_fast

MAX_PEAKS = 64
FP_CLEAN = CleanConfig(min_rel_intensity=0.005, top_k_per_window=8, window_da=50, max_peaks=MAX_PEAKS)
ADDUCT_ID = {a: i for i, a in enumerate(TEST_ADDUCTS)}


@torch.no_grad()
def molecule_logp(models, spectra, device):
    """Mean predicted bit probability over a molecule's spectra and over models.
    Returns (log p(bit=1), log p(bit=0)) as NumPy arrays of length n_bits."""
    n = len(spectra)
    mz, it, npk = np.zeros((n, MAX_PEAKS), np.float32), np.zeros((n, MAX_PEAKS), np.float32), np.zeros(n, np.int64)
    for j, s in enumerate(spectra):
        cm, ci = clean_peaks_fast(s["mz"], s["it"], s["precursor_mz"], FP_CLEAN)
        mz[j, :cm.size], it[j, :cm.size], npk[j] = cm, ci, cm.size
    t = lambda a, dt=torch.float32: torch.as_tensor(np.asarray(a), dtype=dt, device=device)
    x = dict(mz=t(mz), it=t(it), npk=t(npk, torch.long),
             prec=t([s["precursor_mz"] for s in spectra]),
             adduct=t([ADDUCT_ID.get(s["adduct"], 0) for s in spectra], torch.long),
             ce=t([s["ce_ev"] for s in spectra]), positive=t([float(s["mode"] == "positive") for s in spectra]))
    p = torch.stack([torch.sigmoid(m(**x)) for m in models]).mean(dim=(0, 1)).clamp(1e-4, 1 - 1e-4)
    return torch.log(p).cpu().numpy(), torch.log1p(-p).cpu().numpy()
