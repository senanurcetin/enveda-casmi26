"""D3 — spectrum → fingerprint model and candidate scoring.

A small set transformer reads ≤ 64 peaks (no positional encoding: peaks are a set) plus one
global token (precursor mass, adduct, collision energy, polarity) and predicts the 2048-bit
Morgan fingerprint used by the pool, so a prediction can be scored directly against candidates.

m/z values need ~0.01 Da resolution over a 1,000 Da range, which a single linear projection cannot
represent; they are encoded with sinusoids at log-spaced wavelengths (as in DreaMS / MIST).
Each peak is described both by its m/z and by its neutral loss from the precursor.
"""
from dataclasses import asdict, dataclass

import numpy as np
import torch
import torch.nn as nn

N_ADDUCTS = 10


@dataclass
class FPModelConfig:
    n_bits: int = 2048
    d_model: int = 256
    n_layers: int = 4
    n_heads: int = 8
    d_ff: int = 1024
    dropout: float = 0.1
    n_freq: int = 32
    min_wavelength: float = 0.01
    max_wavelength: float = 1000.0


class FourierMass(nn.Module):
    def __init__(self, n_freq, lo, hi):
        super().__init__()
        wl = torch.logspace(np.log10(lo), np.log10(hi), n_freq)
        self.register_buffer("w", 2 * np.pi / wl)

    def forward(self, x):                      # x [..., ] -> [..., 2 * n_freq]
        a = x.unsqueeze(-1) * self.w
        return torch.cat([torch.sin(a), torch.cos(a)], dim=-1)


class FPModel(nn.Module):
    def __init__(self, cfg: FPModelConfig = FPModelConfig()):
        super().__init__()
        self.cfg = cfg
        f = 2 * cfg.n_freq
        self.enc = FourierMass(cfg.n_freq, cfg.min_wavelength, cfg.max_wavelength)
        self.peak_proj = nn.Sequential(nn.Linear(2 * f + 2, cfg.d_model), nn.GELU(), nn.Linear(cfg.d_model, cfg.d_model))
        self.adduct_emb = nn.Embedding(N_ADDUCTS, 32)
        self.glob_proj = nn.Sequential(nn.Linear(f + 32 + 2, cfg.d_model), nn.GELU(), nn.Linear(cfg.d_model, cfg.d_model))
        layer = nn.TransformerEncoderLayer(cfg.d_model, cfg.n_heads, cfg.d_ff, cfg.dropout,
                                           activation="gelu", batch_first=True, norm_first=True)
        self.encoder = nn.TransformerEncoder(layer, cfg.n_layers, enable_nested_tensor=False)
        self.head = nn.Sequential(nn.LayerNorm(cfg.d_model), nn.Linear(cfg.d_model, cfg.d_model),
                                  nn.GELU(), nn.Linear(cfg.d_model, cfg.n_bits))

    def forward(self, mz, it, npk, prec, adduct, ce, positive):
        """mz, it [B, P]; npk [B]; prec, ce, positive [B]; adduct [B] (long). Returns logits [B, n_bits]."""
        B, P = mz.shape
        pad = torch.arange(P, device=mz.device)[None, :] >= npk[:, None]
        loss = (prec[:, None] - mz).clamp(min=0)
        peak = torch.cat([self.enc(mz), self.enc(loss), it.sqrt().unsqueeze(-1), torch.log1p(it * 100).unsqueeze(-1)], -1)
        tok = self.peak_proj(peak)
        ce_n = torch.nan_to_num(ce, nan=30.0) / 100.0
        glob = torch.cat([self.enc(prec), self.adduct_emb(adduct), ce_n[:, None], positive[:, None]], -1)
        x = torch.cat([self.glob_proj(glob)[:, None, :], tok], dim=1)
        mask = torch.cat([torch.zeros(B, 1, dtype=torch.bool, device=mz.device), pad], dim=1)
        h = self.encoder(x, src_key_padding_mask=mask)
        return self.head(h[:, 0])


def save(model, path, extra=None):
    torch.save({"cfg": asdict(model.cfg), "state": model.state_dict(), **(extra or {})}, path)


def load(path, device="cpu"):
    ck = torch.load(path, map_location=device, weights_only=False)
    m = FPModel(FPModelConfig(**ck["cfg"])).to(device)
    m.load_state_dict(ck["state"])
    return m.eval()


def candidate_scores(logp_on, logp_off, cand_bits):
    """Bernoulli log-likelihood of each candidate fingerprint under the predicted bit probabilities.

    logp_on / logp_off: [n_bits] molecule-level log p(bit=1) / log p(bit=0)
    cand_bits: [n_cand, n_bits] {0,1}
    Returns (mean log-likelihood per bit, cosine between probabilities and candidate bits).
    """
    ll = (cand_bits * logp_on + (1 - cand_bits) * logp_off).mean(axis=1)
    p = np.exp(logp_on)
    cos = (cand_bits @ p) / (np.linalg.norm(p) * np.sqrt(cand_bits.sum(axis=1)) + 1e-9)
    return ll, cos
