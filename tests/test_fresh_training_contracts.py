"""Small real neural/ranker fits; no archived weights or competition records."""
import numpy as np
import pandas as pd
import torch

from casmi.channels.fpmodel import FPModel, FPModelConfig
from casmi.rank.ranker import Ranker, RankerConfig


def test_padding_is_not_evidence_and_fresh_network_has_finite_gradients():
    torch.manual_seed(17)
    model = FPModel(FPModelConfig(n_bits=16, d_model=32, n_layers=1, n_heads=4, d_ff=64))
    inputs = dict(mz=torch.tensor([[40., 60., 0., 0.]]), it=torch.tensor([[1., .5, 0., 0.]]),
                  npk=torch.tensor([2]), prec=torch.tensor([101.]), adduct=torch.tensor([0]),
                  ce=torch.tensor([20.]), positive=torch.tensor([1.]))
    model.eval()
    with torch.no_grad():
        expected = model(**inputs)
        changed = {**inputs, "mz": inputs["mz"].clone(), "it": inputs["it"].clone()}
        changed["mz"][0, 2:] = 999.
        changed["it"][0, 2:] = 99.
        torch.testing.assert_close(model(**changed), expected, atol=1e-6, rtol=1e-6)
    model.train()
    before = model.head[-1].weight.detach().clone()
    optimizer = torch.optim.AdamW(model.parameters(), lr=.001)
    loss = torch.nn.functional.binary_cross_entropy_with_logits(model(**inputs), torch.zeros(1, 16))
    loss.backward()
    assert all(torch.isfinite(parameter.grad).all() for parameter in model.parameters() if parameter.grad is not None)
    optimizer.step()
    assert not torch.equal(before, model.head[-1].weight)


def test_fresh_lambdarank_learns_grouped_candidates_and_records_dev_stopping():
    def table(prefix, count):
        return pd.DataFrame([{"qid": f"{prefix}-{q}", "cls": "C2", "label": int(rank == 0),
                              "d1_max": 1.0 if rank == 0 else .1, "abs_ppm": float(rank)}
                             for q in range(count) for rank in range(4)])
    config = RankerConfig(seeds=(7,), num_boost_round=15, early_stopping=3,
                          params={**RankerConfig().params, "min_data_in_leaf": 1,
                                  "num_threads": 1, "deterministic": True, "force_col_wise": True})
    model = Ranker(config).fit(table("train", 12), table("dev", 4), ["d1_max", "abs_ppm"])
    scores = model.predict(table("unseen", 1))
    assert np.isfinite(scores).all() and scores[0] > scores[1:].max()
    assert 1 <= model.best_iters[0] <= 15
