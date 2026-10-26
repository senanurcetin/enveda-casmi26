import numpy as np

from scripts.run_advanced_study import SearchMemo, ablation_columns


class Index:
    def __init__(self):
        self.calls = 0
    def search(self, spectra, mass, excluded_rows):
        self.calls += 1
        return {"candidate": (float(self.calls), 100.0)}


def test_search_memo_cannot_cross_evidence_or_knockout_boundaries():
    index = Index()
    memo = SearchMemo(index)
    spectra = [{"adduct": "[M+H]+", "mode": "positive", "precursor_mz": 101.0,
                "mz_c": np.array([40., 60.]), "it_c": np.array([1., .5])}]
    first = memo.search(spectra, 100., np.array([1, 2]))
    assert memo.search(spectra, 100., np.array([2, 1])) == first
    first.clear()
    assert memo.search(spectra, 100., np.array([1, 2]))
    memo.search(spectra, 100., np.array([1]))
    spectra[0]["it_c"][1] = .7
    memo.search(spectra, 100., np.array([1]))
    assert index.calls == 3
    memo.reset()
    memo.search(spectra, 100., np.array([1]))
    assert index.calls == 4


def test_ablations_remove_channel_context_and_default_source_membership():
    columns = ["d1_max", "d1_max_rank", "an_p4", "an_p4_gap", "d3_ll", "d3_ll_z",
               "d4_int", "gen", "gen_depth", "np_score", "abs_ppm", "prior_log_spectra"]
    assert "prior_log_spectra" not in ablation_columns(columns, "full")
    for name, prefix in (("no_d1", "d1_"), ("no_d2", "an_"), ("no_d3", "d3_"),
                         ("no_d4", "d4_"), ("no_gen_features", "gen")):
        assert not any(c.startswith(prefix) for c in ablation_columns(columns, name))
    assert ablation_columns(columns, "membership_only") == ["prior_log_spectra"]
