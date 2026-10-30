"""Independent assignment oracle and optional-baseline leakage contracts."""
import itertools
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

from casmi.spectra.optimal_cosine import optimal_cosine
from casmi.spectra.similarity import cosine_similarity
from casmi.channels.library import LibraryConfig, LibraryIndex


def test_global_assignment_recovers_peak_lost_by_proximity_greedy():
    ma, mb = np.array([100., 100.03]), np.array([100.01, 99.985])
    intensity = np.ones(2)
    assert optimal_cosine(ma, intensity, mb, intensity) == pytest.approx(1.)
    assert cosine_similarity(np.sort(ma), intensity, np.sort(mb), intensity) < 1.


def test_matches_exhaustive_oracle_and_is_symmetric_with_shift():
    rng = np.random.default_rng(42)
    for n, m in [(2, 3), (3, 3), (4, 2)]:
        for _ in range(12):
            ma, mb = rng.uniform(10, 11, n), rng.uniform(10, 11, m)
            a, b = rng.uniform(.1, 1, n), rng.uniform(.1, 1, m)
            delta = mb[None, :] - ma[:, None]
            w = np.outer(a, b) * ((abs(delta) <= .15) | (abs(delta - .3) <= .15))
            # Enumerate all injections in the smaller side, independently of SciPy.
            if n > m: w = w.T
            best = max(sum(w[i, j] for i, j in enumerate(p))
                       for p in itertools.permutations(range(w.shape[1]), w.shape[0]))
            expected = best / (np.linalg.norm(a) * np.linalg.norm(b))
            assert optimal_cosine(ma, a, mb, b, .15, .3) == pytest.approx(expected)
            assert optimal_cosine(mb, b, ma, a, .15, -.3) == pytest.approx(expected)


def test_zero_empty_scaled_and_unmatched_denominators():
    assert optimal_cosine([], [], [10], [1]) == 0
    assert optimal_cosine([10], [0], [10], [1]) == 0
    assert optimal_cosine([10, 20], [1, 1], [10], [1]) == pytest.approx(1 / np.sqrt(2))
    assert optimal_cosine([10, 20], [1e300, 2e300], [20, 10], [2., 1.]) == pytest.approx(1.)
    # The same pair eligible in direct and shifted graphs must never count twice.
    assert optimal_cosine([10], [1], [10], [1], shift=.01) == 1.


@pytest.mark.parametrize('kwargs', [{'tol_da': -1}, {'shift': np.nan}, {'intensity_power': 0}, {'max_peaks': True}])
def test_invalid_controls_fail(kwargs):
    with pytest.raises(ValueError): optimal_cosine([10], [1], [10], [1], **kwargs)


def test_malformed_peaks_and_resource_limit_fail():
    for mz, intensity in [([10], [1, 2]), ([np.nan], [1]), ([10], [-1]), ([[10]], [[1]])]:
        with pytest.raises(ValueError): optimal_cosine(mz, intensity, [10], [1])
    with pytest.raises(ValueError, match='budget'):
        optimal_cosine(np.arange(257), np.ones(257), [10], [1])


def test_library_baseline_keeps_exclusions_and_mode_boundary():
    meta = pd.DataFrame({'ppm': [0., 0.], 'n_peaks': [3, 3], 'is_timstof': [True, True],
        'adduct': ['[M+H]+', '[M-H]-'], 'ionization_mode': ['positive', 'negative'], 'precursor_mz': [300., 300.]})
    mz, intensity = np.array([50., 100., 150.]), np.ones(3)
    store = SimpleNamespace(meta=meta, peaks=lambda row: (mz, intensity))
    index = LibraryIndex(store, np.array(['key', 'key']), LibraryConfig(similarity='optimal_cosine'))
    query = [{'mode': 'positive', 'adduct': '[M+H]+', 'mz_c': mz, 'it_c': intensity}]
    assert index.score(query, 'key', [])[0] == pytest.approx(1.)
    assert index.score(query, 'key', [0]) == (0., 0., 0)
    assert LibraryConfig().similarity == 'entropy'
    with pytest.raises(ValueError): LibraryConfig(similarity='unknown')
