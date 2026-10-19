import numpy as np
import pytest

from casmi.spectra.clean import CleanConfig, clean_peaks
from casmi.spectra.similarity import cosine_similarity, entropy_similarity

MZ = np.array([100.05, 150.10, 200.20, 250.30], np.float32)
IT = np.array([0.2, 1.0, 0.5, 0.05], np.float32)


def test_identical_spectra_score_one():
    assert entropy_similarity(MZ, IT, MZ, IT) == pytest.approx(1.0, abs=1e-9)
    assert cosine_similarity(MZ, IT, MZ, IT) == pytest.approx(1.0, abs=1e-9)


@pytest.mark.parametrize("similarity", [entropy_similarity, cosine_similarity])
def test_blank_zero_intensity_and_invalid_peaks_never_return_nan(similarity):
    empty = np.array([], dtype=np.float32)
    zero = np.zeros(MZ.size, dtype=np.float32)
    assert similarity(empty, empty, MZ, IT) == 0.0
    assert similarity(MZ, zero, MZ, IT) == 0.0
    assert similarity(MZ, IT, MZ, zero) == 0.0
    for invalid_mz, invalid_it in [
        (np.array([np.nan, 150., 200., 250.], dtype=np.float32), IT),
        (np.array([-3., 150., 200., 250.], dtype=np.float32), IT),
        (MZ, np.array([np.nan, 1., 0.5, 0.1], dtype=np.float32)),
        (MZ, np.array([-1., 1., 0.5, 0.1], dtype=np.float32)),
        (MZ, np.array([np.inf, 1., 0.5, 0.1], dtype=np.float32)),
    ]:
        assert similarity(invalid_mz, invalid_it, MZ, IT) == 0.0
        assert similarity(MZ, IT, invalid_mz, invalid_it) == 0.0


@pytest.mark.parametrize("similarity", [entropy_similarity, cosine_similarity])
def test_valid_zero_intensity_peak_contributes_no_nan(similarity):
    mz = np.array([100., 110.], dtype=np.float32)
    intensity = np.array([1., 0.], dtype=np.float32)
    result = similarity(mz, intensity, mz, intensity)
    assert np.isfinite(result)
    assert result == pytest.approx(1.0, abs=1e-9)


@pytest.mark.parametrize("similarity", [entropy_similarity, cosine_similarity])
def test_mismatched_mass_intensity_lengths_fail_closed(similarity):
    with pytest.raises(ValueError, match="aligned 1D"):
        similarity(MZ, IT[:-1], MZ, IT)


def test_disjoint_spectra_score_zero():
    other = MZ + 7.0
    assert entropy_similarity(MZ, IT, other, IT) == 0.0


def test_similarity_is_symmetric_and_bounded():
    b_mz = np.array([100.05, 150.11, 300.0], np.float32)
    b_it = np.array([1.0, 0.3, 0.4], np.float32)
    s1 = entropy_similarity(MZ, IT, b_mz, b_it)
    s2 = entropy_similarity(b_mz, b_it, MZ, IT)
    assert s1 == pytest.approx(s2, abs=1e-12) and 0.0 < s1 < 1.0


def test_mass_shift_recovers_analog():
    # an analog carrying +14.016 (CH2) on the fragments that hold the substituent
    analog = MZ.copy()
    analog[2:] += 14.01565
    plain = entropy_similarity(MZ, IT, analog, IT)
    shifted = entropy_similarity(MZ, IT, analog, IT, 0.02, 14.01565)
    assert shifted > plain and shifted == pytest.approx(1.0, abs=1e-9)


def test_cleaning_drops_above_precursor_and_noise():
    mz = np.array([50.0, 50.002, 120.0, 180.0, 260.0])
    it = np.array([0.0005, 0.0005, 1.0, 0.3, 0.9])
    cmz, cit = clean_peaks(mz, it, precursor_mz=200.0, cfg=CleanConfig(min_rel_intensity=0.01))
    assert list(np.round(cmz, 3)) == [120.0, 180.0]      # 260 > precursor + 2, 50.x merged and below floor
    assert cit.max() == pytest.approx(1.0)


def test_windowed_top_k():
    mz = np.array([10.0, 11.0, 12.0, 60.0, 61.0])
    it = np.array([0.9, 0.5, 1.0, 0.2, 0.1])
    cmz, _ = clean_peaks(mz, it, None, CleanConfig(min_rel_intensity=0, top_k_per_window=2, window_da=50, merge_tol_da=0))
    assert list(cmz) == [10.0, 12.0, 60.0, 61.0]


@pytest.mark.parametrize("cfg", [
    CleanConfig(),
    CleanConfig(min_rel_intensity=0.01, top_k_per_window=6, window_da=50, max_peaks=64),
    CleanConfig(merge_tol_da=0.0, max_peaks=20),
])
def test_fast_cleaner_matches_reference(cfg):
    from casmi.spectra.clean import clean_peaks_fast
    rng = np.random.default_rng(0)
    for _ in range(300):
        n = int(rng.integers(1, 400))
        mz = np.sort(rng.uniform(40, 600, n)) + rng.normal(0, 0.002, n)
        it = rng.lognormal(0, 1.5, n)
        prec = float(rng.uniform(300, 620))
        a_mz, a_it = clean_peaks(mz, it, prec, cfg)
        b_mz, b_it = clean_peaks_fast(mz, it, prec, cfg)
        np.testing.assert_allclose(a_mz, b_mz, rtol=1e-6)
        np.testing.assert_allclose(a_it, b_it, rtol=1e-5)
