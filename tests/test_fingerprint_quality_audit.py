"""Read-only historical fingerprint audit: aggregate diagnostics only."""
import numpy as np
import pytest

from casmi.pool.fingerprint_audit import fingerprint_quality_summary


def test_zero_vectors_are_counted_without_exposing_structure_rows():
    keys = np.array(["PRIVATE_TEST_MOLECULE_ALPHA", "PRIVATE_TEST_MOLECULE_BETA",
                     "PRIVATE_TEST_MOLECULE_GAMMA", "PRIVATE_TEST_MOLECULE_DELTA"])
    vectors = np.zeros((4, 256), dtype=np.uint8)
    vectors[0, 0] = 1
    vectors[2, 7] = 128
    report = fingerprint_quality_summary(keys, keys.copy(), vectors, batch_rows=2)
    assert report["status"] == "INVALID_ZERO_FINGERPRINTS"
    assert report["candidate_count"] == 4
    assert report["zero_fingerprint_rows"] == 2
    assert report["zero_fraction"] == 0.5
    assert all(key not in str(report) for key in keys)


def test_valid_and_empty_vector_collections_are_stable():
    keys = np.array(["A", "B"])
    fps = np.ones((2, 256), np.uint8)
    result = fingerprint_quality_summary(keys, keys, fps, batch_rows=1)
    assert result["status"] == "PASS_NONZERO_FINGERPRINTS"
    assert result["zero_fraction"] == 0.0
    with pytest.raises(ValueError, match="Empty candidate pools"):
        fingerprint_quality_summary(
            np.array([]), np.array([]), np.zeros((0, 256), np.uint8)
        )


@pytest.mark.parametrize("bad", [0, -1, 1.5, True, 1000001])
def test_invalid_chunk_requests_fail_closed(bad):
    with pytest.raises(ValueError, match="batch_rows"):
        fingerprint_quality_summary([], [], np.zeros((0, 256), np.uint8),
                                    batch_rows=bad)


def test_bad_shapes_key_order_and_types_cannot_claim_validity():
    keys = np.array(["A", "B"])
    fps = np.ones((2, 256), np.uint8)
    with pytest.raises(ValueError, match="key order"):
        fingerprint_quality_summary(keys, keys[::-1], fps)
    with pytest.raises(ValueError, match="row counts"):
        fingerprint_quality_summary(keys[:1], keys, fps)
    with pytest.raises(ValueError, match="uint8"):
        fingerprint_quality_summary(keys, keys, fps.astype(np.float32))
    with pytest.raises(ValueError, match="Morgan 2048"):
        fingerprint_quality_summary(keys, keys, fps[:, :255])
