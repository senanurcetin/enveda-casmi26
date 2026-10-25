import pytest
from packaging.markers import default_environment

from scripts.check_research_lock import lock_entries, check_installed


def test_hash_pins_and_platform_markers_are_enforced(tmp_path):
    path = tmp_path / "lock.txt"
    path.write_text('example==1.2 --hash=sha256:' + 'a' * 64 + '\n'
                    'other==3.0 ; sys_platform == "linux" --hash=sha256:' + 'b' * 64 + '\n')
    entries = lock_entries(path)
    versions = {"example": "1.2"}
    assert check_installed(entries, environment={**default_environment(), "sys_platform": "win32"},
                           version_of=versions.__getitem__) == versions
    with pytest.raises(ValueError, match="differs from lock"):
        check_installed(entries[:1], version_of=lambda _: "9.0")
    for invalid in ('example>=1.2 --hash=sha256:' + 'a' * 64, 'example==1.2',
                    'example==1.2 --hash=sha256:bad'):
        path.write_text(invalid)
        with pytest.raises(ValueError):
            lock_entries(path)
