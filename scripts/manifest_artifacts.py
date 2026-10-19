"""Create or verify a metadata-only input-inventory manifest for CASMI research.

Examples:
  python scripts/manifest_artifacts.py --root /data/casmi/artifacts \
    --file store/meta.parquet --file store/mz.f32 --file store/it.f32 \
    --file store/offsets.i64 --file pool.parquet \
    --output runs/input_snapshot.json

  python scripts/manifest_artifacts.py --root /data/casmi/artifacts \
    --verify runs/input_snapshot.json

The manifest contains only relative file names, byte sizes, SHA-256 hashes,
and a composite inventory digest. It does not include dataset bytes,
absolute local paths, private competition answers or any source URL.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from casmi.validate.source_manifest import create_manifest, verify_manifest  # noqa: E402


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True,
                        help="Parent directory of local CASMI input artifacts")
    action = parser.add_mutually_exclusive_group(required=True)
    action.add_argument("--verify", type=Path, help="Verify an existing manifest")
    action.add_argument("--output", type=Path, help="Create a manifest (exclusive write)")
    parser.add_argument("--file", action="append", default=[],
                        help="Relative file to hash; repeat for each input artifact")
    args = parser.parse_args()
    try:
        if args.verify:
            if args.file:
                parser.error("--file is only valid with --output")
            manifest = json.loads(args.verify.read_text(encoding="utf-8"))
            result = verify_manifest(args.root, manifest)
        else:
            manifest = create_manifest(args.root, args.file)
            args.output.parent.mkdir(parents=True, exist_ok=True)
            with args.output.open("x", encoding="utf-8") as file:
                json.dump(manifest, file, indent=2)
                file.write("\n")
            result = {
                "status": "MANIFEST_CREATED",
                "files_checked": len(manifest["files"]),
                "files_digest_sha256": manifest["files_digest_sha256"],
            }
        print(json.dumps(result, indent=2))
    except (ValueError, OSError, TypeError, KeyError, json.JSONDecodeError) as exc:
        parser.error(str(exc))


if __name__ == "__main__":
    main()
