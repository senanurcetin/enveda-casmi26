"""Check repository-local Markdown and dashboard destinations without a browser."""
from __future__ import annotations

import re
import subprocess
from html.parser import HTMLParser
from pathlib import Path
from urllib.parse import unquote, urlsplit

ROOT = Path(__file__).resolve().parents[1]


class HTMLLinks(HTMLParser):
    def __init__(self):
        super().__init__()
        self.links = []

    def handle_starttag(self, tag, attrs):
        self.links.extend(value for key, value in attrs if key in ("href", "src") and value)


def main():
    names = subprocess.check_output(
        ["git", "ls-files", "--cached", "--others", "--exclude-standard", "--", "*.md", "site/*.html"],
        cwd=ROOT, text=True,
    ).splitlines()
    checked, inspected, failures = 0, 0, []
    for name in sorted(set(names)):
        path = ROOT / name
        if not path.is_file():
            continue  # Explicit removals are still in the index before staging.
        inspected += 1
        content = path.read_text(encoding="utf-8")
        if path.suffix == ".md":
            content = re.sub(r"^```.*?^```[^\n]*", "", content, flags=re.S | re.M)
        parser = HTMLLinks()
        parser.feed(content)
        links = parser.links + re.findall(r"!?\[[^\]\n]*\]\(([^)\n]+)\)", content)
        for destination in links:
            destination = destination.split(' "', 1)[0].strip("<>")
            url = urlsplit(destination)
            if url.scheme or url.netloc or not url.path:
                continue
            target = (path.parent / unquote(url.path)).resolve()
            checked += 1
            if not target.is_relative_to(ROOT) or not target.exists():
                failures.append(f"{name}: {destination}")
    if failures:
        raise SystemExit("Invalid documentation destinations:\n" + "\n".join(failures))
    print(f"Validated {checked} local destinations across {inspected} documentation/dashboard files; remote URLs and anchors not checked")


if __name__ == "__main__":
    main()
