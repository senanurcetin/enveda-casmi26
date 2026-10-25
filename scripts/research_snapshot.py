"""Read embedded executed source as data; never execute a historical notebook."""
from __future__ import annotations

import ast
import hashlib
import json
from pathlib import Path, PurePosixPath


def embedded_sources(notebook_path):
    notebook = json.loads(Path(notebook_path).read_text(encoding="utf-8"))
    found = []
    for cell in notebook["cells"]:
        if cell["cell_type"] != "code":
            continue
        for node in ast.parse("".join(cell["source"])).body:
            if not isinstance(node, ast.Assign) or not any(isinstance(t, ast.Name) and t.id == "sources" for t in node.targets):
                continue
            call = node.value
            if (not isinstance(call, ast.Call) or not isinstance(call.func, ast.Attribute)
                    or not isinstance(call.func.value, ast.Name) or call.func.value.id != "json"
                    or call.func.attr != "loads" or len(call.args) != 1):
                raise ValueError("Historical source must be a literal JSON assignment")
            found.append(json.loads(ast.literal_eval(call.args[0])))
    if len(found) != 1 or not isinstance(found[0], dict) or not found[0]:
        raise ValueError("Exactly one embedded source inventory is required")
    for name, source in found[0].items():
        path = PurePosixPath(name)
        if (not isinstance(source, str) or path.is_absolute() or ".." in path.parts
                or "\\" in name or path.as_posix() != name or not name.endswith(".py")):
            raise ValueError("Invalid embedded source path or content")
    return found[0]


def snapshot_hashes(path):
    return {name: hashlib.sha256(source.encode()).hexdigest() for name, source in embedded_sources(path).items()}
