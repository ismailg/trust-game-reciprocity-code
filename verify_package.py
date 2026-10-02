#!/usr/bin/env python3
"""Verify the exact code-only distribution without opening private inputs."""
import hashlib
import json
from pathlib import Path, PurePosixPath

ROOT = Path(__file__).resolve().parent
MANIFEST = "RELEASE_MANIFEST.json"


def verify(root=ROOT):
    root = Path(root)
    manifest = json.loads((root / MANIFEST).read_text())
    expected = manifest["files"]
    if not expected or len({item["path"] for item in expected}) != len(expected):
        raise ValueError("Empty or duplicate release file list")
    expected_paths = {item["path"] for item in expected} | {MANIFEST}
    actual = set()
    for path in root.rglob("*"):
        # A normal GitHub clone has its own Git metadata. It is never a
        # distribution file and is not copied into private analysis workspaces.
        if path.relative_to(root).parts[0] == ".git":
            if (root / ".git").is_symlink() or not (root / ".git").is_dir():
                raise ValueError("Git metadata must be a local directory")
            continue
        if path.is_symlink():
            raise ValueError(f"Symbolic link is not permitted: {path.relative_to(root)}")
        if path.is_file():
            actual.add(path.relative_to(root).as_posix())
    if actual != expected_paths:
        raise ValueError(f"File inventory mismatch: missing={sorted(expected_paths-actual)}, "
                         f"unexpected={sorted(actual-expected_paths)}")
    for item in expected:
        name = item["path"]
        relative = PurePosixPath(name)
        if relative.is_absolute() or ".." in relative.parts or "\\" in name:
            raise ValueError("Unsafe manifest path")
        if relative.suffix.lower() in {".rds", ".rdata", ".rda", ".pdf", ".docx", ".zip", ".csv"}:
            if name not in {"data_dictionary.csv", "output_map.csv"}:
                raise ValueError(f"Non-code payload is not permitted: {name}")
        if relative.parts[0] in {".git", "results", "modsData", "output", "private_fits"}:
            raise ValueError(f"Private working material is not permitted: {name}")
        if relative.parts[0] == "Data" and name != "Data/README.md":
            raise ValueError(f"Input data are not permitted: {name}")
        raw = (root / name).read_bytes()
        raw.decode("utf-8")  # Every distributed file is inspectable text.
        if len(raw) != item["bytes"] or hashlib.sha256(raw).hexdigest() != item["sha256"]:
            raise ValueError(f"Content checksum mismatch: {name}")
    return manifest


if __name__ == "__main__":
    result = verify()
    print(f"Verified {len(result['files'])} declared text files plus the manifest; no unexpected files.")
