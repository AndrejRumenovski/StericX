"""Verify the published portion of the completed, immutable engine campaign.

The full seal also names workstation-only raw measurements under .stericx.
This CI gate verifies the published campaign files and archived Python syntax;
it does not claim to rerun science or verify raw files absent from GitHub.
"""

import ast
import hashlib
import json
import subprocess
import tomllib
from pathlib import Path, PurePosixPath

ROOT = Path(__file__).resolve().parents[1]
CAMPAIGN = PurePosixPath("docs/ultimate_engine")
SEAL_SHA256 = "0495e25e2bba1e3eea01c30f4f28c76b821456fe6f7a3f12354b32389f5cd2a4"


def published_records(data):
    """Select exactly the campaign paths, without requiring private raw files."""
    records = []
    seen = set()
    for row in data["records"]:
        path = PurePosixPath(row["path"])
        if path.is_absolute() or ".." in path.parts:
            raise ValueError(f"Unsafe manifest path: {path}")
        if path in seen:
            raise ValueError(f"Duplicate manifest path: {path}")
        seen.add(path)
        if path.is_relative_to(CAMPAIGN):
            records.append(row)
    if not records:
        raise ValueError("No published campaign records")
    return records


def verify_records(root, records, exclusions):
    """Reject changed files and lint exclusions beyond the sealed Python names."""
    archived_python = set()
    for row in records:
        relative = PurePosixPath(row["path"])
        path = root / relative
        if path.is_symlink():
            raise ValueError(f"Archived file replaced with symlink: {relative}")
        content = path.read_bytes()
        if len(content) != row["bytes"]:
            raise ValueError(f"Size mismatch: {relative}")
        if hashlib.sha256(content).hexdigest() != row["sha256"]:
            raise ValueError(f"Hash mismatch: {relative}")
        if relative.suffix == ".py":
            if relative.parent != CAMPAIGN:
                raise ValueError(f"Unexpected nested archived Python: {relative}")
            ast.parse(content, filename=str(relative))
            archived_python.add(relative.name)
    if len(exclusions) != len(set(exclusions)) or set(exclusions) != archived_python:
        raise ValueError("Ruff exclusions must be exactly the sealed Python filenames")
    return len(archived_python)


def main():
    directory = ROOT / CAMPAIGN
    seal = directory / "final_campaign_evidence.json"
    content = seal.read_bytes()
    if hashlib.sha256(content).hexdigest() != SEAL_SHA256:
        raise ValueError("Original campaign seal changed")
    records = published_records(json.loads(content))
    config = tomllib.loads((directory / ".ruff.toml").read_text())
    if set(config) != {"extend", "extend-exclude"}:
        raise ValueError("Archive Ruff config may only inherit and name sealed files")
    if config["extend"] != "../../pyproject.toml":
        raise ValueError("New campaign helpers must inherit repository Ruff rules")
    archived_python = verify_records(ROOT, records, config["extend-exclude"])
    tracked = set(
        subprocess.check_output(
            ["git", "ls-files", "-z", "--", str(CAMPAIGN)], cwd=ROOT
        )
        .decode()
        .split("\0")
    )
    if any(row["path"] not in tracked for row in records):
        raise ValueError("A sealed published file is missing from the Git index")
    print(
        json.dumps(
            {
                "passed": True,
                "published_files": len(records),
                "archived_python_syntax_checks": archived_python,
                "seal_sha256": SEAL_SHA256,
                "scope": "Published campaign integrity; private raw runs not verified",
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
