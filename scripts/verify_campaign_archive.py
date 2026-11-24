"""Verify the published portions of completed, immutable engine campaigns.

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
ARCHIVES = [
    (CAMPAIGN, "final_campaign_evidence.json", SEAL_SHA256),
    (
        PurePosixPath("docs/occupancy_precision"),
        "published_archive.json",
        "6678d84c94f64f65c7b771413760dc294f25c484a3e2d3be553ec2ff69c384d7",
    ),
    (
        PurePosixPath("docs/experimental_validation"),
        "campaign_manifest.json",
        "7d24348d5a555f29c41ceed641f678fcf7e993f571927aba46c5371df729c6bc",
    ),
]


def published_records(data, campaign=CAMPAIGN):
    """Select exactly the campaign paths, without requiring private raw files."""
    records = []
    seen = set()
    source_records = data["records"]
    if isinstance(source_records, dict):
        source_records = [
            dict(metadata, path=path) for path, metadata in source_records.items()
        ]
    for row in source_records:
        path = PurePosixPath(row["path"])
        if path.is_absolute() or ".." in path.parts:
            raise ValueError(f"Unsafe manifest path: {path}")
        if path in seen:
            raise ValueError(f"Duplicate manifest path: {path}")
        seen.add(path)
        if path.is_relative_to(campaign):
            records.append(row)
    if not records:
        raise ValueError("No published campaign records")
    return records


def verify_records(root, records, exclusions, campaign=CAMPAIGN):
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
            if relative.parent != campaign:
                raise ValueError(f"Unexpected nested archived Python: {relative}")
            ast.parse(content, filename=str(relative))
            archived_python.add(relative.name)
    if len(exclusions) != len(set(exclusions)) or set(exclusions) != archived_python:
        raise ValueError("Ruff exclusions must be exactly the sealed Python filenames")
    return len(archived_python)


def verify_index(records, tracked, exclusions):
    """Require archived files in Git and prevent hiding new nested helpers."""
    archived = {row["path"] for row in records}
    if not archived.issubset(tracked):
        raise ValueError("A sealed published file is missing from the Git index")
    for name in tracked - archived:
        if PurePosixPath(name).name in exclusions:
            raise ValueError(f"Archive exclusion would hide a new helper: {name}")


def verify_archive(root, campaign, seal_name, expected_hash):
    directory = root / campaign
    seal = directory / seal_name
    content = seal.read_bytes()
    if hashlib.sha256(content).hexdigest() != expected_hash:
        raise ValueError("Original campaign seal changed")
    records = published_records(json.loads(content), campaign)
    config = tomllib.loads((directory / ".ruff.toml").read_text())
    if set(config) != {"extend", "extend-exclude"}:
        raise ValueError("Archive Ruff config may only inherit and name sealed files")
    if config["extend"] != "../../pyproject.toml":
        raise ValueError("New campaign helpers must inherit repository Ruff rules")
    archived_python = verify_records(root, records, config["extend-exclude"], campaign)
    tracked = set(
        subprocess.check_output(
            ["git", "ls-files", "-z", "--", str(campaign)], cwd=root
        )
        .decode()
        .split("\0")
    )
    verify_index(records, tracked, config["extend-exclude"])
    return {
        "campaign": str(campaign),
        "passed": True,
        "published_files": len(records),
        "archived_python_syntax_checks": archived_python,
        "seal_sha256": expected_hash,
        "scope": "Published campaign integrity; private raw runs not verified",
    }


def main():
    print(
        json.dumps(
            [verify_archive(ROOT, *archive) for archive in ARCHIVES],
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
