#!/usr/bin/env python3
"""Losslessly archive every performance artifact with deterministic duplicate links."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import io
import json
import shutil
import tarfile
from pathlib import Path


def digest(data):
    return hashlib.sha256(data).hexdigest()


def save(path, value):
    path.write_text(json.dumps(value, indent=2) + "\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--review", type=Path, required=True)
    parser.add_argument("--quiet-window", type=Path, required=True)
    parser.add_argument("--runner-log", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    source = args.input.resolve()
    args.output.mkdir(parents=True, exist_ok=False)
    paths = sorted(p for p in source.rglob("*") if p.is_file())
    manifest = {}
    first = {}
    archive = args.output / "raw_evidence.tar.gz"
    with archive.open("wb") as stream:
        with gzip.GzipFile(fileobj=stream, mode="wb", filename="", mtime=0) as zipped:
            with tarfile.open(
                fileobj=zipped, mode="w|", format=tarfile.PAX_FORMAT
            ) as tar:
                for path in paths:
                    if path.is_symlink():
                        raise ValueError(f"Unexpected symbolic link: {path}")
                    name = str(path.relative_to(source))
                    data = path.read_bytes()
                    mode = path.stat().st_mode & 0o777
                    sha = digest(data)
                    manifest[name] = {"sha256": sha, "bytes": len(data), "mode": mode}
                    info = tarfile.TarInfo(name)
                    info.mode = mode
                    key = (sha, mode)
                    if key in first:
                        info.type = tarfile.LNKTYPE
                        info.linkname = first[key]
                        tar.addfile(info)
                    else:
                        first[key] = name
                        info.size = len(data)
                        tar.addfile(info, io.BytesIO(data))
    save(args.output / "raw_evidence_manifest.json", manifest)
    verified = set()
    with tarfile.open(archive) as tar:
        for member in tar:
            if member.name in verified or member.name not in manifest:
                raise ValueError(f"Unexpected or repeated member: {member.name}")
            data = tar.extractfile(member).read()
            expected = manifest[member.name]
            if (
                len(data) != expected["bytes"]
                or digest(data) != expected["sha256"]
                or member.mode != expected["mode"]
            ):
                raise ValueError(f"Archive differs from original: {member.name}")
            if digest((source / member.name).read_bytes()) != expected["sha256"]:
                raise ValueError(f"Original changed during archive: {member.name}")
            verified.add(member.name)
    if verified != set(manifest):
        raise ValueError("Archive coverage differs from original")
    for path in sorted(args.review.iterdir()):
        if path.is_file():
            shutil.copyfile(path, args.output / path.name)
    shutil.copyfile(args.quiet_window, args.output / "quiet_window.json")
    shutil.copyfile(args.runner_log, args.output / "runner.log")
    save(
        args.output / "archive_verification.json",
        {
            "source_directory": str(source),
            "archiver_sha256": digest(Path(__file__).read_bytes()),
            "files_verified": len(verified),
            "logical_bytes": sum(row["bytes"] for row in manifest.values()),
            "distinct_content_and_mode_entries": len(first),
            "archive_bytes": archive.stat().st_size,
            "archive_sha256": digest(archive.read_bytes()),
            "all_original_bytes_rechecked": True,
            "all_archive_entries_rechecked": True,
            "format": (
                "Deterministic gzip/PAX tar; byte-identical duplicate files use "
                "ordinary hardlink entries. Standard tar extraction restores every "
                "path. No sample or raw-output bytes are discarded."
            ),
        },
    )
    save(
        args.output / "manifest.json",
        {
            p.name: {"sha256": digest(p.read_bytes()), "bytes": p.stat().st_size}
            for p in sorted(args.output.iterdir())
            if p.is_file()
        },
    )
    print(
        json.dumps(
            {
                "files_verified": len(verified),
                "archive_bytes": archive.stat().st_size,
                "archive_sha256": digest(archive.read_bytes()),
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
