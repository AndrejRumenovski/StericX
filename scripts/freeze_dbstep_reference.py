#!/usr/bin/env python3
"""Seal an unmodified official DBSTEP checkout and installed reference runtime."""

import argparse
import hashlib
import importlib.metadata
import json
import platform
import subprocess
import sys
from datetime import UTC, datetime
from pathlib import Path

PIN = "4365353e0a2f24d05336e619759e350bc630a3cb"


def sha(path):
    with Path(path).open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--checkout", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    checkout, output = args.checkout.resolve(), args.output.resolve()
    if (output / "environment.json").exists() or (
        output / "requirements.lock"
    ).exists():
        raise FileExistsError("reference seal already exists")

    def git(*words):
        return subprocess.check_output(
            ["git", "-C", str(checkout), *words], text=True
        ).strip()

    if git("rev-parse", "HEAD") != PIN or git("status", "--porcelain"):
        raise ValueError("reference checkout must be the clean pinned official source")
    import dbstep

    installed = Path(dbstep.__file__).parent
    packages = sorted(
        (dist.metadata["Name"], dist.version)
        for dist in importlib.metadata.distributions()
    )
    source_hashes, installed_hashes = {}, {}
    for name in git("ls-files").splitlines():
        path = checkout / name
        if path.is_file():
            source_hashes[name] = sha(path)
        if name.startswith("dbstep/") and name.endswith(".py"):
            target = installed / Path(name).relative_to("dbstep")
            installed_hashes[str(target)] = sha(target)
            if installed_hashes[str(target)] != source_hashes[name]:
                raise ValueError(
                    f"installed DBSTEP differs from official source: {name}"
                )
    distribution = importlib.metadata.distribution("dbstep")
    record = {
        "created_utc": datetime.now(UTC).isoformat(),
        "command": sys.argv,
        "script_sha256": sha(__file__),
        "official_repository": "https://github.com/patonlab/DBSTEP",
        "commit": PIN,
        "version": dbstep.__version__,
        "checkout": str(checkout),
        "checkout_clean": True,
        "python": sys.version,
        "python_executable": sys.executable,
        "python_sha256": sha(sys.executable),
        "platform": platform.platform(),
        "packages": dict(packages),
        "source_hashes": source_hashes,
        "installed_dbstep_python_hashes": installed_hashes,
        "installed_source_matches_official": True,
        "installation_direct_url": json.loads(
            distribution.read_text("direct_url.json") or "null"
        ),
    }
    output.mkdir(parents=True, exist_ok=True)
    (output / "environment.json").write_text(json.dumps(record, indent=2) + "\n")
    lock = [
        f"# Python {platform.python_version()}; source and runtime receipt: environment.json"  # noqa: E501
    ]
    for name, version in packages:
        lock.append(
            f"dbstep @ git+https://github.com/patonlab/DBSTEP@{PIN}"
            if name.lower() == "dbstep"
            else f"{name}=={version}"
        )
    (output / "requirements.lock").write_text("\n".join(lock) + "\n")
    print(
        json.dumps(
            {
                "commit": PIN,
                "version": dbstep.__version__,
                "packages": len(packages),
                "installed_source_matches_official": True,
                "tracked_reference_files": len(source_hashes),
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
