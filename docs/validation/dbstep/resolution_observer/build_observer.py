#!/usr/bin/env python3
"""Build a separately targeted observation adapter with immutable provenance."""

import argparse
import difflib
import hashlib
import json
import os
import subprocess
import sys
from pathlib import Path

ADAPTER = Path(__file__).resolve().parent
ROOT = ADAPTER.parents[3]


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline", type=Path, required=True)
    parser.add_argument("--oracle", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--evidence", type=Path, required=True)
    args = parser.parse_args()
    evidence = args.evidence.resolve()
    evidence.mkdir(parents=True, exist_ok=False)
    baseline = json.loads((args.baseline / "manifest.json").read_text())
    baseline_source = Path(baseline["source"])
    source_equivalence = []
    for record in baseline["sources"]:
        relative = Path(record["path"]).relative_to(baseline_source)
        current = ROOT / relative
        assert sha(current) == record["sha256"], relative
        source_equivalence.append({"path": str(relative), "sha256": sha(current)})
    original = ROOT / "docs/scientific_remediation/observer"
    delta = []
    for old in sorted(original.rglob("*")):
        if not old.is_file():
            continue
        relative = old.relative_to(original)
        new = ADAPTER / relative
        if relative != Path("src/main.rs"):
            assert old.read_bytes() == new.read_bytes(), relative
        delta.extend(
            difflib.unified_diff(
                old.read_text().splitlines(keepends=True),
                new.read_text().splitlines(keepends=True),
                fromfile=str(Path("original") / relative),
                tofile=str(Path("resolution") / relative),
            )
        )
    for new in sorted(ADAPTER.rglob("*")):
        if new.is_file() and not (original / new.relative_to(ADAPTER)).exists():
            delta.extend(
                difflib.unified_diff(
                    [],
                    new.read_text().splitlines(keepends=True),
                    fromfile="/dev/null",
                    tofile=str(Path("resolution") / new.relative_to(ADAPTER)),
                )
            )
    diff = evidence / "adapter.diff"
    diff.write_text("".join(delta))
    command = [
        sys.executable,
        str(ROOT / "scripts/replay_scientific_remediation.py"),
        "build",
        "--oracle",
        str(args.oracle.resolve()),
        "--adapter",
        str(ADAPTER),
        "--output",
        str(args.output.resolve()),
    ]
    receipt = {
        "purpose": "SUT observation only; no independent scientific pass inferred",
        "baseline_manifest": {
            "path": str((args.baseline / "manifest.json").resolve()),
            "sha256": sha(args.baseline / "manifest.json"),
        },
        "source_equivalence": source_equivalence,
        "adapter_diff": {"path": str(diff), "sha256": sha(diff)},
        "command": command,
        "controller_sha256": sha(__file__),
    }
    (evidence / "started.json").write_text(json.dumps(receipt, indent=2) + "\n")
    with (
        (evidence / "stdout.txt").open("w") as stdout,
        (evidence / "stderr.txt").open("w") as stderr,
    ):
        result = subprocess.run(
            command,
            cwd=ROOT,
            stdout=stdout,
            stderr=stderr,
            env={**os.environ, "PYTHONDONTWRITEBYTECODE": "1"},
            check=False,
        )
    receipt["returncode"] = result.returncode
    if result.returncode == 0:
        built = json.loads((args.output / "manifest.json").read_text())
        prefix = (ROOT / "src/geometry/buried_volume.rs").read_bytes()
        observed = (args.output / "observer/src/observed_buried_volume.rs").read_bytes()
        assert observed[: len(prefix)] == prefix
        assert (
            built["bv_prefix_sha256"]
            == baseline["bv_prefix_sha256"]
            == sha(ROOT / "src/geometry/buried_volume.rs")
        )
        receipt["prefix_bytes"] = len(prefix)
        receipt["prefix_sha256"] = built["bv_prefix_sha256"]
        receipt["build_manifest_sha256"] = sha(args.output / "manifest.json")
        receipt["observer_binary"] = built["binaries"]["observer"]
    (evidence / "receipt.json").write_text(json.dumps(receipt, indent=2) + "\n")
    print(json.dumps(receipt, indent=2))
    raise SystemExit(result.returncode)


if __name__ == "__main__":
    main()
