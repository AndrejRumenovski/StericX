"""Publish sealed forensic evidence in small, verified Git archive parts."""

import gzip
import hashlib
import json
import tarfile
import tempfile
from datetime import UTC, datetime
from pathlib import Path

AUDIT = Path(__file__).resolve().parents[1]
REPO = AUDIT.parent.parent
PART_BYTES = 48 * 1024 * 1024


def digest(path):
    result = hashlib.sha256()
    with path.open("rb") as source:
        for block in iter(lambda: source.read(1024 * 1024), b""):
            result.update(block)
    return result.hexdigest()


def record(path):
    return {
        "path": path.relative_to(AUDIT).as_posix(),
        "bytes": path.stat().st_size,
        "sha256": digest(path),
    }


def local_only(name):
    if (
        name == "frozen/kraken_dft_source.py.txt"
        or (name.startswith("conventions/") and name.endswith(".py.txt"))
        or (name.startswith("historical_dependencies/") and name.endswith(".py.body"))
        or name == "historical_dependencies/journal_to_current_methods.diff"
    ):
        return (
            "Downloaded third-party source or source excerpt; retained locally under "
            "the repository's existing source-publication policy. Acquisition URLs "
            "and original hashes remain published."
        )
    return None


def direct(name):
    path = Path(name)
    return (
        len(path.parts) == 1
        or name.startswith("publication/")
        or (name.startswith("scripts/") and path.suffix == ".py")
        or name.startswith("adapter/")
        or name.startswith("historical_target/final/plots/")
        or name.startswith("historical_target/baseline/plots/")
        or name.startswith("family_plots_v2/")
        or name == "convergence/analysis/headline_grid_convergence.png"
        or name == "results/exact_available_input_algorithm_metrics.csv"
        or name == "results/exact_historical_input_subset.json"
        or name == "historical_target/conditional_ceilings.csv"
        or name == "historical_target/stage_metrics.csv"
        or name == "historical_target/final/metrics.csv"
        or name == "historical_target/baseline/metrics.csv"
    )


def main():
    pub = AUDIT / "publication"
    if (pub / "manifest.json").exists():
        raise SystemExit("Refusing to replace an existing publication manifest.")
    seal = AUDIT / "BUNDLE_MANIFEST.json"
    seal_hash = digest(seal)
    original = json.loads(seal.read_text())
    # Validate the complete scientific seal before changing only Git packaging.
    for name, info in original["files"].items():
        path = AUDIT / name
        if path.stat().st_size != info["bytes"] or digest(path) != info["sha256"]:
            raise SystemExit(f"Sealed evidence changed: {name}")
    visible, bundled, withheld = [], [], []
    for path in sorted(AUDIT.rglob("*")):
        rel = path.relative_to(AUDIT)
        if not path.is_file() or {"target", "__pycache__"}.intersection(rel.parts):
            continue
        if path.is_symlink():
            raise SystemExit(f"Unexpected symlink: {rel}")
        name = rel.as_posix()
        reason = local_only(name)
        if reason:
            withheld.append({**record(path), "reason": reason})
        elif direct(name):
            visible.append(path)
        else:
            bundled.append(path)
    print(
        json.dumps(
            {
                "direct": len(visible),
                "bundled": len(bundled),
                "local_only": len(withheld),
                "uncompressed_bytes": sum(p.stat().st_size for p in bundled),
            }
        ),
        flush=True,
    )
    with tempfile.TemporaryDirectory(prefix="stericx-forensics-publication-") as tmp:
        archive = Path(tmp) / "evidence.tar.gz"
        with (
            archive.open("wb") as stream,
            gzip.GzipFile(
                fileobj=stream, mode="wb", filename="", mtime=0, compresslevel=6
            ) as compressed,
            tarfile.open(
                fileobj=compressed, mode="w|", format=tarfile.PAX_FORMAT
            ) as tar,
        ):
            for path in bundled:
                member = tar.gettarinfo(
                    str(path), arcname=path.relative_to(AUDIT).as_posix()
                )
                member.uid = member.gid = 0
                member.uname = member.gname = ""
                member.mtime = 0
                member.mode = 0o755 if path.stat().st_mode & 0o111 else 0o644
                with path.open("rb") as stream:
                    tar.addfile(member, stream)
        archive_hash = digest(archive)
        parts = []
        with archive.open("rb") as source:
            number = 1
            while block := source.read(PART_BYTES):
                path = pub / f"evidence.tar.gz.part{number:03d}"
                with path.open("xb") as output:
                    output.write(block)
                parts.append(record(path))
                number += 1
    manifest = {
        "schema": 1,
        "created_utc": datetime.now(UTC).isoformat(),
        "frozen_commit": original["source_commit"],
        "original_sealed_manifest_sha256": seal_hash,
        "archive_sha256": archive_hash,
        "bundles": parts,
        "direct_files": [record(p) for p in visible],
        "bundled_files": [record(p) for p in bundled],
        "local_only": withheld,
        "policy": (
            "All original scientific bytes and seals are unchanged. Authored reports, "
            "scientific data, outputs and failures are published directly or bundled. "
            "Only explicitly listed third-party source captures remain local, "
            "following the existing scientific-audit publication policy."
        ),
    }
    (pub / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    # Keep the already sealed audit/.gitignore unchanged; apply publication
    # visibility at the repository root instead.
    names = [p.relative_to(AUDIT).as_posix() for p in visible]
    names += [p["path"] for p in parts] + ["publication/manifest.json"]
    prefix = "/docs/kraken_residual_forensics/"
    ignore = REPO / ".gitignore"
    marker = "# Sealed Kraken forensic publication (raw evidence restores from parts)."
    text = ignore.read_text()
    if marker in text:
        raise SystemExit("Publication ignore block already exists.")
    ignore.write_text(
        text.rstrip()
        + "\n"
        + marker
        + "\n"
        + prefix
        + "**\n!"
        + prefix
        + "**/\n"
        + "".join("!" + prefix + name + "\n" for name in sorted(set(names)))
    )
    assert digest(seal) == seal_hash
    print(
        json.dumps(
            {
                "parts": len(parts),
                "compressed_bytes": sum(p["bytes"] for p in parts),
                "seal_unchanged": True,
            }
        ),
        flush=True,
    )


if __name__ == "__main__":
    main()
