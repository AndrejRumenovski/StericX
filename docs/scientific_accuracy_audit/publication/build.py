"""Package the preserved audit for Git without changing its original evidence.

Downloaded publisher documents and selected third-party code remain local.
Scientific data, frozen native outputs and derived results are compressed
losslessly. Run only after all audit/publication files have been finalized.
"""

import gzip
import hashlib
import json
import tarfile
import tempfile
from datetime import UTC, datetime
from pathlib import Path

A = Path(__file__).resolve().parents[1]
PUB = A / "publication"
PART_BYTES = 48 * 1024 * 1024


def sha(p):
    h = hashlib.sha256()
    with p.open("rb") as f:
        for b in iter(lambda: f.read(1048576), b""):
            h.update(b)
    return h.hexdigest()


def local_only(rel):
    p = Path(rel)
    if (
        rel.startswith("frozen/reference_sources/")
        or rel == "geometry/alignment/glam_0.30.10_sse2_quat.rs"
    ):
        return (
            "Downloaded third-party source; use the pinned package/dependency "
            "and verify the original source hash."
        )
    if rel.startswith("geometry/sources/") and p.suffix == ".body":
        return (
            "Captured third-party publication/documentation response; "
            "original URL and hash remain in the source manifest."
        )
    if rel.startswith("kinetics/sources/") and p.name not in {
        "manifest.json",
        "nist_2022_constants.txt",
    }:
        return (
            "Captured third-party source, rendered page or transcription; "
            "original URL and hash remain in the source manifest."
        )
    if rel.startswith("models/sources/") and (
        p.suffix in {".pdf", ".html", ".txt", ".ipynb", ".zip"}
        or p.name.endswith("_notebook.py")
        or "/rendered/" in rel
    ):
        return (
            "Downloaded publication, documentation, notebook or supporting archive; "
            "reacquire locally from recorded provenance."
        )
    if rel.startswith("kraken/raw_sources/") and (
        p.name
        in {
            "kraken_paper_pdf.body",
            "kraken_si_pdf.body",
            "kraken_si_zip.body",
            "descriptors.xlsx",
        }
        or p.suffix in {".txt", ".png"}
        or (
            p.name.startswith("official_")
            and p.suffix == ".body"
            and p.name not in {"official_tree.body", "official_tree_master.body"}
        )
    ):
        return (
            "Downloaded paper/supporting information or third-party source; "
            "metadata and derived numerical evidence are published."
        )
    if (
        rel.startswith("kraken/energy_source_search/")
        and p.suffix == ".body"
        and (p.name.startswith("framework_") or "readme" in p.name.lower())
    ):
        return (
            "Downloaded third-party code/documentation; "
            "retrieval metadata and investigation results are published."
        )
    return None


def direct(rel):
    p = Path(rel)
    if rel.startswith("frozen/"):
        return False
    return (
        len(p.parts) == 1
        or p.suffix in {".md", ".py", ".rs", ".png", ".pdf"}
        or rel.startswith("observer/")
        or p.name in {"metrics.json", "expanded_matched_default_metrics.json"}
    )


def info(p):
    return {
        "path": p.relative_to(A).as_posix(),
        "bytes": p.stat().st_size,
        "sha256": sha(p),
    }


def main():
    if (PUB / "manifest.json").exists():
        raise SystemExit(
            "Refusing to replace an existing publication; "
            "preserve its manifest before rebuilding."
        )
    files = []
    for p in A.rglob("*"):
        if not p.is_file() or "__pycache__" in p.parts or p.suffix in {".pyc", ".pyo"}:
            continue
        rel = p.relative_to(A).as_posix()
        if rel.startswith("publication/") and p.name not in {"build.py", "manage.py"}:
            continue
        if rel == ".gitignore":
            continue
        if p.is_symlink():
            raise SystemExit("Unexpected symlink: " + rel)
        files.append(p)
    visible, bundled, withheld = [], [], []
    for p in sorted(files):
        rel = p.relative_to(A).as_posix()
        reason = local_only(rel)
        if reason:
            withheld.append({**info(p), "reason": reason})
        elif direct(rel):
            visible.append(p)
        else:
            bundled.append(p)
    initial = json.loads((A / "manifest_initial.json").read_text())
    original_manifest_hash = sha(A / "manifest_final.json")
    print(
        json.dumps(
            {
                "direct_files": len(visible),
                "bundled_files": len(bundled),
                "local_only_files": len(withheld),
                "uncompressed_bundle_bytes": sum(p.stat().st_size for p in bundled),
            }
        ),
        flush=True,
    )
    records = [info(p) for p in bundled]
    with tempfile.TemporaryDirectory(prefix="stericx-audit-package-") as tmp:
        archive = Path(tmp) / "evidence.tar.gz"
        with (
            archive.open("wb") as out,
            gzip.GzipFile(
                fileobj=out, mode="wb", compresslevel=6, mtime=0, filename=""
            ) as compressed,
            tarfile.open(
                fileobj=compressed, mode="w|", format=tarfile.PAX_FORMAT
            ) as tar,
        ):
            for p in bundled:
                member = tar.gettarinfo(str(p), arcname=p.relative_to(A).as_posix())
                member.uid = member.gid = 0
                member.uname = member.gname = ""
                member.mtime = 0
                member.mode = 0o755 if p.stat().st_mode & 0o111 else 0o644
                with p.open("rb") as f:
                    tar.addfile(member, f)
        archive_hash = sha(archive)
        parts = []
        with archive.open("rb") as src:
            number = 1
            while block := src.read(PART_BYTES):
                target = PUB / f"evidence.tar.gz.part{number:03d}"
                with target.open("xb") as dst:
                    dst.write(block)
                parts.append(info(target))
                number += 1
    # Local raw data remain intact and ignored. Only the explicit publication
    # surface is staged; no blanket force-add of downloaded captures is needed.
    names = [p.relative_to(A).as_posix() for p in visible]
    names += [r["path"] for r in parts] + ["publication/manifest.json", ".gitignore"]
    ignore = A / ".gitignore"
    ignore.write_text(
        "# Restored audit evidence and downloaded source captures remain local.\n"
        "# Publication files are explicitly versioned; see README.md.\n*\n!*/\n"
        + "".join("!/" + n + "\n" for n in sorted(set(names)))
    )
    visible.append(ignore)
    manifest = {
        "schema": 1,
        "created_utc": datetime.now(UTC).isoformat(),
        "frozen_commit": initial["git_commit"],
        "original_sealed_manifest_sha256": original_manifest_hash,
        "policy": (
            "Original evidence is unchanged. This Git publication contains authored "
            "reports/harnesses, numeric/geometry inputs, metadata, "
            "frozen native outputs and derived results. "
            "Publisher/source captures explicitly listed under "
            "local_only are preserved locally and omitted from Git under "
            "repository policy."
        ),
        "archive_format": (
            "Concatenate bundles in listed order to reconstruct "
            "one gzip-compressed tar; paths relative to the audit root."
        ),
        "archive_sha256": archive_hash,
        "bundles": parts,
        "direct_files": [info(p) for p in sorted(visible)],
        "bundled_files": records,
        "local_only": withheld,
    }
    (PUB / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    assert sha(A / "manifest_final.json") == original_manifest_hash
    print(
        json.dumps(
            {
                "bundles": len(parts),
                "compressed_bytes": sum(r["bytes"] for r in parts),
                "direct_files": len(visible),
                "bundled_files": len(records),
                "local_only_files": len(withheld),
                "original_seal_unchanged": True,
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
