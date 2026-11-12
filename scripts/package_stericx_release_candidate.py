#!/usr/bin/env python3
"""Copy admitted RC artifacts into a local, checksummed, unpublished bundle."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import io
import json
import shutil
import tarfile
from pathlib import Path

import validate_dbstep_rc_local as rc
from review_dbstep_resolution_release import RC_COMMIT, WINDOWS_WHEEL_SHA, require


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ["review", "versions", "output", "archive"]:
        parser.add_argument("--" + name, type=Path, required=True)
    args = parser.parse_args()
    review = json.loads(args.review.read_text())
    require(
        review["release_candidate_admitted"]
        and review["rc_source_commit"] == RC_COMMIT
        and not review["publication_performed"],
        "Release candidate has not passed admission",
    )
    for record in review["evidence_rechecked"]:
        require(rc.sha(record["path"]) == record["sha256"], "Admitted evidence changed")
    require(
        rc.record(args.versions) in review["evidence_rechecked"],
        "Artifact inventory is not bound to the release review",
    )
    versions = json.loads(args.versions.read_text())
    selected = [
        row
        for row in versions["artifacts"]
        if row["kind"] == "sdist" or row.get("tested_python") == "3.11"
    ]
    require(
        len(selected) == 4
        and {row.get("wheel_platform_tag") for row in selected}
        == {None, "manylinux_2_28_x86_64", "macosx_11_0_arm64", "win_amd64"},
        "Expected one admitted wheel per platform and one source archive",
    )
    groups = versions["platform_payload_groups"]
    require(
        len(groups) == 3
        and {g["wheel_platform_tag"] for g in groups}
        == {"manylinux_2_28_x86_64", "macosx_11_0_arm64", "win_amd64"},
        "Platform payload coverage incomplete",
    )
    for group in groups:
        tag = group["wheel_platform_tag"]
        rows = [a for a in versions["artifacts"] if a.get("wheel_platform_tag") == tag]
        require(
            len(rows) == 4
            and len(group["artifacts"]) == 4
            and set(group["artifacts"]) == {a["path"] for a in rows}
            and sorted(a["tested_python"] for a in rows)
            == ["3.11", "3.12", "3.13", "3.14"],
            "Payload group does not match the tested artifacts",
        )
        if tag != "win_amd64":
            require(
                group["identical_native_bytes"]
                and group[
                    "identical_runtime_payload_after_python_newline_normalization"
                ]
                and sorted(group["tested_python_versions"])
                == ["3.11", "3.12", "3.13", "3.14"]
                and all(
                    a["runtime_entry_hashes"] == rows[0]["runtime_entry_hashes"]
                    and a["native_entry_hashes"] == rows[0]["native_entry_hashes"]
                    for a in rows
                ),
                "Selected runtime payload lacks all-interpreter validation",
            )
    require(
        next(a for a in selected if a.get("wheel_platform_tag") == "win_amd64")[
            "sha256"
        ]
        == WINDOWS_WHEEL_SHA,
        "Selected Windows wheel differs from the single-artifact test",
    )
    require(review["selected_windows_wheel_jobs"] == 4, "Windows artifact gate missing")
    args.output.mkdir(parents=True, exist_ok=False)
    args.archive.mkdir(parents=True, exist_ok=False)
    artifacts = {}
    for row in selected:
        source = Path(row["absolute_path"])
        require(rc.sha(source) == row["sha256"], "Selected artifact changed")
        destination = args.output / source.name
        shutil.copyfile(source, destination)
        require(rc.sha(destination) == row["sha256"], "Artifact copy differs")
        artifacts[source.name] = {
            "sha256": row["sha256"],
            "bytes": destination.stat().st_size,
            "source_artifact": row["path"],
        }
    receipt = {
        "version": "0.4.0rc1",
        "source_commit": RC_COMMIT,
        "release_review": rc.record(args.review),
        "controller": rc.record(__file__),
        "artifacts": artifacts,
        "publication_performed": False,
    }
    (args.output / "manifest.json").write_text(json.dumps(receipt, indent=2) + "\n")
    (args.output / "SHA256SUMS").write_text(
        "".join(f"{r['sha256']}  {name}\n" for name, r in sorted(artifacts.items()))
    )
    (args.output / "README.txt").write_text(
        "StericX Python 0.4.0rc1; Rust engine 0.3.0. Unpublished release candidate.\n"
        "CPython 3.11-3.14: Linux x86-64 (glibc >=2.28), Windows x86-64, macOS arm64.\n"
        "Install the wheel for your platform: python -m pip install ./<wheel>.whl\n"
        "The source archive requires a Rust build toolchain.\n"
        "manifest.json binds exact tested artifacts and the release review.\n"
        "See docs/validation/dbstep/FINAL_RESOLUTION.md for scientific scope,\n"
        "conventions, finite-grid accuracy limits, and preserved original failures.\n"
        "Publication to PyPI or GitHub Releases requires separate user approval.\n"
    )
    archive = args.archive / "candidate_artifacts.tar.gz"
    members = {}
    with archive.open("wb") as stream:
        with gzip.GzipFile(fileobj=stream, mode="wb", filename="", mtime=0) as zipped:
            with tarfile.open(
                fileobj=zipped, mode="w|", format=tarfile.PAX_FORMAT
            ) as tar:
                for path in sorted(args.output.iterdir()):
                    data = path.read_bytes()
                    members[path.name] = hashlib.sha256(data).hexdigest()
                    info = tarfile.TarInfo(path.name)
                    info.size = len(data)
                    info.mode = 0o644
                    tar.addfile(info, io.BytesIO(data))
    with tarfile.open(archive) as tar:
        require(
            {m.name for m in tar} == set(members), "Bundle archive coverage changed"
        )
        for member in tar.getmembers():
            require(
                hashlib.sha256(tar.extractfile(member).read()).hexdigest()
                == members[member.name],
                "Bundle archive member changed",
            )
    for name in ["manifest.json", "SHA256SUMS", "README.txt"]:
        shutil.copyfile(args.output / name, args.archive / name)
    (args.archive / "archive_verification.json").write_text(
        json.dumps(
            {"archive": rc.record(archive), "verified_members": members}, indent=2
        )
        + "\n"
    )
    print(json.dumps({"artifacts": len(artifacts), "archive": rc.record(archive)}))


if __name__ == "__main__":
    main()
