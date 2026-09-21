"""Restore and verify the published scientific audit using only the standard library."""

import argparse
import hashlib
import io
import json
import os
import re
import sys
import tarfile
import tempfile
from pathlib import Path, PurePosixPath

BLOCK = 1024 * 1024


def relative_path(value):
    if not isinstance(value, str) or not value or "\\" in value or "\0" in value:
        raise ValueError(f"Invalid relative path: {value!r}")
    path = PurePosixPath(value)
    if (
        path.is_absolute()
        or ".." in path.parts
        or path.as_posix() != value
        or value == "."
    ):
        raise ValueError(f"Unsafe relative path: {value!r}")
    return path


def safe_path(root, value):
    path = root
    for part in relative_path(value).parts:
        path /= part
        if path.is_symlink():
            raise ValueError(f"Symlinks are not allowed: {path}")
    return path


def digest(path):
    result = hashlib.sha256()
    with path.open("rb") as source:
        for chunk in iter(lambda: source.read(BLOCK), b""):
            result.update(chunk)
    return result.hexdigest()


def check_file(root, row):
    path = safe_path(root, row["path"])
    if not path.is_file():
        raise ValueError(f"Missing regular file: {row['path']}")
    if path.stat().st_size != row["bytes"] or digest(path) != row["sha256"]:
        raise ValueError(f"Size or SHA-256 mismatch: {row['path']}")
    return path


def read_manifest(audit):
    manifest = json.loads(safe_path(audit, "publication/manifest.json").read_text())
    if manifest.get("schema") != 1:
        raise ValueError("Unsupported publication manifest schema")
    for key in ("archive_sha256", "original_sealed_manifest_sha256"):
        if not re.fullmatch(r"[0-9a-f]{64}", manifest[key]):
            raise ValueError(f"Invalid SHA-256: {key}")
    seen = set()
    for group in ("direct_files", "bundled_files", "local_only", "bundles"):
        for row in manifest[group]:
            relative_path(row["path"])
            if row["path"] in seen:
                raise ValueError(f"Duplicate manifest path: {row['path']}")
            seen.add(row["path"])
            if type(row["bytes"]) is not int or row["bytes"] < 0:
                raise ValueError(f"Invalid byte count: {row['path']}")
            if not re.fullmatch(r"[0-9a-f]{64}", row["sha256"]):
                raise ValueError(f"Invalid SHA-256: {row['path']}")
    return manifest


def check_bundles(audit, manifest):
    combined = hashlib.sha256()
    paths = []
    for row in manifest["bundles"]:
        path = safe_path(audit, row["path"])
        part_hash = hashlib.sha256()
        length = 0
        with path.open("rb") as source:
            for chunk in iter(lambda: source.read(BLOCK), b""):
                combined.update(chunk)
                part_hash.update(chunk)
                length += len(chunk)
        if length != row["bytes"] or part_hash.hexdigest() != row["sha256"]:
            raise ValueError(f"Bundle size or SHA-256 mismatch: {row['path']}")
        paths.append(path)
    if combined.hexdigest() != manifest["archive_sha256"]:
        raise ValueError("Combined archive SHA-256 mismatch")
    return paths


class BundleStream(io.RawIOBase):
    """Read ordered split archive parts without joining them on disk or in memory."""

    def __init__(self, paths):
        super().__init__()
        self.paths = iter(paths)
        self.current = None

    def readable(self):
        return True

    def readinto(self, buffer):
        while True:
            if self.current is None:
                path = next(self.paths, None)
                if path is None:
                    return 0
                self.current = path.open("rb")
            count = self.current.readinto(buffer)
            if count:
                return count
            self.current.close()
            self.current = None

    def close(self):
        if self.current is not None:
            self.current.close()
        super().close()


def restore(audit, manifest):
    paths = check_bundles(audit, manifest)
    expected = {row["path"]: row for row in manifest["bundled_files"]}
    # Reject an existing mismatch before restoring any files.
    for row in expected.values():
        path = safe_path(audit, row["path"])
        if path.exists():
            check_file(audit, row)
    seen = set()
    restored = 0
    with BundleStream(paths) as raw, io.BufferedReader(raw) as stream:
        with tarfile.open(fileobj=stream, mode="r|gz") as archive:
            for member in archive:
                relative_path(member.name)
                if member.type not in (tarfile.REGTYPE, tarfile.AREGTYPE):
                    raise ValueError(
                        f"Archive member is not a regular file: {member.name}"
                    )
                if member.name not in expected or member.name in seen:
                    raise ValueError(
                        f"Unlisted or duplicate archive member: {member.name}"
                    )
                seen.add(member.name)
                row = expected[member.name]
                if member.size != row["bytes"]:
                    raise ValueError(f"Archive member size mismatch: {member.name}")
                destination = safe_path(audit, member.name)
                destination.parent.mkdir(parents=True, exist_ok=True)
                temporary = None
                try:
                    with archive.extractfile(member) as source:
                        # Even an already-restored file must match the archived bytes.
                        with tempfile.NamedTemporaryFile(
                            dir=destination.parent,
                            prefix=".audit-restore-",
                            delete=False,
                        ) as output:
                            temporary = Path(output.name)
                            checksum = hashlib.sha256()
                            length = 0
                            for chunk in iter(lambda: source.read(BLOCK), b""):
                                output.write(chunk)
                                checksum.update(chunk)
                                length += len(chunk)
                    if length != row["bytes"] or checksum.hexdigest() != row["sha256"]:
                        raise ValueError(
                            f"Archive member SHA-256 mismatch: {member.name}"
                        )
                    os.chmod(temporary, member.mode & 0o777)
                    # Installing a hard link is atomic and cannot replace a file.
                    safe_path(audit, member.name)
                    try:
                        os.link(temporary, destination)
                        restored += 1
                    except FileExistsError:
                        check_file(audit, row)
                finally:
                    if temporary is not None:
                        temporary.unlink(missing_ok=True)
    missing = sorted(set(expected) - seen)
    if missing:
        raise ValueError(
            f"Archive is missing {len(missing)} listed files: {missing[:5]}"
        )
    return {
        "command": "restore",
        "scope": "published evidence",
        "passed": True,
        "restored": restored,
        "already_present": len(seen) - restored,
        "next_step": (
            "Run verify to check the restored evidence and current repository files."
        ),
    }


def verify(audit, repo, manifest):
    failures = []
    missing_local = []
    checked = 0

    def check(root, row):
        nonlocal checked
        try:
            check_file(root, row)
            checked += 1
        except (OSError, ValueError) as error:
            failures.append(str(error))

    try:
        check_bundles(audit, manifest)
    except (OSError, ValueError) as error:
        failures.append(str(error))
    for group in ("direct_files", "bundled_files"):
        for row in manifest[group]:
            check(audit, row)
    for row in manifest["local_only"]:
        try:
            path = safe_path(audit, row["path"])
            if path.exists():
                check(audit, row)
            else:
                missing_local.append(
                    {
                        "path": row["path"],
                        "reason": row.get("reason", "local-only source capture"),
                    }
                )
        except ValueError as error:
            failures.append(str(error))

    repository_checked = 0
    try:
        final_path = safe_path(audit, "manifest_final.json")
        if digest(final_path) != manifest["original_sealed_manifest_sha256"]:
            raise ValueError("Original sealed manifest SHA-256 mismatch")
        original = json.loads(final_path.read_text())
        published = {
            row["path"]: row
            for group in ("direct_files", "bundled_files", "local_only")
            for row in manifest[group]
        }
        for row in original["files"]:
            current = published.get(row["path"], {})
            if any(current.get(key) != row[key] for key in ("bytes", "sha256")):
                failures.append(
                    f"Publication differs from original sealed inventory: {row['path']}"
                )
        initial = json.loads(safe_path(audit, "manifest_initial.json").read_text())
        if (
            initial["git_commit"] != manifest["frozen_commit"]
            or original["initial_commit"] != manifest["frozen_commit"]
        ):
            raise ValueError("Frozen commit metadata differs between manifests")
        for name, row in initial["repository_files"].items():
            check(repo, {**row, "path": name})
            repository_checked += 1
    except (OSError, ValueError, KeyError, TypeError) as error:
        failures.append(f"Original manifest or repository verification failed: {error}")
    return {
        "command": "verify",
        "scope": "published evidence",
        "passed": not failures,
        "frozen_commit": manifest["frozen_commit"],
        "files_checked_successfully": checked,
        "repository_files_checked": repository_checked,
        "local_only_missing_count": len(missing_local),
        "local_only_missing": missing_local,
        "local_only_present_checked": len(manifest["local_only"]) - len(missing_local),
        "scope_note": (
            "Missing local-only source captures are intentional "
            "publication exclusions. "
            "This does not certify the full original source corpus "
            "or repeat the scientific analyses."
        ),
        "failures": failures,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=("restore", "verify"))
    parser.add_argument(
        "--audit-dir", type=Path, default=Path(__file__).resolve().parents[1]
    )
    parser.add_argument(
        "--repo-dir",
        type=Path,
        help="Repository to check (default: two levels above the audit directory)",
    )
    args = parser.parse_args()
    audit = args.audit_dir.resolve()
    repo = args.repo_dir.resolve() if args.repo_dir else audit.parent.parent
    try:
        manifest = read_manifest(audit)
        result = (
            restore(audit, manifest)
            if args.command == "restore"
            else verify(audit, repo, manifest)
        )
    except (
        OSError,
        ValueError,
        KeyError,
        TypeError,
        tarfile.TarError,
        EOFError,
    ) as error:
        result = {
            "command": args.command,
            "scope": "published evidence",
            "passed": False,
            "failures": [str(error)],
        }
    print(json.dumps(result, indent=2))
    return 0 if result["passed"] else 1


if __name__ == "__main__":
    sys.exit(main())
