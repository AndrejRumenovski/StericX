"""Restore and verify the forensic publication using standard-library helpers."""

import argparse
import importlib.util
import json
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[3]
SPEC = importlib.util.spec_from_file_location(
    "scientific_publication",
    REPO / "docs/scientific_accuracy_audit/publication/manage.py",
)
HELPER = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(HELPER)


def verify(audit, repo, manifest):
    failures, missing_local = [], []
    checked = 0

    def check(root, row):
        nonlocal checked
        try:
            HELPER.check_file(root, row)
            checked += 1
        except (OSError, ValueError) as error:
            failures.append(str(error))

    try:
        HELPER.check_bundles(audit, manifest)
    except (OSError, ValueError) as error:
        failures.append(str(error))
    published = {}
    for group in ("direct_files", "bundled_files", "local_only"):
        for row in manifest[group]:
            published[row["path"]] = row
            if group == "local_only" and not (audit / row["path"]).exists():
                missing_local.append({"path": row["path"], "reason": row["reason"]})
            else:
                check(audit, row)
    seal = HELPER.safe_path(audit, "BUNDLE_MANIFEST.json")
    if HELPER.digest(seal) != manifest["original_sealed_manifest_sha256"]:
        failures.append("Original scientific seal changed")
    original = json.loads(seal.read_text())
    for name, row in original["files"].items():
        current = published.get(name, {})
        if any(current.get(k) != row[k] for k in ("bytes", "sha256")):
            failures.append(f"Publication differs from sealed evidence: {name}")
    frozen = json.loads((audit / "frozen/manifest.json").read_text())
    if original["source_commit"] != manifest["frozen_commit"]:
        failures.append("Frozen source commit mismatch")
    expected = {
        name: row
        for name, row in frozen["source_files"].items()
        if name.startswith(("src/", "tests/", "scripts/", "studies/"))
        or name in {"Cargo.toml", "Cargo.lock", "pyproject.toml", "uv.lock"}
    }
    for name, checksum in original["source_changes"].items():
        path = HELPER.safe_path(repo, name)
        if not path.is_file() or HELPER.digest(path) != checksum:
            failures.append(f"Accepted source/test differs: {name}")
        expected.pop(name, None)
    for name, row in expected.items():
        check(repo, {"path": name, **row})
    return {
        "passed": not failures,
        "scope": "published evidence and recorded repository scientific source",
        "checked_files": checked,
        "local_only_missing": missing_local,
        "failures": failures,
        "note": (
            "Missing local-only source captures are declared publication exclusions. "
            "This is not a new scientific replay or a claim of recovered "
            "historical inputs."
        ),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=("restore", "verify"))
    parser.add_argument(
        "--audit-dir", type=Path, default=Path(__file__).resolve().parents[1]
    )
    parser.add_argument("--repo-dir", type=Path, default=REPO)
    args = parser.parse_args()
    try:
        audit = args.audit_dir.resolve()
        manifest = HELPER.read_manifest(audit)
        result = (
            HELPER.restore(audit, manifest)
            if args.command == "restore"
            else verify(audit, args.repo_dir.resolve(), manifest)
        )
    except (OSError, ValueError, KeyError, TypeError, HELPER.tarfile.TarError) as error:
        result = {"passed": False, "failures": [str(error)]}
    print(json.dumps(result, indent=2))
    return 0 if result["passed"] else 1


if __name__ == "__main__":
    sys.exit(main())
