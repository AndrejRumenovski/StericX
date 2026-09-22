"""Verify the previously sealed audit identity before executing reference code."""

import hashlib
import importlib.util
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
AUDIT = ROOT / "docs/scientific_accuracy_audit"
MANIFEST = AUDIT / "manifest_final.json"
# Independently verified v3 audit seal, preserved before remediation.
MANIFEST_SHA256 = "c31b655641a352402b14dc2d4535261a9de9c9004979c5e571d78e45cac82f03"
REFERENCE = AUDIT / "scripts/geometry_reference.py"


def sha(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def verify(path):
    if sha(MANIFEST) != MANIFEST_SHA256:
        raise ValueError(
            "Original audit manifest no longer matches its pre-remediation seal"
        )
    relative = str(path.resolve().relative_to(AUDIT.resolve()))
    entries = [
        row
        for row in json.loads(MANIFEST.read_text())["files"]
        if row["path"] == relative
    ]
    if (
        len(entries) != 1
        or sha(path) != entries[0]["sha256"]
        or path.stat().st_size != entries[0]["bytes"]
    ):
        raise ValueError(
            f"Reference file differs from immutable audit manifest: {relative}"
        )
    return entries[0]


def load_reference():
    verify(REFERENCE)
    spec = importlib.util.spec_from_file_location(
        "sealed_geometry_reference", REFERENCE
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def refuse_audit_output(path):
    if path.resolve().is_relative_to(AUDIT.resolve()):
        raise ValueError("The original scientific audit is immutable")
