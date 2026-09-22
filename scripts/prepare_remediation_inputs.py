#!/usr/bin/env python3
"""Create a new Ni-hDA input bundle with corrected response-temperature metadata.

Preserve the historical CSV, targets, geometries and supplied populations. This
does not rerun conformer searches or reinterpret populations at a new temperature.
The output must be a new directory; native descriptors must be recomputed later.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
AUDIT = REPO / "docs/scientific_accuracy_audit"
# Exact inventory of the immutable audit which authorizes this migration.
AUDIT_MANIFEST_SHA256 = (
    "c31b655641a352402b14dc2d4535261a9de9c9004979c5e571d78e45cac82f03"
)
EXPECTED_IDS = {
    "401",
    "498",
    "723",
    "724",
    "785",
    "1057",
    "1058",
    "2062",
    "2063",
    "2064",
    "2067",
}


def record(path: Path) -> dict:
    data = path.read_bytes()
    return {
        "path": str(path.resolve()),
        "bytes": len(data),
        "sha256": hashlib.sha256(data).hexdigest(),
    }


class SealedAudit:
    """Bind every historical input to the original audit, then recheck at exit."""

    def __init__(self) -> None:
        self.manifest = record(AUDIT / "manifest_final.json")
        if self.manifest["sha256"] != AUDIT_MANIFEST_SHA256:
            raise ValueError("Original audit inventory identity changed")
        inventory = json.loads(Path(self.manifest["path"]).read_text())["files"]
        self.files = {entry["path"]: entry for entry in inventory}
        if len(self.files) != len(inventory):
            raise ValueError("Duplicate paths in original audit inventory")
        self.used: dict[str, dict] = {}

    def check(self, relative: str) -> dict:
        expected = self.files.get(relative)
        if expected is None:
            raise ValueError(f"Not an immutable audited input: {relative}")
        actual = record(AUDIT / relative)
        if any(actual[key] != expected[key] for key in ("sha256", "bytes")):
            raise ValueError(f"Immutable audited input changed: {relative}")
        self.used[relative] = actual
        return actual

    def repository_file(self, path: Path) -> dict:
        resolved = path.resolve(strict=True)
        if not resolved.is_relative_to(REPO):
            raise ValueError(f"Input resolves outside the audited repository: {path}")
        relative = "frozen/repository/" + str(resolved.relative_to(REPO))
        sealed = self.check(relative)
        current = record(resolved)
        if any(current[key] != sealed[key] for key in ("sha256", "bytes")):
            raise ValueError(f"Current input differs from original audit: {resolved}")
        return current

    def recheck(self, live: list[dict]) -> None:
        for expected in [self.manifest, *self.used.values(), *live]:
            if record(Path(expected["path"])) != expected:
                raise ValueError(f"Input changed during migration: {expected['path']}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    sealed = SealedAudit()
    helper_identity = record(Path(__file__))
    source = REPO / "data/reactions_raw.csv"
    source_identity = sealed.repository_file(source)
    primary_targets = sealed.repository_file(REPO / "data/official/ni_hda_kraken.csv")
    primary_provenance = sealed.repository_file(REPO / "data/official/provenance.json")
    if (
        json.loads(Path(primary_provenance["path"]).read_text())["sha256"]
        != primary_targets["sha256"]
    ):
        raise ValueError("Published target table does not match source provenance")
    temperature_evidence = [
        sealed.check(name)
        for name in (
            "models/sources/correction_si.pdf",
            "models/sources/correction_si.txt",
            "models/inputs/corrected_ni_hda_tableS3_transcription_v2.json",
            "kinetics/results/ni_hda_target_temperature.csv",
        )
    ]
    with source.open(newline="") as stream:
        reader = csv.DictReader(stream)
        columns = list(reader.fieldnames or [])
        rows = list(reader)
    if len(rows) != 11 or {row["Source_ID"] for row in rows} != EXPECTED_IDS:
        raise ValueError(
            "This correction is scoped to the eleven audited Ni-hDA records"
        )
    if any(row["Temp_K"] != "298.15" for row in rows):
        raise ValueError(
            "Historical temperature changed; review provenance before migration"
        )
    output = args.output.resolve()
    if output.is_relative_to(REPO / "docs/scientific_accuracy_audit"):
        raise ValueError("Historical audit evidence is immutable")
    output.mkdir(parents=True, exist_ok=False)
    geometries = {}
    migrated = []
    for row in rows:
        new = dict(row)
        new["Temp_K"] = "353.15"
        new["Historical_Recorded_Temp_K"] = row["Temp_K"]
        new["Conformer_Population_Source"] = (
            "Supplied historical MMFF94 populations at 298.15 K; "
            "retained without reweighting"
        )
        new["Target_Source_Note"] = (
            "Published ddG_abs retained; source ee/DDG inconsistency remains unresolved"
            if row["Source_ID"] == "2064"
            else "Published ddG_abs retained"
        )
        for field in ("Ligand_XYZ_Path", "Conformer_XYZ_Paths"):
            paths = []
            for value in row[field].split(";"):
                path = (source.parent / value).resolve(strict=True)
                geometries[str(path)] = sealed.repository_file(path)
                paths.append(str(path))
            new[field] = ";".join(paths)
        migrated.append(new)
    destination = output / "reactions.csv"
    with destination.open("x", newline="") as stream:
        writer = csv.DictWriter(
            stream,
            fieldnames=[
                *columns,
                "Historical_Recorded_Temp_K",
                "Conformer_Population_Source",
                "Target_Source_Note",
            ],
        )
        writer.writeheader()
        writer.writerows(migrated)
    preserved = [
        name
        for name in columns
        if name not in {"Temp_K", "Ligand_XYZ_Path", "Conformer_XYZ_Paths"}
    ]
    if any(
        old[key] != new[key]
        for old, new in zip(rows, migrated, strict=True)
        for key in preserved
    ):
        raise AssertionError("An unrelated scientific input changed")
    claims_identity = sealed.check("claims.json")
    frozen_helper = output / "harness.py"
    frozen_helper.write_bytes(Path(__file__).read_bytes())
    if record(frozen_helper)["sha256"] != helper_identity["sha256"]:
        raise ValueError("Migration helper changed while being archived")
    sealed.recheck(
        [
            source_identity,
            primary_targets,
            primary_provenance,
            helper_identity,
            *geometries.values(),
        ]
    )
    data = {
        "kind": "Ni_hDA_response_temperature_metadata_correction",
        "source": source_identity,
        "output": record(destination),
        "helper": record(frozen_helper),
        "executed_helper": helper_identity,
        "audit_manifest": sealed.manifest,
        "audit_claims": claims_identity,
        "primary_target_table": primary_targets,
        "primary_target_provenance": primary_provenance,
        "temperature_primary_evidence": temperature_evidence,
        "sealed_inputs": list(sealed.used.values()),
        "pre_post_input_identities_verified": True,
        "authority": ["C21", "M67", "M68"],
        "temperature_k": 353.15,
        "conformer_populations_recomputed": False,
        "conformer_population_source_temperature_k": 298.15,
        "preserved_columns_verified": preserved,
        "geometry_inputs": list(geometries.values()),
        "not_a_new_experiment": True,
        "next_step": (
            "Recompute packed descriptors with the corrected native implementation; "
            "fit a new model with supplied_weight_mean "
            "and response temperature 353.15 K"
        ),
    }
    path = output / "manifest.json"
    path.write_text(json.dumps(data, indent=2) + "\n")
    (output / "manifest.sha256").write_text(record(path)["sha256"] + "\n")
    print(
        json.dumps(
            {"manifest": str(path), "records": len(rows), "geometries": len(geometries)}
        )
    )


if __name__ == "__main__":
    main()
