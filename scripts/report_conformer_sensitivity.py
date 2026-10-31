#!/usr/bin/env python3
"""Report the retained repository MMFF ensembles through the Rust descriptor API.

This is orchestration/provenance only. No chemistry or descriptor reductions
are implemented in Python. Historical Kraken energies are never inferred.
"""

import argparse
import csv
import hashlib
import json
import subprocess
import sys
from datetime import UTC, datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--engine", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--source", type=Path, default=ROOT / "data/reactions_raw.csv")
    args = parser.parse_args()
    engine, source, output = (
        args.engine.resolve(),
        args.source.resolve(),
        args.output.resolve(),
    )
    output.mkdir(parents=True, exist_ok=False)
    manifest = {
        "created_utc": datetime.now(UTC).isoformat(),
        "command": sys.argv,
        "script_sha256": sha(__file__),
        "engine_sha256": sha(engine),
        "source": str(source.relative_to(ROOT)),
        "source_sha256": sha(source),
        "interpretation": "Retained-state MMFF-energy approximation, not experimental or historical Kraken populations.",  # noqa: E501
        "ligands": [],
        "failures": [],
    }
    rows = list(csv.DictReader(source.open()))
    tables = []
    for row in rows:
        identity = row["Reaction_ID"]
        directory = output / identity
        directory.mkdir()
        record = {"ligand_id": identity, "status": "pending", "inputs": []}
        manifest["ligands"].append(record)
        try:
            paths = [
                source.parent / name for name in row["Conformer_XYZ_Paths"].split(";")
            ]
            energies = [
                float(value)
                for value in row["Conformer_Relative_Energies_kcal_mol"].split(";")
            ]
            if len(paths) != len(energies) or len(paths) != int(row["Conformer_Count"]):
                raise ValueError("inconsistent recorded conformer membership/energies")
            texts = []
            for path in paths:
                text = path.read_text()
                if "MMFF94/ETKDGv3 ensemble" not in text.splitlines()[1]:
                    raise ValueError(
                        f"source does not declare the expected MMFF model: {path}"
                    )
                texts.append(text if text.endswith("\n") else text + "\n")
                record["inputs"].append(
                    {"path": str(path.relative_to(ROOT)), "sha256": sha(path)}
                )
            (directory / "ensemble.xyz").write_text("".join(texts))
            metadata = [
                {
                    "value": value,
                    "unit": "kcal_mol",
                    "kind": "force_field",
                    "model": "MMFF94 retained ETKDGv3 geometries; repository prepare_data workflow",  # noqa: E501
                    "provenance": f"{source.relative_to(ROOT)}#{identity}/Conformer_Relative_Energies_kcal_mol[{i}]; geometry={path.relative_to(ROOT)}",  # noqa: E501
                    "temperature_k": None,
                    "degeneracy": 1,
                }
                for i, (path, value) in enumerate(zip(paths, energies, strict=True))
            ]
            (directory / "energies.json").write_text(
                json.dumps(metadata, indent=2) + "\n"
            )
            options = {
                "axis": "coordination",
                "donor_index": int(row["Attach_Atom_Idx"]),
            }
            (directory / "options.json").write_text(
                json.dumps(options, indent=2) + "\n"
            )
            command = [
                str(engine),
                str(directory / "ensemble.xyz"),
                str(directory / "report"),
                "--ligand-id",
                identity,
                "--options",
                str(directory / "options.json"),
                "--energies",
                str(directory / "energies.json"),
                "--temperature",
                row["Temp_K"],
                "--population-assumptions",
                "Retained records at unit multiplicity; no merging or symmetry correction. MMFF electronic-like energy approximation; completeness, distinct-state membership and thermochemical entropy not established.",  # noqa: E501
            ]
            process = subprocess.run(
                command, text=True, capture_output=True, check=False
            )
            record.update(command=command, returncode=process.returncode)
            (directory / "stdout.txt").write_text(process.stdout)
            (directory / "stderr.txt").write_text(process.stderr)
            if process.returncode:
                raise RuntimeError(
                    f"engine exited {process.returncode}; see stderr.txt"
                )
            table = list(csv.DictReader((directory / "report.csv").open()))
            tables.extend(table)
            record.update(
                status="success",
                conformer_count=len(paths),
                files={
                    path.name: sha(path)
                    for path in sorted(directory.iterdir())
                    if path.is_file()
                },
            )
        except (OSError, ValueError, RuntimeError) as error:
            record.update(status="error", error=str(error))
            manifest["failures"].append({"ligand_id": identity, "error": str(error)})
    if tables:
        with (output / "sensitivity.csv").open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(tables[0]))
            writer.writeheader()
            writer.writerows(tables)
        manifest["sensitivity_csv_sha256"] = sha(output / "sensitivity.csv")
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(
        json.dumps(
            {
                "ligands": len(rows),
                "table_rows": len(tables),
                "failures": manifest["failures"],
            },
            indent=2,
        )
    )
    return bool(manifest["failures"])


if __name__ == "__main__":
    raise SystemExit(main())
