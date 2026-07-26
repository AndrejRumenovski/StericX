#!/usr/bin/env python3
"""Prepare xTB-LMO or CREST/xTB quantum inputs for StericX Study 003.

The default ``lmo`` mode preserves the existing conformer ensemble and replaces
the geometrically inferred virtual center with Kraken's xTB localized-orbital
selection. ``crest`` mode additionally replaces the ensemble and Boltzmann
weights with a pinned CREST 2.12/GFN2-xTB calculation.
"""

from __future__ import annotations

import argparse
import hashlib
import os
import re
import shutil
import sys
from collections.abc import Iterable
from datetime import UTC, datetime
from pathlib import Path

import pandas as pd

from stericx_quantum import (
    QuantumBackend,
    QuantumBackendError,
    QuantumConfig,
    atomic_write_json,
    sha256_file,
)


def parse_args(argv: Iterable[str] | None = None) -> argparse.Namespace:
    root = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(
        description="Generate quantum-derived conformers and coordination centers."
    )
    parser.add_argument(
        "--input-csv",
        type=Path,
        default=root / "data" / "reactions_raw.csv",
    )
    parser.add_argument("--xyz-root", type=Path, default=root / "data")
    parser.add_argument(
        "--output-csv",
        type=Path,
        default=root / "data" / "reactions_quantum.csv",
    )
    parser.add_argument(
        "--provenance",
        type=Path,
        default=root / "data" / "quantum" / "provenance.json",
    )
    parser.add_argument(
        "--cache-dir",
        type=Path,
        default=root / ".stericx" / "cache",
    )
    parser.add_argument(
        "--mode",
        choices=("lmo", "crest"),
        default="lmo",
        help="Derive LMO centers only or replace the complete ensemble.",
    )
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument(
        "--lmo-workers",
        type=int,
        default=4,
        help="Maximum concurrent xTB LMO processes (CPU-bounded automatically).",
    )
    parser.add_argument(
        "--stale-lock-seconds",
        type=float,
        default=86_400.0,
        help="Age threshold for reclaiming locks owned by another host.",
    )
    parser.add_argument("--charge", type=int, default=0)
    parser.add_argument("--uhf", type=int, default=0)
    parser.add_argument("--temperature", type=float, default=298.15)
    parser.add_argument("--energy-window", type=float, default=6.0)
    parser.add_argument("--center-distance", type=float, default=2.1)
    parser.add_argument(
        "--quick",
        action="store_true",
        help="Use CREST's reduced quick search; recorded as non-production.",
    )
    parser.add_argument("--max-records", type=int, default=None)
    args = parser.parse_args(argv)
    if args.max_records is not None and args.max_records <= 0:
        parser.error("--max-records must be positive")
    return args


def safe_stem(value: str) -> str:
    """Return a deterministic portable directory name."""
    stem = re.sub(r"[^A-Za-z0-9._-]+", "_", value).strip("._-") or "reaction"
    digest = hashlib.blake2s(value.encode(), digest_size=4).hexdigest()
    return f"{stem[:70]}_{digest}"


def semicolon_values(value: object) -> list[str]:
    """Split a non-empty semicolon-delimited CSV cell."""
    fields = [field.strip() for field in str(value).split(";") if field.strip()]
    if not fields:
        raise ValueError("conformer path list is empty")
    return fields


def atomic_write_csv(path: Path, frame: pd.DataFrame) -> None:
    """Atomically replace one tabular artifact."""
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    frame.to_csv(temporary, index=False, float_format="%.10g")
    os.replace(temporary, path)


def copy_immutable(source: Path, destination: Path) -> None:
    """Materialize one cached conformer without overwriting different data."""
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.is_file():
        if sha256_file(source) != sha256_file(destination):
            raise ValueError(f"refusing to overwrite different file: {destination}")
        return
    temporary = destination.with_name(destination.name + ".tmp")
    shutil.copy2(source, temporary)
    os.replace(temporary, destination)


def format_center(center: tuple[float, float, float]) -> str:
    """Serialize one center for the Rust CSV contract."""
    return ",".join(f"{coordinate:.10f}" for coordinate in center)


def main(argv: Iterable[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        source_bytes = args.input_csv.read_bytes()
        frame = pd.read_csv(args.input_csv)
        required = {
            "Reaction_ID",
            "Ligand_XYZ_Path",
            "Attach_Atom_Idx",
            "Conformer_XYZ_Paths",
            "Conformer_Boltzmann_Weights",
        }
        missing = sorted(required.difference(frame.columns))
        if missing:
            raise ValueError(f"input CSV lacks: {', '.join(missing)}")
        if args.max_records is not None:
            frame = frame.iloc[: args.max_records].copy()

        backend = QuantumBackend(
            QuantumConfig(
                cache_dir=args.cache_dir,
                threads=args.threads,
                charge=args.charge,
                uhf=args.uhf,
                temperature_k=args.temperature,
                energy_window_kcal_mol=args.energy_window,
                center_distance_angstrom=args.center_distance,
                quick=args.quick,
                lmo_workers=args.lmo_workers,
                stale_lock_seconds=args.stale_lock_seconds,
            )
        )
        output_rows: list[dict[str, object]] = []
        cache_hits = 0
        cache_misses = 0
        conformer_total = 0
        calculation_manifests: list[str] = []
        crest_cache_hits = 0
        crest_cache_misses = 0
        lmo_cache_hits = 0
        lmo_cache_misses = 0
        resumed_jobs = 0
        effective_lmo_workers: list[int] = []

        for position, (_, row) in enumerate(frame.iterrows(), start=1):
            reaction_id = str(row["Reaction_ID"])
            donor_idx = int(row["Attach_Atom_Idx"])
            output = row.to_dict()
            print(
                f"quantum_start={position}/{len(frame)},"
                f"reaction_id={reaction_id},mode={args.mode}",
                flush=True,
            )
            if args.mode == "lmo":
                relative_paths = semicolon_values(row["Conformer_XYZ_Paths"])
                centers: list[str] = []
                lmo_keys: list[str] = []
                batch = backend.lmo_centers(
                    [args.xyz_root / path for path in relative_paths],
                    donor_idx,
                )
                for result in batch.results:
                    centers.append(format_center(result.center_angstrom))
                    lmo_keys.append(result.cache_key)
                    calculation_manifests.append(result.manifest_path)
                cache_hits += batch.cache_hits
                cache_misses += batch.cache_misses
                lmo_cache_hits += batch.cache_hits
                lmo_cache_misses += batch.cache_misses
                effective_lmo_workers.append(batch.effective_workers)
                output["Conformer_Coordination_Centers_Angstrom"] = ";".join(centers)
                output["Coordination_Center_Method"] = "xtb_lmo_kraken"
                output["Quantum_LMO_Cache_Keys"] = ";".join(lmo_keys)
                output["Quantum_Ensemble_Cache_Key"] = ""
                output["Quantum_CREST_Cache_Key"] = ""
                conformer_total += len(relative_paths)
            else:
                input_xyz = args.xyz_root / str(row["Ligand_XYZ_Path"])
                ensemble = backend.conformer_ensemble(input_xyz, donor_idx)
                destination_root = (
                    args.output_csv.parent
                    / "quantum"
                    / "conformers"
                    / safe_stem(reaction_id)
                    / ensemble.cache_key[:16]
                )
                relative_paths = []
                centers = []
                energies = []
                weights = []
                lmo_keys = []
                for conformer_index, conformer in enumerate(ensemble.conformers):
                    destination = destination_root / f"conf_{conformer_index:03d}.xyz"
                    copy_immutable(Path(conformer.xyz_path), destination)
                    relative_paths.append(
                        destination.relative_to(args.output_csv.parent).as_posix()
                    )
                    centers.append(
                        format_center(conformer.coordination_center_angstrom)
                    )
                    energies.append(f"{conformer.relative_energy_kcal_mol:.10f}")
                    weights.append(f"{conformer.boltzmann_weight:.12f}")
                    lmo_keys.append(conformer.lmo_cache_key)
                output["Ligand_XYZ_Path"] = relative_paths[0]
                output["Conformer_XYZ_Paths"] = ";".join(relative_paths)
                output["Conformer_Relative_Energies_kcal_mol"] = ";".join(energies)
                output["Conformer_Boltzmann_Weights"] = ";".join(weights)
                output["Conformer_Count"] = len(relative_paths)
                output["Ensemble_Energy_Span_kcal_mol"] = max(
                    conformer.relative_energy_kcal_mol
                    for conformer in ensemble.conformers
                )
                output["Conformer_Coordination_Centers_Angstrom"] = ";".join(centers)
                output["Coordination_Center_Method"] = "crest212_xtb640_lmo_kraken"
                output["Quantum_LMO_Cache_Keys"] = ";".join(lmo_keys)
                output["Quantum_Ensemble_Cache_Key"] = ensemble.cache_key
                output["Quantum_CREST_Cache_Key"] = ensemble.crest_cache_key
                calculation_manifests.append(ensemble.manifest_path)
                cache_hits += int(ensemble.cache_hit)
                cache_misses += int(not ensemble.cache_hit)
                crest_cache_hits += int(ensemble.crest_cache_hit)
                crest_cache_misses += int(not ensemble.crest_cache_hit)
                lmo_cache_hits += ensemble.lmo_cache_hits
                lmo_cache_misses += ensemble.lmo_cache_misses
                resumed_jobs += int(ensemble.job_resumed)
                effective_lmo_workers.append(ensemble.effective_lmo_workers)
                conformer_total += len(relative_paths)
            output["Quantum_xTB_Version"] = backend.xtb.version
            output["Quantum_CREST_Version"] = backend.crest.version
            output_rows.append(output)
            print(
                f"[{position:03d}/{len(frame):03d}] {reaction_id}: "
                f"mode={args.mode}, cumulative_conformers={conformer_total}"
            )

        output_frame = pd.DataFrame.from_records(output_rows)
        atomic_write_csv(args.output_csv, output_frame)
        provenance = {
            "schema_version": 2,
            "generated_at_utc": datetime.now(UTC).isoformat(),
            "mode": args.mode,
            "production_profile": args.mode == "crest" and not args.quick,
            "input_csv": str(args.input_csv),
            "input_csv_sha256": hashlib.sha256(source_bytes).hexdigest(),
            "output_csv": str(args.output_csv),
            "output_csv_sha256": sha256_file(args.output_csv),
            "records": len(output_frame),
            "conformers": conformer_total,
            "cache_hits": cache_hits,
            "cache_misses": cache_misses,
            "cache_stages": {
                "crest_hits": crest_cache_hits,
                "crest_misses": crest_cache_misses,
                "lmo_hits": lmo_cache_hits,
                "lmo_misses": lmo_cache_misses,
                "resumed_jobs": resumed_jobs,
            },
            "xtb": {
                "version": backend.xtb.version,
                "path": backend.xtb.path,
                "sha256": backend.xtb.sha256,
            },
            "crest": {
                "version": backend.crest.version,
                "path": backend.crest.path,
                "sha256": backend.crest.sha256,
            },
            "settings": {
                "threads": args.threads,
                "charge": args.charge,
                "uhf": args.uhf,
                "temperature_k": args.temperature,
                "energy_window_kcal_mol": args.energy_window,
                "center_distance_angstrom": args.center_distance,
                "solvent": "toluene",
                "quick": args.quick,
                "lmo_workers_requested": args.lmo_workers,
                "lmo_workers_effective_max": max(effective_lmo_workers, default=0),
                "stale_lock_seconds": args.stale_lock_seconds,
            },
            "calculation_manifests": sorted(set(calculation_manifests)),
        }
        atomic_write_json(args.provenance, provenance)
        print("\nQuantum data preparation complete")
        print(f"  mode={args.mode}")
        print(f"  records={len(output_frame)}")
        print(f"  conformers={conformer_total}")
        print(f"  cache_hits={cache_hits}")
        print(f"  cache_misses={cache_misses}")
        print(f"  crest_cache_hits={crest_cache_hits}")
        print(f"  crest_cache_misses={crest_cache_misses}")
        print(f"  lmo_cache_hits={lmo_cache_hits}")
        print(f"  lmo_cache_misses={lmo_cache_misses}")
        print(f"  resumed_jobs={resumed_jobs}")
        print(f"  lmo_workers_effective_max={max(effective_lmo_workers, default=0)}")
        print(f"  output={args.output_csv}")
        print(f"  provenance={args.provenance}")
        return 0
    except (
        FileNotFoundError,
        OSError,
        QuantumBackendError,
        ValueError,
    ) as exc:
        print(f"Quantum data preparation failed: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
