#!/usr/bin/env python3
"""Prepare provenance-rich ligand ensembles for the StericX engine.

The default source is the Sigman Group's public Ni-catalyzed homo-Diels-Alder
Kraken table.  The complete 1,566-ligand library is preserved under
``data/official/`` while rows with measured ee/ΔΔG‡ are converted into
conformer ensembles and the fixed-width reaction matrix.  If the public source
cannot be downloaded or normalized, a deterministic 100-reaction synthetic
benchmark embedded in this file is used instead.

Dependencies:
    rdkit, pandas, requests, numpy

Run:
    python scripts/prepare_data.py

Outputs:
    data/xyz/<Reaction_ID>.xyz
    data/conformers/<Reaction_ID>/conf_<N>.xyz
    data/reactions_raw.csv
    data/official/ni_hda_kraken.csv
    data/official/provenance.json
    data/etl_failures.csv             # only when individual rows fail
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import os
import re
import sys
import time
from collections.abc import Iterable
from dataclasses import dataclass
from datetime import UTC, datetime
from pathlib import Path
from typing import Final

try:
    import numpy as np
    import pandas as pd
    import requests
    from rdkit import Chem, RDLogger
    from rdkit.Chem import AllChem
    from rdkit.Chem.Scaffolds import MurckoScaffold
except ImportError as exc:  # pragma: no cover - depends on the host environment.
    missing = getattr(exc, "name", "required package")
    raise SystemExit(
        f"Missing dependency `{missing}`. Install rdkit, pandas, requests, and "
        "numpy before running scripts/prepare_data.py."
    ) from exc


PUBLIC_DATASET_URLS: Final[tuple[str, ...]] = (
    (
        "https://raw.githubusercontent.com/SigmanGroup/"
        "Ni-Catalyzed-hDA/main/data/kraken.csv"
    ),
)
DEFAULT_SUBSTRATE_SMILES: Final[str] = "C=CC(=O)C"
DEFAULT_TEMPERATURE_K: Final[float] = 298.15
DEFAULT_IR_FREQUENCY_CM1: Final[float] = 1650.0
GAS_CONSTANT_KCAL_MOL_K: Final[float] = 0.00198720425864083
PUBLISHED_TRAIN_IDS: Final[frozenset[str]] = frozenset(
    {"401", "498", "724", "785", "1057", "1058", "2062", "2063", "2064", "2067"}
)
HISTORICAL_BLIND_IDS: Final[frozenset[str]] = frozenset({"723"})
OUTPUT_COLUMNS: Final[list[str]] = [
    "Reaction_ID",
    "Ligand_XYZ_Path",
    "Attach_Atom_Idx",
    "Primary_Bond_Vector_Idx",
    "NBO_Charge",
    "IR_Frequency",
    "Temp_K",
    "Exp_ddG_kcal_mol",
    "Conformer_XYZ_Paths",
    "Conformer_Relative_Energies_kcal_mol",
    "Conformer_Boltzmann_Weights",
    "Conformer_Count",
    "Ensemble_Energy_Span_kcal_mol",
    "Ligand_SMILES",
    "Ligand_Group",
    "Dataset_Split",
    "Source_ID",
    "Source_URL",
]

# These structures span P- and N-donor ligand families and are intentionally
# repeated under different reaction conditions/substrates by the fallback
# builder.  The resulting FALLBACK_CSV contains exactly 100 data rows.
_FALLBACK_LIGANDS: Final[tuple[str, ...]] = (
    "P(c1ccccc1)(c1ccccc1)c1ccccc1",
    "P(C1CCCCC1)(C1CCCCC1)C1CCCCC1",
    "CCP(CC)CC",
    "C1CCP(CC1)c1ccccc1",
    "COc1cccc(P(c2ccccc2)c2ccccc2)c1",
    "CC(C)c1cc(C(C)C)c(P(c2ccccc2)c2ccccc2)c(C(C)C)c1",
    "CC(C)(C)c1cc(P(c2ccccc2)c2ccccc2)cc(C(C)(C)C)c1",
    "CP(C)c1ccccc1",
    "n1ccccc1",
    "Cc1cccc(C)n1",
    "CC(C)c1cccc(C(C)C)n1",
    "CN(C)c1ccccc1",
    "N1CCCCC1",
    "CN1CCCCC1",
    "CCN(CC)CC",
    "c1ccc2ncccc2c1",
    "N[C@@H](C)c1ccccc1",
    "O=C(O)[C@@H](N)Cc1ccccc1",
    "C1(C2=N[C@@H](Cc3ccccc3)CO2)=N[C@@H](Cc2ccccc2)CO1",
    "C1(C2=N[C@@H](C(C)C)CO2)=N[C@@H](C(C)C)CO1",
)
_FALLBACK_SUBSTRATES: Final[tuple[str, ...]] = (
    "C=CC(=O)C",
    "C=CC(=O)OC",
    "C=CC(=O)N",
    "O=C(C)c1ccccc1",
    "O=C1CCCCC1",
    "CC(=O)c1ccc(F)cc1",
    "N#Cc1ccccc1",
    "O=C(CCl)c1ccccc1",
    "CCOC(=O)C=C",
    "O=C(c1ccccc1)c1ccccc1",
)
_FALLBACK_NBO_BASE: Final[tuple[float, ...]] = (
    0.82,
    0.77,
    0.69,
    0.74,
    0.86,
    0.91,
    0.88,
    0.71,
    -0.48,
    -0.46,
    -0.44,
    -0.51,
    -0.57,
    -0.55,
    -0.59,
    -0.43,
    -0.62,
    -0.64,
    -0.49,
    -0.52,
)


def _build_fallback_csv() -> str:
    """Build the embedded deterministic 100-row CSV benchmark."""
    buffer = io.StringIO()
    writer = csv.writer(buffer, lineterminator="\n")
    writer.writerow(
        [
            "Reaction_ID",
            "Ligand_SMILES",
            "Substrate_SMILES",
            "NBO_Charge_N",
            "IR_Stretching_Freq",
            "Temp_K",
            "Experimental_ddG",
        ]
    )
    temperatures = (243.15, 258.15, 273.15, 288.15, 298.15, 313.15, 323.15)
    for index in range(100):
        ligand_index = index % len(_FALLBACK_LIGANDS)
        replicate = index // len(_FALLBACK_LIGANDS)
        nbo = _FALLBACK_NBO_BASE[ligand_index] + 0.0125 * (replicate - 2)
        ir_frequency = 1582.0 + float((index * 37 + ligand_index * 11) % 157)
        temperature = temperatures[(index * 3 + ligand_index) % len(temperatures)]
        donor_term = -0.36 * nbo
        spectral_term = 0.0095 * (ir_frequency - 1650.0)
        periodic_term = 0.78 * np.sin((index + 1) * 0.71)
        ddg = float(
            np.clip(0.72 + donor_term + spectral_term + periodic_term, -2.4, 2.8)
        )
        writer.writerow(
            [
                f"SYN-{index + 1:03d}",
                _FALLBACK_LIGANDS[ligand_index],
                _FALLBACK_SUBSTRATES[(index * 3) % len(_FALLBACK_SUBSTRATES)],
                f"{nbo:.5f}",
                f"{ir_frequency:.2f}",
                f"{temperature:.2f}",
                f"{ddg:.6f}",
            ]
        )
    return buffer.getvalue()


FALLBACK_CSV: Final[str] = _build_fallback_csv()


@dataclass(frozen=True)
class GeometryResult:
    """Indices, conformer paths, and energies for one ligand ensemble."""

    xyz_relative_path: str
    conformer_relative_paths: tuple[str, ...]
    relative_energies_kcal_mol: tuple[float, ...]
    boltzmann_weights: tuple[float, ...]
    attach_atom_idx: int
    primary_bond_vector_idx: int
    mmff_statuses: tuple[int, ...]


@dataclass(frozen=True)
class ConformerEnsemble:
    """One hydrogenated molecule and its retained MMFF94 conformers."""

    molecule: Chem.Mol
    conformer_ids: tuple[int, ...]
    relative_energies_kcal_mol: tuple[float, ...]
    boltzmann_weights: tuple[float, ...]
    mmff_statuses: tuple[int, ...]


@dataclass(frozen=True)
class LoadedDataset:
    """Normalized modeling rows plus the complete downloaded source table."""

    reactions: pd.DataFrame
    complete_catalog: pd.DataFrame | None
    complete_catalog_bytes: bytes | None
    source_name: str


@dataclass(frozen=True)
class ProcessingFailure:
    """A source row that could not be converted into a usable geometry."""

    reaction_id: str
    ligand_smiles: str
    reason: str


class DatasetError(RuntimeError):
    """Raised when a downloaded table cannot be used as a reaction dataset."""


def parse_args(argv: Iterable[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Download an asymmetric-catalysis benchmark, preserve its provenance, "
            "generate MMFF94 ligand conformer ensembles, and export StericX inputs."
        )
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("data"),
        help="Output directory (default: data).",
    )
    parser.add_argument(
        "--dataset-url",
        action="append",
        dest="dataset_urls",
        help="CSV URL to try before the built-in Sigman source; may be repeated.",
    )
    parser.add_argument(
        "--offline",
        action="store_true",
        help="Skip network access and use the embedded 100-row benchmark.",
    )
    parser.add_argument(
        "--timeout",
        type=float,
        default=30.0,
        help="HTTP timeout in seconds (default: 30).",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=20260725,
        help="Base ETKDG random seed for reproducible conformers.",
    )
    parser.add_argument(
        "--conformers",
        type=int,
        default=20,
        help="Maximum ETKDGv3 conformers generated per ligand (default: 20).",
    )
    parser.add_argument(
        "--prune-rms",
        type=float,
        default=0.5,
        help="ETKDG heavy-atom RMS pruning threshold in Å (default: 0.5).",
    )
    parser.add_argument(
        "--energy-window",
        type=float,
        default=5.0,
        help="Retain conformers within this many kcal/mol of the minimum (default: 5).",
    )
    parser.add_argument(
        "--max-records",
        type=int,
        default=None,
        help="Optionally process only the first N normalized records.",
    )
    parser.add_argument(
        "--strict",
        action="store_true",
        help="Abort on the first invalid SMILES, embedding, or MMFF failure.",
    )
    args = parser.parse_args(argv)
    if args.timeout <= 0.0:
        parser.error("--timeout must be positive")
    if args.max_records is not None and args.max_records <= 0:
        parser.error("--max-records must be positive")
    if args.conformers <= 0:
        parser.error("--conformers must be positive")
    if args.prune_rms < 0.0:
        parser.error("--prune-rms must be non-negative")
    if args.energy_window <= 0.0:
        parser.error("--energy-window must be positive")
    return args


def download_csv(
    url: str, timeout: float, attempts: int = 3
) -> tuple[pd.DataFrame, bytes]:
    """Download one CSV with bounded retries and basic response validation."""
    headers = {
        "Accept": "text/csv,text/plain;q=0.9,*/*;q=0.1",
        "User-Agent": "steric-x-data-preparation/1.0",
    }
    last_error: Exception | None = None
    for attempt in range(1, attempts + 1):
        try:
            response = requests.get(
                url,
                headers=headers,
                timeout=(min(10.0, timeout), timeout),
            )
            response.raise_for_status()
            if not response.content.strip():
                raise DatasetError("server returned an empty response")
            frame = pd.read_csv(io.BytesIO(response.content))
            if frame.empty:
                raise DatasetError("downloaded CSV has no data rows")
            return frame, response.content
        except (requests.RequestException, pd.errors.ParserError, DatasetError) as exc:
            last_error = exc
            if attempt < attempts:
                time.sleep(float(attempt))
    raise DatasetError(f"download failed after {attempts} attempts: {last_error}")


def ee_to_ddg(ee_percent: pd.Series, temperature_k: float) -> pd.Series:
    """Convert absolute percent ee to a positive ΔΔG‡ magnitude."""
    ee = pd.to_numeric(ee_percent, errors="coerce").abs().clip(lower=0.0, upper=99.999)
    major = (100.0 + ee) / 2.0
    minor = (100.0 - ee) / 2.0
    return GAS_CONSTANT_KCAL_MOL_K * temperature_k * np.log(major / minor)


def normalize_public_sigman(frame: pd.DataFrame) -> pd.DataFrame:
    """Normalize the official Ni-hDA/Kraken table to the ETL source schema."""
    columns = {str(column).strip(): column for column in frame.columns}
    smiles_column = columns.get("smiles")
    if smiles_column is None:
        raise DatasetError("public table has no `smiles` column")

    ddg = pd.to_numeric(frame.get("ddG_abs"), errors="coerce")
    if ddg.isna().all() and "ee" not in frame.columns:
        raise DatasetError("public table has neither experimental `ddG_abs` nor `ee`")
    if "ee" in frame.columns:
        ddg = ddg.fillna(ee_to_ddg(frame["ee"], DEFAULT_TEMPERATURE_K))

    source_ids = (
        frame.iloc[:, 0].astype(str).str.replace(r"\.0$", "", regex=True).str.strip()
    )
    source_ids = source_ids.where(source_ids.ne(""), frame.index.astype(str))
    donor_charge = pd.to_numeric(frame.get("nbo_P_boltz"), errors="coerce")
    dataset_split = source_ids.map(
        lambda source_id: (
            "train"
            if source_id in PUBLISHED_TRAIN_IDS
            else "blind"
            if source_id in HISTORICAL_BLIND_IDS
            else "unlabeled"
        )
    )
    normalized = pd.DataFrame(
        {
            "Reaction_ID": "SIG-NIHDA-" + source_ids,
            "Source_ID": source_ids,
            "Ligand_SMILES": frame[smiles_column].astype(str).str.strip(),
            "Substrate_SMILES": DEFAULT_SUBSTRATE_SMILES,
            # The source contains P-donor NBO charge. It maps to the engine's
            # generic donor-charge slot despite the historical `_N` CSV name.
            "NBO_Charge_N": donor_charge,
            # This benchmark does not report IR values or per-row temperature.
            # Explicit benchmark defaults keep the fixed-width matrix complete.
            "IR_Stretching_Freq": DEFAULT_IR_FREQUENCY_CM1,
            "Temp_K": DEFAULT_TEMPERATURE_K,
            "Experimental_ddG": ddg,
            "Dataset_Split": dataset_split,
        }
    )
    labeled_count = int(normalized["Experimental_ddG"].notna().sum())
    if labeled_count == 0:
        raise DatasetError("public table contains no measured ee/ΔΔG‡ rows")
    print(
        "  Source note: mapped `nbo_P_boltz` to the generic donor-charge field; "
        f"the {DEFAULT_IR_FREQUENCY_CM1:.1f} cm^-1 IR value is an explicit "
        "missing-data placeholder and is excluded from scientific feature "
        f"selection; temperature is {DEFAULT_TEMPERATURE_K:.2f} K."
    )
    return normalized


def load_source_dataset(
    offline: bool,
    custom_urls: list[str] | None,
    timeout: float,
) -> LoadedDataset:
    """Load the first usable public table, otherwise return embedded fallback."""
    if not offline:
        urls = tuple(custom_urls or ()) + PUBLIC_DATASET_URLS
        for url in urls:
            print(f"Fetching public benchmark: {url}")
            try:
                downloaded, source_bytes = download_csv(url, timeout)
                normalized = normalize_public_sigman(downloaded)
                labeled = normalized.loc[
                    normalized["Experimental_ddG"].notna()
                ].reset_index(drop=True)
                print(
                    f"  Downloaded {len(downloaded):,} rows; "
                    f"{len(labeled):,} have experimental selectivity."
                )
                return LoadedDataset(labeled, downloaded, source_bytes, url)
            except (DatasetError, KeyError, TypeError, ValueError) as exc:
                print(f"  Warning: source unavailable or incompatible: {exc}")
    else:
        print("Offline mode selected; skipping public dataset download.")

    fallback = pd.read_csv(io.StringIO(FALLBACK_CSV))
    if len(fallback) != 100:
        raise DatasetError(
            "embedded benchmark invariant failed: "
            f"expected 100 rows, got {len(fallback)}"
        )
    print("Using embedded 100-reaction synthetic benchmark.")
    fallback["Source_ID"] = fallback["Reaction_ID"]
    fallback["Dataset_Split"] = "train"
    return LoadedDataset(
        fallback,
        None,
        None,
        "embedded://steric_x-synthetic-100",
    )


def validate_source_frame(frame: pd.DataFrame) -> pd.DataFrame:
    """Validate and coerce the normalized source schema."""
    required = {
        "Reaction_ID",
        "Source_ID",
        "Ligand_SMILES",
        "Substrate_SMILES",
        "NBO_Charge_N",
        "IR_Stretching_Freq",
        "Temp_K",
        "Experimental_ddG",
        "Dataset_Split",
    }
    missing = sorted(required.difference(frame.columns))
    if missing:
        raise DatasetError(
            f"normalized dataset is missing columns: {', '.join(missing)}"
        )

    clean = frame.copy()
    clean["Reaction_ID"] = clean["Reaction_ID"].astype(str).str.strip()
    clean["Source_ID"] = clean["Source_ID"].astype(str).str.strip()
    clean["Ligand_SMILES"] = clean["Ligand_SMILES"].astype(str).str.strip()
    clean["Substrate_SMILES"] = clean["Substrate_SMILES"].astype(str).str.strip()
    clean["Dataset_Split"] = clean["Dataset_Split"].astype(str).str.strip()
    numeric_columns = [
        "NBO_Charge_N",
        "IR_Stretching_Freq",
        "Temp_K",
        "Experimental_ddG",
    ]
    for column in numeric_columns:
        clean[column] = pd.to_numeric(clean[column], errors="coerce")

    finite_mask = np.isfinite(clean[numeric_columns].to_numpy(dtype=float)).all(axis=1)
    content_mask = (
        clean["Reaction_ID"].ne("")
        & clean["Source_ID"].ne("")
        & clean["Ligand_SMILES"].ne("")
        & clean["Substrate_SMILES"].ne("")
        & clean["Dataset_Split"].isin({"train", "external", "blind"})
        & clean["Temp_K"].gt(0.0)
    )
    valid_mask = finite_mask & content_mask.to_numpy()
    dropped = int((~valid_mask).sum())
    if dropped:
        print(
            f"Warning: dropping {dropped} rows with missing/non-finite required values."
        )
    clean = clean.loc[valid_mask].reset_index(drop=True)
    if clean.empty:
        raise DatasetError("no valid rows remain after source validation")

    duplicate_mask = clean["Reaction_ID"].duplicated(keep=False)
    if duplicate_mask.any():
        duplicate_ids = clean.loc[duplicate_mask, "Reaction_ID"].unique()
        raise DatasetError(
            "reaction identifiers must be unique; duplicates include "
            + ", ".join(map(str, duplicate_ids[:5]))
        )
    return clean


def safe_filename(reaction_id: str) -> str:
    """Create a portable filename while retaining a collision-resistant suffix."""
    stem = re.sub(r"[^A-Za-z0-9._-]+", "_", reaction_id).strip("._-") or "reaction"
    digest = hashlib.blake2s(reaction_id.encode("utf-8"), digest_size=4).hexdigest()
    return f"{stem[:80]}_{digest}.xyz"


def branch_size(start_idx: int, blocked_idx: int, molecule: Chem.Mol) -> int:
    """Count atoms reachable from one neighbor without crossing the donor atom."""
    seen = {blocked_idx}
    stack = [start_idx]
    size = 0
    while stack:
        atom_idx = stack.pop()
        if atom_idx in seen:
            continue
        seen.add(atom_idx)
        size += 1
        atom = molecule.GetAtomWithIdx(atom_idx)
        stack.extend(neighbor.GetIdx() for neighbor in atom.GetNeighbors())
    return size


def identify_donor_and_axis(molecule: Chem.Mol) -> tuple[int, int]:
    """Select an N/P donor and a bonded heavy atom defining the primary axis."""
    candidates: list[tuple[tuple[int, int, int, int], Chem.Atom]] = []
    for atom in molecule.GetAtoms():
        atomic_number = atom.GetAtomicNum()
        if atomic_number not in (7, 15):
            continue
        heavy_neighbors = [
            neighbor for neighbor in atom.GetNeighbors() if neighbor.GetAtomicNum() > 1
        ]
        if not heavy_neighbors:
            continue
        priority = 200 if atomic_number == 15 else 100
        if atom.GetFormalCharge() > 0:
            priority -= 80
        if atom.GetIsAromatic():
            priority += 10
        if atom.GetTotalNumHs() > 0:
            priority -= 5
        score = (
            priority,
            len(heavy_neighbors),
            atom.GetDegree(),
            -atom.GetIdx(),
        )
        candidates.append((score, atom))
    if not candidates:
        raise ValueError("no bonded nitrogen or phosphorus donor atom was found")

    donor = max(candidates, key=lambda item: item[0])[1]
    heavy_neighbors = [
        neighbor for neighbor in donor.GetNeighbors() if neighbor.GetAtomicNum() > 1
    ]
    axis_neighbor = max(
        heavy_neighbors,
        key=lambda atom: (
            branch_size(atom.GetIdx(), donor.GetIdx(), molecule),
            atom.GetAtomicNum(),
            -atom.GetIdx(),
        ),
    )
    return donor.GetIdx(), axis_neighbor.GetIdx()


def embed_and_optimize(
    smiles: str,
    seed: int,
    conformer_count: int,
    prune_rms: float,
    energy_window: float,
    temperature_k: float,
) -> ConformerEnsemble:
    """Generate, MMFF94-optimize, and Boltzmann-weight a ligand ensemble."""
    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None:
        raise ValueError("RDKit could not parse ligand SMILES")
    molecule = Chem.AddHs(molecule)

    params = AllChem.ETKDGv3()
    params.randomSeed = int(seed % 2_147_483_647)
    params.pruneRmsThresh = float(prune_rms)
    params.numThreads = 1
    params.useSmallRingTorsions = True
    params.useMacrocycleTorsions = True
    conformer_ids = tuple(
        int(conformer_id)
        for conformer_id in AllChem.EmbedMultipleConfs(
            molecule,
            numConfs=int(conformer_count),
            params=params,
        )
    )
    if not conformer_ids:
        # A deterministic random-coordinate retry rescues difficult macrocycles.
        retry = AllChem.ETKDGv3()
        retry.randomSeed = int((seed + 7_919) % 2_147_483_647)
        retry.useRandomCoords = True
        retry.pruneRmsThresh = float(prune_rms)
        retry.numThreads = 1
        retry.useSmallRingTorsions = True
        retry.useMacrocycleTorsions = True
        conformer_ids = tuple(
            int(conformer_id)
            for conformer_id in AllChem.EmbedMultipleConfs(
                molecule,
                numConfs=int(conformer_count),
                params=retry,
            )
        )
    if not conformer_ids:
        raise ValueError("ETKDGv3 conformer embedding failed")
    if not AllChem.MMFFHasAllMoleculeParams(molecule):
        raise ValueError("MMFF94 parameters are unavailable for this ligand")

    optimization = AllChem.MMFFOptimizeMoleculeConfs(
        molecule,
        numThreads=0,
        maxIters=1_000,
        mmffVariant="MMFF94",
        nonBondedThresh=100.0,
    )
    if len(optimization) != len(conformer_ids):
        raise ValueError("MMFF94 returned an inconsistent conformer result count")

    optimized: list[tuple[int, int, float]] = []
    for conformer_id, (status, energy) in zip(conformer_ids, optimization, strict=True):
        status = int(status)
        energy = float(energy)
        if status < 0 or not np.isfinite(energy):
            continue
        optimized.append((conformer_id, status, energy))
    if not optimized:
        raise ValueError("MMFF94 failed for every embedded conformer")

    optimized.sort(key=lambda item: (item[2], item[0]))
    minimum_energy = optimized[0][2]
    retained = [item for item in optimized if item[2] - minimum_energy <= energy_window]
    relative_energies = np.asarray(
        [item[2] - minimum_energy for item in retained],
        dtype=float,
    )
    thermal_energy = GAS_CONSTANT_KCAL_MOL_K * temperature_k
    raw_weights = np.exp(-relative_energies / thermal_energy)
    weights = raw_weights / raw_weights.sum()

    return ConformerEnsemble(
        molecule=molecule,
        conformer_ids=tuple(item[0] for item in retained),
        relative_energies_kcal_mol=tuple(float(value) for value in relative_energies),
        boltzmann_weights=tuple(float(value) for value in weights),
        mmff_statuses=tuple(item[1] for item in retained),
    )


def write_xyz(
    molecule: Chem.Mol,
    destination: Path,
    reaction_id: str,
    attach_idx: int,
    axis_idx: int,
    conformer_id: int,
    relative_energy: float,
    boltzmann_weight: float,
) -> None:
    """Write one conformer as a standard XYZ file via atomic replacement."""
    conformer = molecule.GetConformer(conformer_id)
    lines = [
        str(molecule.GetNumAtoms()),
        (
            f"{reaction_id} | MMFF94/ETKDGv3 ensemble | "
            f"attach_idx={attach_idx} axis_idx={axis_idx} "
            f"relative_energy_kcal_mol={relative_energy:.8f} "
            f"boltzmann_weight={boltzmann_weight:.10f}"
        ),
    ]
    for atom in molecule.GetAtoms():
        position = conformer.GetAtomPosition(atom.GetIdx())
        coordinates = np.asarray([position.x, position.y, position.z], dtype=float)
        if not np.isfinite(coordinates).all():
            raise ValueError(f"non-finite coordinate at atom {atom.GetIdx()}")
        lines.append(
            f"{atom.GetSymbol():<3s} "
            f"{position.x: .10f} {position.y: .10f} {position.z: .10f}"
        )
    temporary = destination.with_suffix(destination.suffix + ".tmp")
    temporary.write_text("\n".join(lines) + "\n", encoding="utf-8")
    os.replace(temporary, destination)


def process_reactions(
    source: pd.DataFrame,
    output_dir: Path,
    base_seed: int,
    strict: bool,
    conformer_count: int,
    prune_rms: float,
    energy_window: float,
    source_name: str,
) -> tuple[pd.DataFrame, list[ProcessingFailure], int]:
    """Generate ensembles and assemble the provenance-rich feature matrix."""
    xyz_dir = output_dir / "xyz"
    conformer_root = output_dir / "conformers"
    xyz_dir.mkdir(parents=True, exist_ok=True)
    conformer_root.mkdir(parents=True, exist_ok=True)
    records: list[dict[str, object]] = []
    failures: list[ProcessingFailure] = []
    nonconverged = 0
    total = len(source)

    for position, row in source.iterrows():
        reaction_id = str(row["Reaction_ID"])
        ligand_smiles = str(row["Ligand_SMILES"])
        try:
            ensemble = embed_and_optimize(
                ligand_smiles,
                base_seed + int(position) * 104_729,
                conformer_count,
                prune_rms,
                energy_window,
                float(row["Temp_K"]),
            )
            molecule = ensemble.molecule
            attach_idx, axis_idx = identify_donor_and_axis(molecule)
            filename = safe_filename(reaction_id)
            destination = xyz_dir / filename
            write_xyz(
                molecule,
                destination,
                reaction_id,
                attach_idx,
                axis_idx,
                ensemble.conformer_ids[0],
                ensemble.relative_energies_kcal_mol[0],
                ensemble.boltzmann_weights[0],
            )

            ensemble_directory = conformer_root / Path(filename).stem
            ensemble_directory.mkdir(parents=True, exist_ok=True)
            conformer_paths: list[str] = []
            for conformer_position, (
                conformer_id,
                relative_energy,
                boltzmann_weight,
            ) in enumerate(
                zip(
                    ensemble.conformer_ids,
                    ensemble.relative_energies_kcal_mol,
                    ensemble.boltzmann_weights,
                    strict=True,
                )
            ):
                conformer_filename = f"conf_{conformer_position:03d}.xyz"
                conformer_destination = ensemble_directory / conformer_filename
                write_xyz(
                    molecule,
                    conformer_destination,
                    reaction_id,
                    attach_idx,
                    axis_idx,
                    conformer_id,
                    relative_energy,
                    boltzmann_weight,
                )
                conformer_paths.append(
                    conformer_destination.relative_to(output_dir).as_posix()
                )

            geometry = GeometryResult(
                xyz_relative_path=(Path("xyz") / filename).as_posix(),
                conformer_relative_paths=tuple(conformer_paths),
                relative_energies_kcal_mol=ensemble.relative_energies_kcal_mol,
                boltzmann_weights=ensemble.boltzmann_weights,
                attach_atom_idx=attach_idx,
                primary_bond_vector_idx=axis_idx,
                mmff_statuses=ensemble.mmff_statuses,
            )
            nonconverged += sum(status > 0 for status in geometry.mmff_statuses)
            scaffold = MurckoScaffold.MurckoScaffoldSmiles(mol=Chem.RemoveHs(molecule))
            if not scaffold:
                scaffold = Chem.MolToSmiles(Chem.RemoveHs(molecule), canonical=True)
            records.append(
                {
                    "Reaction_ID": reaction_id,
                    "Ligand_XYZ_Path": geometry.xyz_relative_path,
                    "Attach_Atom_Idx": geometry.attach_atom_idx,
                    "Primary_Bond_Vector_Idx": geometry.primary_bond_vector_idx,
                    "NBO_Charge": float(row["NBO_Charge_N"]),
                    "IR_Frequency": float(row["IR_Stretching_Freq"]),
                    "Temp_K": float(row["Temp_K"]),
                    "Exp_ddG_kcal_mol": float(row["Experimental_ddG"]),
                    "Conformer_XYZ_Paths": ";".join(geometry.conformer_relative_paths),
                    "Conformer_Relative_Energies_kcal_mol": ";".join(
                        f"{value:.8f}" for value in geometry.relative_energies_kcal_mol
                    ),
                    "Conformer_Boltzmann_Weights": ";".join(
                        f"{value:.10f}" for value in geometry.boltzmann_weights
                    ),
                    "Conformer_Count": len(geometry.conformer_relative_paths),
                    "Ensemble_Energy_Span_kcal_mol": max(
                        geometry.relative_energies_kcal_mol
                    ),
                    "Ligand_SMILES": ligand_smiles,
                    "Ligand_Group": scaffold,
                    "Dataset_Split": str(row["Dataset_Split"]),
                    "Source_ID": str(row["Source_ID"]),
                    "Source_URL": source_name,
                }
            )
            print(
                f"[{position + 1:04d}/{total:04d}] {reaction_id}: "
                f"retained {len(geometry.conformer_relative_paths)} conformer(s), "
                f"ΔE={max(geometry.relative_energies_kcal_mol):.3f} kcal/mol"
            )
        except Exception as exc:
            failure = ProcessingFailure(reaction_id, ligand_smiles, str(exc))
            failures.append(failure)
            print(
                f"[{position + 1:04d}/{total:04d}] {reaction_id}: FAILED - {exc}",
                file=sys.stderr,
            )
            if strict:
                raise RuntimeError(
                    f"strict mode stopped at reaction {reaction_id}: {exc}"
                ) from exc

    result = pd.DataFrame.from_records(records, columns=OUTPUT_COLUMNS)
    return result, failures, nonconverged


def atomic_write_csv(frame: pd.DataFrame, destination: Path) -> None:
    """Write a CSV without exposing a partially written destination."""
    temporary = destination.with_suffix(destination.suffix + ".tmp")
    frame.to_csv(temporary, index=False, float_format="%.8g")
    os.replace(temporary, destination)


def preserve_official_catalog(
    catalog: pd.DataFrame,
    source_bytes: bytes,
    source_url: str,
    output_dir: Path,
) -> None:
    """Persist the complete source table and a reproducibility manifest."""
    official_dir = output_dir / "official"
    official_dir.mkdir(parents=True, exist_ok=True)
    catalog_path = official_dir / "ni_hda_kraken.csv"
    temporary_catalog = catalog_path.with_suffix(".csv.tmp")
    temporary_catalog.write_bytes(source_bytes)
    os.replace(temporary_catalog, catalog_path)
    ddg = pd.to_numeric(catalog.get("ddG_abs"), errors="coerce")
    ee = pd.to_numeric(catalog.get("ee"), errors="coerce")
    manifest = {
        "schema_version": 1,
        "source_url": source_url,
        "retrieved_at_utc": datetime.now(UTC).isoformat(),
        "sha256": hashlib.sha256(source_bytes).hexdigest(),
        "rows": len(catalog),
        "columns": len(catalog.columns),
        "experimental_ddg_rows": int(ddg.notna().sum()),
        "experimental_ee_rows": int(ee.notna().sum()),
        "published_training_source_ids": sorted(
            PUBLISHED_TRAIN_IDS, key=lambda value: int(value)
        ),
        "historical_blind_source_ids": sorted(
            HISTORICAL_BLIND_IDS, key=lambda value: int(value)
        ),
        "scientific_note": (
            "The complete Kraken ligand library is unlabeled except for the "
            "reported experimental subsets; unlabeled entries must not be "
            "treated as measured reaction outcomes."
        ),
    }
    manifest_path = official_dir / "provenance.json"
    temporary = manifest_path.with_suffix(".json.tmp")
    temporary.write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    os.replace(temporary, manifest_path)
    print(f"  Preserved complete official catalog: {catalog_path}")
    print(f"  Provenance manifest: {manifest_path}")


def print_summary(
    result: pd.DataFrame,
    failures: list[ProcessingFailure],
    nonconverged: int,
    source_name: str,
    output_path: Path,
) -> None:
    """Print completion counts and compact numerical summary statistics."""
    print("\nDataset preparation complete")
    print(f"  Source: {source_name}")
    print(f"  Successful geometries: {len(result):,}")
    print(f"  Retained conformers: {int(result['Conformer_Count'].sum()):,}")
    print(f"  Failed rows: {len(failures):,}")
    print(f"  MMFF94 iteration-limit conformers retained: {nonconverged:,}")
    print(f"  Feature matrix: {output_path}")
    if result.empty:
        return
    summary_columns = [
        "NBO_Charge",
        "IR_Frequency",
        "Temp_K",
        "Exp_ddG_kcal_mol",
    ]
    print("\nSummary statistics")
    print(result[summary_columns].describe().round(4).to_string())


def main(argv: Iterable[str] | None = None) -> int:
    args = parse_args(argv)
    RDLogger.DisableLog("rdApp.warning")
    args.output_dir.mkdir(parents=True, exist_ok=True)

    try:
        loaded = load_source_dataset(
            offline=args.offline,
            custom_urls=args.dataset_urls,
            timeout=args.timeout,
        )
        if (
            loaded.complete_catalog is not None
            and loaded.complete_catalog_bytes is not None
        ):
            preserve_official_catalog(
                loaded.complete_catalog,
                loaded.complete_catalog_bytes,
                loaded.source_name,
                args.output_dir,
            )
        source = validate_source_frame(loaded.reactions)
        if args.max_records is not None:
            source = source.head(args.max_records).reset_index(drop=True)
        print(f"Preparing {len(source):,} reaction record(s).")

        result, failures, nonconverged = process_reactions(
            source=source,
            output_dir=args.output_dir,
            base_seed=args.seed,
            strict=args.strict,
            conformer_count=args.conformers,
            prune_rms=args.prune_rms,
            energy_window=args.energy_window,
            source_name=loaded.source_name,
        )
        if result.empty:
            raise RuntimeError(
                "all ligand geometries failed; no feature matrix was written"
            )

        output_path = args.output_dir / "reactions_raw.csv"
        atomic_write_csv(result, output_path)
        failure_path = args.output_dir / "etl_failures.csv"
        if failures:
            failure_frame = pd.DataFrame(
                [
                    {
                        "Reaction_ID": failure.reaction_id,
                        "Ligand_SMILES": failure.ligand_smiles,
                        "Reason": failure.reason,
                    }
                    for failure in failures
                ]
            )
            atomic_write_csv(failure_frame, failure_path)
            print(f"Failure details: {failure_path}", file=sys.stderr)
        elif failure_path.exists():
            failure_path.unlink()

        print_summary(
            result,
            failures,
            nonconverged,
            loaded.source_name,
            output_path,
        )
        return 0
    except (DatasetError, OSError, RuntimeError, ValueError) as exc:
        print(f"Fatal ETL error: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
