#!/usr/bin/env python3
"""Freeze a diverse, target-free Ni-hDA ligand deck for prospective testing.

The deck is a computational recommendation for expert experimental triage, not
an instruction to perform chemistry. No experimental target is joined or
consulted for any selected ligand.
"""

from __future__ import annotations

import argparse
import hashlib
import io
import json
import os
import sys
from collections.abc import Iterable
from datetime import UTC, datetime
from pathlib import Path
from typing import Final

import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem.Scaffolds import MurckoScaffold

TRAIN_IDS: Final[tuple[int, ...]] = (
    401,
    498,
    724,
    785,
    1057,
    1058,
    2062,
    2063,
    2064,
    2067,
)
HISTORICAL_IDS: Final[frozenset[int]] = frozenset({*TRAIN_IDS, 723})
FEATURE: Final[str] = "vbur_max_delta_qvbur_min"
DIVERSITY_FEATURES: Final[tuple[str, ...]] = (
    FEATURE,
    "nbo_P_boltz",
    "fmo_e_homo_boltz",
    "fmo_e_lumo_boltz",
    "dipolemoment_boltz",
    "vbur_near_vbur_boltz",
    "vbur_far_vbur_boltz",
    "pyr_P_vburminconf",
)
GAS_CONSTANT_KCAL_MOL_K: Final[float] = 0.00198720425864083
TEMPERATURE_K: Final[float] = 298.15


def parse_args(argv: Iterable[str] | None = None) -> argparse.Namespace:
    root = Path(__file__).resolve().parent.parent
    parser = argparse.ArgumentParser(
        description="Freeze a target-free diverse ligand deck for new experiments."
    )
    parser.add_argument(
        "--catalog",
        type=Path,
        default=root / "data" / "official" / "ni_hda_kraken.csv",
    )
    parser.add_argument(
        "--model",
        type=Path,
        default=root / "docs" / "study_001" / "published_model.json",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=root / "docs" / "study_003" / "prospective_ligand_deck.csv",
    )
    parser.add_argument(
        "--manifest",
        type=Path,
        default=root / "docs" / "study_003" / "prospective_deck_manifest.json",
    )
    parser.add_argument(
        "--interpolation-count",
        type=int,
        default=8,
    )
    args = parser.parse_args(argv)
    if args.interpolation_count < 4 or args.interpolation_count > 10:
        parser.error("--interpolation-count must be between 4 and 10")
    return args


def scaffold_smiles(smiles: str) -> str:
    """Return a deterministic scaffold label used only for deck diversity."""
    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None:
        return ""
    scaffold = MurckoScaffold.MurckoScaffoldSmiles(mol=molecule)
    return scaffold or Chem.MolToSmiles(molecule, canonical=True)


def robust_scaler(values: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Learn a robust library median and interquartile scale."""
    median = np.median(values, axis=0)
    first, third = np.quantile(values, [0.25, 0.75], axis=0)
    scale = third - first
    scale[scale <= np.finfo(float).eps] = 1.0
    return median, scale


def greedy_bin_selection(
    candidates: pd.DataFrame,
    standardized: np.ndarray,
    training_standardized: np.ndarray,
    count: int,
) -> list[int]:
    """Select across prediction quantiles with scaffold-aware maximin diversity."""
    prediction = candidates["Predicted_ddG_kcal_mol"].to_numpy(dtype=float)
    quantile_edges = np.quantile(prediction, np.linspace(0.0, 1.0, 5))
    selected: list[int] = []
    used_scaffolds: set[str] = set()
    reference = training_standardized.copy()
    target_per_bin = [count // 4 + int(index < count % 4) for index in range(4)]
    for bin_index, bin_target in enumerate(target_per_bin):
        lower = quantile_edges[bin_index]
        upper = quantile_edges[bin_index + 1]
        mask = (prediction >= lower) & (
            prediction <= upper if bin_index == 3 else prediction < upper
        )
        pool = np.flatnonzero(mask).tolist()
        for _ in range(bin_target):
            scored: list[tuple[float, int, int]] = []
            for position in pool:
                if position in selected:
                    continue
                scaffold = str(candidates.iloc[position]["Scaffold"])
                if not scaffold or scaffold in used_scaffolds:
                    continue
                distance = float(
                    np.linalg.norm(reference - standardized[position], axis=1).min()
                )
                source_id = int(candidates.iloc[position]["Source_ID"])
                scored.append((distance, -source_id, position))
            if not scored:
                raise ValueError(
                    "could not satisfy scaffold diversity in prediction "
                    f"bin {bin_index}"
                )
            _, _, choice = max(scored)
            selected.append(choice)
            used_scaffolds.add(str(candidates.iloc[choice]["Scaffold"]))
            reference = np.vstack([reference, standardized[choice]])
    return selected


def boundary_candidate(
    candidates: pd.DataFrame,
    feature_limit: float,
    direction: str,
    used_scaffolds: set[str],
) -> int:
    """Choose the nearest unlabeled descriptor just outside one training bound."""
    feature = candidates[FEATURE].to_numpy(dtype=float)
    if direction == "lower":
        eligible = np.flatnonzero(feature < feature_limit)
        order = eligible[np.argsort(feature_limit - feature[eligible])]
    else:
        eligible = np.flatnonzero(feature > feature_limit)
        order = eligible[np.argsort(feature[eligible] - feature_limit)]
    for position in order:
        scaffold = str(candidates.iloc[position]["Scaffold"])
        if scaffold and scaffold not in used_scaffolds:
            return int(position)
    raise ValueError(f"no unique-scaffold {direction} boundary candidate exists")


def frozen_csv_bytes(frame: pd.DataFrame) -> bytes:
    """Serialize the target-free deck deterministically."""
    buffer = io.StringIO()
    frame.to_csv(buffer, index=False, float_format="%.10g", lineterminator="\n")
    return buffer.getvalue().encode()


def write_frozen(path: Path, content: bytes) -> None:
    """Write once, accepting only byte-identical reruns."""
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.is_file():
        if path.read_bytes() != content:
            raise FileExistsError(
                f"frozen deck already exists with different content: {path}"
            )
        return
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_bytes(content)
    os.replace(temporary, path)


def main(argv: Iterable[str] | None = None) -> int:
    args = parse_args(argv)
    RDLogger.DisableLog("rdApp.*")
    try:
        catalog_bytes = args.catalog.read_bytes()
        catalog = pd.read_csv(args.catalog, index_col=0)
        catalog.index = pd.to_numeric(catalog.index).astype(int)
        required = {"smiles", "ddG_abs", *DIVERSITY_FEATURES}
        missing = sorted(required.difference(catalog.columns))
        if missing:
            raise ValueError(f"catalog lacks: {', '.join(missing)}")
        model_bytes = args.model.read_bytes()
        model = json.loads(model_bytes)
        if model.get("feature") != FEATURE:
            raise ValueError(
                "model does not use the preregistered buried-volume feature"
            )
        intercept = float(model["intercept"])
        slope = float(model["slope"])
        training_min = float(model["training_feature_minimum"])
        training_max = float(model["training_feature_maximum"])

        unlabeled = catalog.loc[catalog["ddG_abs"].isna()].copy()
        unlabeled = unlabeled.loc[~unlabeled.index.isin(HISTORICAL_IDS)]
        unlabeled = unlabeled.dropna(subset=["smiles", *DIVERSITY_FEATURES])
        unlabeled.insert(0, "Source_ID", unlabeled.index.astype(int))
        unlabeled["Scaffold"] = unlabeled["smiles"].map(scaffold_smiles)
        unlabeled = unlabeled.loc[unlabeled["Scaffold"] != ""].copy()
        unlabeled["Predicted_ddG_kcal_mol"] = intercept + slope * unlabeled[
            FEATURE
        ].to_numpy(dtype=float)
        ratio = np.exp(
            unlabeled["Predicted_ddG_kcal_mol"].to_numpy(dtype=float)
            / (GAS_CONSTANT_KCAL_MOL_K * TEMPERATURE_K)
        )
        unlabeled["Predicted_ee_percent"] = 100.0 * (ratio - 1.0) / (ratio + 1.0)
        unlabeled["Applicability_Domain"] = np.where(
            unlabeled[FEATURE].between(training_min, training_max),
            "inside_training_feature_range",
            np.where(
                unlabeled[FEATURE] < training_min,
                "below_training_feature_range",
                "above_training_feature_range",
            ),
        )

        diversity_values = unlabeled[list(DIVERSITY_FEATURES)].to_numpy(dtype=float)
        training_values = catalog.loc[
            list(TRAIN_IDS), list(DIVERSITY_FEATURES)
        ].to_numpy(dtype=float)
        median, scale = robust_scaler(diversity_values)
        candidate_standardized = (diversity_values - median) / scale
        training_standardized = (training_values - median) / scale

        interpolation = unlabeled.loc[
            unlabeled["Applicability_Domain"] == "inside_training_feature_range"
        ].copy()
        interpolation_positions = unlabeled.index.get_indexer(interpolation.index)
        selected_local = greedy_bin_selection(
            interpolation,
            candidate_standardized[interpolation_positions],
            training_standardized,
            args.interpolation_count,
        )
        selected_indices = [
            int(interpolation_positions[position]) for position in selected_local
        ]
        used_scaffolds = {
            str(unlabeled.iloc[position]["Scaffold"]) for position in selected_indices
        }
        lower = boundary_candidate(
            unlabeled,
            training_min,
            "lower",
            used_scaffolds,
        )
        used_scaffolds.add(str(unlabeled.iloc[lower]["Scaffold"]))
        upper = boundary_candidate(
            unlabeled,
            training_max,
            "upper",
            used_scaffolds,
        )
        selected_indices.extend([lower, upper])

        deck = unlabeled.iloc[selected_indices].copy()
        deck["Selection_Stratum"] = [
            *["interpolation_maximin"] * args.interpolation_count,
            "lower_boundary_challenge",
            "upper_boundary_challenge",
        ]
        deck["Experimental_Target_Accessed"] = False
        deck["Measurement_Status"] = "pending"
        deck["Model_Frozen"] = args.model.name
        deck = deck[
            [
                "Source_ID",
                "smiles",
                "Scaffold",
                FEATURE,
                "Predicted_ddG_kcal_mol",
                "Predicted_ee_percent",
                "Applicability_Domain",
                "Selection_Stratum",
                "Experimental_Target_Accessed",
                "Measurement_Status",
                "Model_Frozen",
            ]
        ].rename(columns={"smiles": "Ligand_SMILES"})
        content = frozen_csv_bytes(deck)
        write_frozen(args.output, content)
        deck_sha = hashlib.sha256(content).hexdigest()
        manifest = {
            "schema_version": 1,
            "frozen_at_utc": datetime.now(UTC).isoformat(),
            "status": "predictions_frozen_measurements_pending",
            "deck_path": str(args.output),
            "deck_sha256": deck_sha,
            "deck_size": len(deck),
            "interpolation_candidates": args.interpolation_count,
            "boundary_challenges": 2,
            "source_catalog": str(args.catalog),
            "source_catalog_sha256": hashlib.sha256(catalog_bytes).hexdigest(),
            "model_path": str(args.model),
            "model_sha256": hashlib.sha256(model_bytes).hexdigest(),
            "training_source_ids": list(TRAIN_IDS),
            "excluded_historical_source_ids": sorted(HISTORICAL_IDS),
            "experimental_targets_accessed_for_candidates": False,
            "selection": {
                "primary_feature": FEATURE,
                "diversity_features": list(DIVERSITY_FEATURES),
                "scaling": "complete-candidate median and interquartile range",
                "interpolation": (
                    "four prediction quantile bins; scaffold-unique greedy "
                    "maximin distance from training and already selected ligands"
                ),
                "boundary": (
                    "nearest unique-scaffold unlabeled ligand immediately "
                    "outside each training feature bound"
                ),
            },
            "limitations": [
                (
                    "Candidates require expert review for availability, "
                    "stability, and safety."
                ),
                (
                    "Predictions are historical-model extrapolations, "
                    "not measured outcomes."
                ),
                "The deck must not be refit after outcomes are known.",
            ],
        }
        args.manifest.parent.mkdir(parents=True, exist_ok=True)
        if args.manifest.is_file():
            existing = json.loads(args.manifest.read_text(encoding="utf-8"))
            if existing.get("deck_sha256") != deck_sha:
                raise FileExistsError(
                    f"frozen manifest refers to a different deck: {args.manifest}"
                )
        else:
            temporary = args.manifest.with_name(args.manifest.name + ".tmp")
            temporary.write_text(
                json.dumps(manifest, indent=2, sort_keys=True) + "\n",
                encoding="utf-8",
            )
            os.replace(temporary, args.manifest)
        print("Prospective ligand deck frozen")
        print(f"  candidates={len(deck)}")
        print(f"  deck_sha256={deck_sha}")
        print(f"  output={args.output}")
        print(f"  manifest={args.manifest}")
        return 0
    except (
        FileExistsError,
        FileNotFoundError,
        KeyError,
        OSError,
        ValueError,
    ) as exc:
        print(f"Prospective deck failed: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
