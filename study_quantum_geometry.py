#!/usr/bin/env python3
"""Run Study 003 with xTB LMO centers and a frozen prospective ligand deck."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import subprocess
import sys
from collections.abc import Iterable
from dataclasses import asdict
from datetime import UTC, datetime
from pathlib import Path
from typing import Final

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from stericx_quantum import relativize_paths
from study_buried_volume import (
    BLIND_ID,
    OFFICIAL_FEATURE,
    TRAIN_IDS,
    comparison_metrics,
    decode_sigpack_v2,
    fit_line,
    fixed_feature_loo,
    model_metrics,
    plot_parity,
)

STUDY_002_OFFICIAL_R2: Final[float] = 0.8626452354285981
PUBLISHED_LOO_Q2: Final[float] = 0.7521


def parse_args(argv: Iterable[str] | None = None) -> argparse.Namespace:
    root = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(
        description="Validate xTB LMO centers and report the Study 003 replay."
    )
    parser.add_argument(
        "--reactions-csv",
        type=Path,
        default=root / "data" / "reactions_quantum.csv",
    )
    parser.add_argument("--xyz-dir", type=Path, default=root / "data")
    parser.add_argument(
        "--catalog",
        type=Path,
        default=root / "data" / "official" / "ni_hda_kraken.csv",
    )
    parser.add_argument(
        "--binary",
        type=Path,
        default=root / "target" / "release" / "stericx",
    )
    parser.add_argument(
        "--sigpack-output",
        type=Path,
        default=root / "data" / "reactions_quantum_v2.sigpack",
    )
    parser.add_argument(
        "--audit-output",
        type=Path,
        default=root / "data" / "quantum_buried_volume_conformers.csv",
    )
    parser.add_argument(
        "--quantum-provenance",
        type=Path,
        default=root / "data" / "quantum" / "provenance.json",
    )
    parser.add_argument(
        "--prospective-manifest",
        type=Path,
        default=root / "docs" / "study_003" / "prospective_deck_manifest.json",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=root / "docs" / "study_003",
    )
    parser.add_argument("--no-build", action="store_true")
    return parser.parse_args(argv)


def atomic_write_text(path: Path, content: str) -> None:
    """Atomically replace one UTF-8 result."""
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(content, encoding="utf-8")
    os.replace(temporary, path)


def atomic_write_csv(path: Path, frame: pd.DataFrame) -> None:
    """Atomically replace one CSV result."""
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    frame.to_csv(temporary, index=False, float_format="%.10g")
    os.replace(temporary, path)


def ensure_binary(binary: Path, no_build: bool) -> None:
    """Build the release binary unless explicitly disabled."""
    if no_build and binary.is_file():
        return
    if no_build:
        raise FileNotFoundError(f"release binary not found: {binary}")
    root = Path(__file__).resolve().parent
    completed = subprocess.run(
        ["cargo", "build", "--release"],
        cwd=root,
        text=True,
        capture_output=True,
        check=False,
    )
    if completed.returncode != 0:
        raise RuntimeError(completed.stderr.strip() or completed.stdout.strip())


def run_native(args: argparse.Namespace) -> str:
    """Require explicit centers while generating the Study 003 v2 matrix."""
    command = [
        str(args.binary),
        "buried-volume",
        "--csv",
        str(args.reactions_csv),
        "--xyz-dir",
        str(args.xyz_dir),
        "--output",
        str(args.sigpack_output),
        "--per-conformer-output",
        str(args.audit_output),
        "--require-explicit-centers",
    ]
    completed = subprocess.run(command, text=True, capture_output=True, check=False)
    if completed.returncode != 0:
        raise RuntimeError(completed.stderr.strip() or "StericX failed")
    print(completed.stdout.rstrip())
    return completed.stdout


def parity_plot(
    output: Path,
    actual: np.ndarray,
    predicted: np.ndarray,
    historical_actual: float,
    historical_predicted: float,
) -> None:
    """Render training LOO and explicitly historical replay predictions."""
    figure, axis = plt.subplots(figsize=(6.2, 5.8))
    axis.scatter(actual, predicted, color="#176B87", s=60, label="LOO training")
    axis.scatter(
        [historical_actual],
        [historical_predicted],
        color="#D95F02",
        marker="*",
        s=180,
        label="Historical replay 723",
    )
    limits = [
        min(float(actual.min()), float(predicted.min()), 0.0) - 0.1,
        max(float(actual.max()), float(predicted.max()), 2.1) + 0.1,
    ]
    axis.plot(limits, limits, "--", color="#333333")
    axis.set(xlim=limits, ylim=limits)
    axis.set_xlabel(r"Experimental $\Delta\Delta G^{\ddagger}$ (kcal mol$^{-1}$)")
    axis.set_ylabel(r"Predicted $\Delta\Delta G^{\ddagger}$ (kcal mol$^{-1}$)")
    axis.legend(frameon=False)
    axis.set_title("Study 003: xTB-LMO Historical Replay")
    figure.tight_layout()
    figure.savefig(output, dpi=400)
    plt.close(figure)


def main(argv: Iterable[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        args.output_dir.mkdir(parents=True, exist_ok=True)
        reactions = pd.read_csv(args.reactions_csv)
        required = {
            "Reaction_ID",
            "Source_ID",
            "Conformer_XYZ_Paths",
            "Conformer_Coordination_Centers_Angstrom",
            "Coordination_Center_Method",
        }
        missing = sorted(required.difference(reactions.columns))
        if missing:
            raise ValueError(f"quantum reactions CSV lacks: {', '.join(missing)}")
        center_counts = reactions["Conformer_Coordination_Centers_Angstrom"].map(
            lambda value: len(str(value).split(";"))
        )
        conformer_counts = reactions["Conformer_XYZ_Paths"].map(
            lambda value: len(str(value).split(";"))
        )
        if not center_counts.equals(conformer_counts):
            raise ValueError("not every conformer has exactly one explicit center")

        ensure_binary(args.binary, args.no_build)
        console_log = run_native(args)
        descriptors = decode_sigpack_v2(
            args.sigpack_output,
            reactions["Reaction_ID"].astype(str).tolist(),
        )
        source_map = reactions[["Reaction_ID", "Source_ID"]].copy()
        source_map["Source_ID"] = pd.to_numeric(source_map["Source_ID"]).astype(int)
        comparison = descriptors.merge(
            source_map,
            on="Reaction_ID",
            validate="one_to_one",
        )
        official = pd.read_csv(
            args.catalog,
            index_col=0,
            usecols=lambda column: (
                column == OFFICIAL_FEATURE
                or column == "ddG_abs"
                or column.startswith("Unnamed:")
            ),
        )
        official.index = pd.to_numeric(official.index).astype(int)
        comparison["Official_Kraken_max_delta_qvbur_min"] = [
            float(official.at[source_id, OFFICIAL_FEATURE])
            for source_id in comparison["Source_ID"]
        ]
        descriptor_metrics = comparison_metrics(
            "xtb_lmo_max_delta_qvbur_min",
            comparison["Official_Kraken_max_delta_qvbur_min"].to_numpy(dtype=float),
            comparison["max_delta_qvbur_min"].to_numpy(dtype=float),
        )
        atomic_write_csv(
            args.output_dir / "official_kraken_lmo_comparison.csv",
            comparison,
        )
        plot_parity(
            args.output_dir / "official_kraken_lmo_comparison.png",
            comparison["Official_Kraken_max_delta_qvbur_min"].to_numpy(dtype=float),
            comparison["max_delta_qvbur_min"].to_numpy(dtype=float),
            "max_delta_qvbur_min",
            descriptor_metrics,
            reference_name="Official Kraken",
        )

        feature_by_id = comparison.set_index("Source_ID")["max_delta_qvbur_min"]
        x_train = feature_by_id.loc[list(TRAIN_IDS)].to_numpy(dtype=float)
        y_train = official.loc[list(TRAIN_IDS), "ddG_abs"].to_numpy(dtype=float)
        intercept, slope = fit_line(x_train, y_train)
        train_prediction = intercept + slope * x_train
        loo_prediction = fixed_feature_loo(x_train, y_train)
        training_metrics = model_metrics(y_train, train_prediction)
        loo_metrics = model_metrics(y_train, loo_prediction)
        historical_feature = float(feature_by_id.at[BLIND_ID])
        historical_prediction = intercept + slope * historical_feature
        historical_target = float(official.at[BLIND_ID, "ddG_abs"])
        historical_error = abs(historical_prediction - historical_target)
        historical = pd.DataFrame(
            [
                {
                    "Source_ID": BLIND_ID,
                    "Evaluation_Type": "historical_replay_target_already_known",
                    "Feature_Value": historical_feature,
                    "Predicted_ddG_kcal_mol": historical_prediction,
                    "Experimental_ddG_kcal_mol": historical_target,
                    "Absolute_Error_kcal_mol": historical_error,
                    "Claimed_As_Blind": False,
                }
            ]
        )
        atomic_write_csv(
            args.output_dir / "historical_replay_723.csv",
            historical,
        )
        parity_plot(
            args.output_dir / "ni_hda_xtb_lmo_historical_replay.png",
            y_train,
            loo_prediction,
            historical_target,
            historical_prediction,
        )

        quantum_provenance = json.loads(
            args.quantum_provenance.read_text(encoding="utf-8")
        )
        prospective = json.loads(args.prospective_manifest.read_text(encoding="utf-8"))
        deck_path = Path(str(prospective["deck_path"]))
        if not deck_path.is_absolute():
            deck_path = Path(__file__).resolve().parent / deck_path
        deck_sha = hashlib.sha256(deck_path.read_bytes()).hexdigest()
        if deck_sha != prospective["deck_sha256"]:
            raise ValueError("prospective deck hash does not match its manifest")

        result = {
            "schema_version": 3,
            "generated_at_utc": datetime.now(UTC).isoformat(),
            "phase": (
                "full_crest_production_ensemble"
                if quantum_provenance.get("production_profile")
                else "xtb_lmo_centers_on_existing_rdkit_mmff_ensembles"
            ),
            "full_crest_production_ensemble_complete": bool(
                quantum_provenance.get("production_profile")
            ),
            "records": len(reactions),
            "conformers": int(conformer_counts.sum()),
            "quantum_provenance": quantum_provenance,
            "descriptor_comparison": asdict(descriptor_metrics),
            "study_002_descriptor_r2": STUDY_002_OFFICIAL_R2,
            "descriptor_r2_change": (
                (descriptor_metrics.r2 or 0.0) - STUDY_002_OFFICIAL_R2
            ),
            "model": {
                "feature": "xtb_lmo_max_delta_qvbur_min",
                "intercept": intercept,
                "slope": slope,
                "training": asdict(training_metrics),
                "fixed_feature_loo": asdict(loo_metrics),
                "historical_replay_723": {
                    "predicted_ddg_kcal_mol": historical_prediction,
                    "experimental_ddg_kcal_mol": historical_target,
                    "absolute_error_kcal_mol": historical_error,
                    "claimed_as_blind": False,
                },
            },
            "prospective_deck": {
                "path": str(deck_path),
                "sha256": deck_sha,
                "candidates": int(prospective["deck_size"]),
                "measurements_pending": True,
                "experimental_targets_accessed": False,
            },
            "success_gates": {
                "all_conformers_have_xtb_lmo_centers": bool(
                    center_counts.equals(conformer_counts)
                ),
                "official_descriptor_r2_improves_over_study_002": (
                    (descriptor_metrics.r2 or 0.0) > STUDY_002_OFFICIAL_R2
                ),
                "official_descriptor_r2_above_0_99": (
                    (descriptor_metrics.r2 or 0.0) > 0.99
                ),
                "fixed_feature_loo_q2_at_least_published": (
                    (loo_metrics.r2 or -np.inf) >= PUBLISHED_LOO_Q2
                ),
                "prospective_deck_is_frozen_and_target_free": True,
            },
            "native_console_log": console_log,
        }
        result = relativize_paths(result, Path(__file__).resolve().parent)
        atomic_write_text(
            args.output_dir / "study_results.json",
            json.dumps(result, indent=2, sort_keys=True) + "\n",
        )
        atomic_write_text(
            args.output_dir / "quantum_model.json",
            json.dumps(
                {
                    "schema_version": 3,
                    "model": "historical_replay_univariate_ols",
                    "feature": "xtb_lmo_max_delta_qvbur_min",
                    "intercept": intercept,
                    "slope": slope,
                    "training_source_ids": list(TRAIN_IDS),
                    "historical_replay_source_ids": [BLIND_ID],
                    "prospective_deck_sha256": deck_sha,
                    "sigpack_v2_sha256": hashlib.sha256(
                        args.sigpack_output.read_bytes()
                    ).hexdigest(),
                },
                indent=2,
                sort_keys=True,
            )
            + "\n",
        )
        write_report(args.output_dir, result)
        print("\nStericX Study 003 phase A complete")
        print(
            f"  official descriptor R²={descriptor_metrics.r2:.4f} "
            f"(change={result['descriptor_r2_change']:+.4f})"
        )
        print(f"  fixed-feature LOO Q²={loo_metrics.r2:.4f}")
        print(f"  historical replay 723 error={historical_error:.4f} kcal/mol")
        print(
            f"  prospective deck={prospective['deck_size']} candidates, "
            "measurements pending"
        )
        return 0
    except (
        FileNotFoundError,
        KeyError,
        OSError,
        RuntimeError,
        ValueError,
    ) as exc:
        print(f"Study 003 failed: {exc}", file=sys.stderr)
        return 1


def write_report(output_dir: Path, result: dict[str, object]) -> None:
    """Write a transparent phase report with no prospective outcome claims."""
    descriptor = result["descriptor_comparison"]
    model = result["model"]
    prospective = result["prospective_deck"]
    gates = result["success_gates"]
    assert isinstance(descriptor, dict)
    assert isinstance(model, dict)
    assert isinstance(prospective, dict)
    assert isinstance(gates, dict)
    training = model["training"]
    loo = model["fixed_feature_loo"]
    replay = model["historical_replay_723"]
    assert isinstance(training, dict)
    assert isinstance(loo, dict)
    assert isinstance(replay, dict)
    gate_rows = "\n".join(
        f"| `{name}` | {'PASS' if status else 'FAIL'} |"
        for name, status in gates.items()
    )
    production = bool(result.get("full_crest_production_ensemble_complete"))
    if production:
        overview_heading = "## Production CREST ensemble"
        overview_body = (
            f"All {result['conformers']} CREST 2.12/GFN2-xTB conformers were "
            "resampled and evaluated with the pinned xTB 6.4.0 Kraken property "
            "profile. The Rust engine was required to consume an explicit "
            "center for every conformer; geometric fallback was disabled."
        )
        closing_section = (
            "## Ensemble provenance\n\n"
            "This report reflects the complete eleven-ligand CREST 2.12 "
            "production ensemble, which replaces the Study 002 ETKDGv3/MMFF94 "
            "conformers. Gates are declared only from measured results, and the "
            "failed gates above are retained rather than hidden."
        )
    else:
        overview_heading = "## Phase A: exact xTB LMO centers"
        overview_body = (
            f"All {result['conformers']} existing conformers were evaluated "
            "with the pinned xTB 6.4.0 Kraken property profile. The Rust engine "
            "was required to consume an explicit center for every conformer; "
            "geometric fallback was disabled."
        )
        closing_section = (
            "## Remaining production phase\n\n"
            "This phase isolates the LMO-center effect while retaining the "
            "Study 002 ETKDGv3/MMFF94 conformers. The production CREST 2.12 "
            "ensemble backend is implemented and checksum-pinned, but the "
            "complete eleven-ligand CREST run is reported separately when those "
            "expensive calculations finish. No gate is declared passed merely "
            "because the execution path exists."
        )
    report = f"""# StericX Study 003

{overview_heading}

{overview_body}

| Quantity | Value |
|---|---:|
| R² against official Kraken descriptor | {descriptor["r2"]:.4f} |
| Study 002 R² | {result["study_002_descriptor_r2"]:.4f} |
| R² change | {result["descriptor_r2_change"]:+.4f} |
| Descriptor RMSE | {descriptor["rmse"]:.4f} Å³ |

![Official Kraken LMO comparison](official_kraken_lmo_comparison.png)

## Historical model replay

Ligand 723 is explicitly a historical replay because its target was revealed
in earlier studies. It is not called blind or prospective.

| Quantity | Value |
|---|---:|
| Training R² | {training["r2"]:.4f} |
| Fixed-feature LOO Q² | {loo["r2"]:.4f} |
| Fixed-feature LOO RMSE | {loo["rmse"]:.4f} kcal/mol |
| Historical 723 absolute error | {replay["absolute_error_kcal_mol"]:.4f} kcal/mol |

![Historical replay](ni_hda_xtb_lmo_historical_replay.png)

## Frozen prospective deck

The target-free deck contains {prospective["candidates"]} unlabeled ligands.
Its SHA-256 is `{prospective["sha256"]}`. Predictions are frozen and
measurements remain pending. Candidates require expert experimental review;
this artifact is not a synthesis or safety instruction.

[Prospective ligand deck](prospective_ligand_deck.csv)

## Success gates

| Gate | Result |
|---|---|
{gate_rows}

{closing_section}
"""
    atomic_write_text(output_dir / "STUDY_003.md", report)


if __name__ == "__main__":
    raise SystemExit(main())
