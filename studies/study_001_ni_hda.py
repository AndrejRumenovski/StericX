#!/usr/bin/env python3
"""Reproduce and stress-test the published Ni-hDA enantioselectivity model.

This workflow deliberately separates model fitting from target revelation:

1. Download or load the complete official 1,566-ligand Kraken table.
2. Fit the published ten-ligand, one-descriptor model.
3. Write and hash the frozen prediction for source ligand 723.
4. Only then reveal its reported ΔΔG‡ and create the scored evaluation.
5. Compare OLS with nested-LOO ridge/LASSO baselines and emit a model card.

Dependencies:
    numpy, pandas, scipy, scikit-learn, matplotlib, requests
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
from collections.abc import Iterable
from dataclasses import asdict, dataclass
from datetime import UTC, datetime
from pathlib import Path
from typing import Final

try:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd
    import requests
    from scipy import stats
    from sklearn.linear_model import Lasso, LinearRegression, Ridge
    from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
    from sklearn.preprocessing import StandardScaler
except ImportError as exc:  # pragma: no cover - environment dependent
    raise SystemExit(
        "Missing study dependency. Install numpy, pandas, scipy, "
        "scikit-learn, matplotlib, and requests."
    ) from exc


SOURCE_URL: Final[str] = (
    "https://raw.githubusercontent.com/SigmanGroup/"
    "Ni-Catalyzed-hDA/main/data/kraken.csv"
)
PUBLISHED_NOTEBOOK_URL: Final[str] = (
    "https://github.com/SigmanGroup/Ni-Catalyzed-hDA/"
    "blob/main/Enantioselectivity_Model.ipynb"
)
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
BLIND_IDS: Final[tuple[int, ...]] = (723,)
PUBLISHED_FEATURE: Final[str] = "vbur_max_delta_qvbur_min"
CURATED_FEATURES: Final[tuple[str, ...]] = (
    PUBLISHED_FEATURE,
    "nbo_P_boltz",
    "fmo_e_homo_boltz",
    "fmo_e_lumo_boltz",
    "dipolemoment_boltz",
    "vbur_near_vbur_boltz",
    "vbur_far_vbur_boltz",
    "pyr_P_vburminconf",
)
RIDGE_GRID: Final[np.ndarray] = np.logspace(-4, 3, 20)
LASSO_GRID: Final[np.ndarray] = np.logspace(-4, 0, 20)


@dataclass(frozen=True)
class Metrics:
    """Regression accuracy for one fixed partition."""

    count: int
    r2: float | None
    mae: float
    rmse: float


@dataclass(frozen=True)
class BaselineResult:
    """Nested-LOO result and final full-training fit."""

    model: str
    alpha: float
    metrics: Metrics
    predictions: list[float]
    coefficients_standardized: dict[str, float]


def parse_args(argv: Iterable[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run the reproducible StericX Ni-hDA model study."
    )
    parser.add_argument(
        "--catalog",
        type=Path,
        default=Path("data/official/ni_hda_kraken.csv"),
        help="Cached official Kraken CSV (downloaded when absent).",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("docs/study_001"),
        help="Study artifact directory (default: docs/study_001).",
    )
    parser.add_argument(
        "--bootstrap",
        type=int,
        default=2_000,
        help="Bootstrap replicates (default: 2000).",
    )
    parser.add_argument(
        "--permutations",
        type=int,
        default=2_000,
        help="Response permutations (default: 2000).",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=20260725,
        help="Deterministic resampling seed.",
    )
    parser.add_argument(
        "--offline",
        action="store_true",
        help="Require the cached catalog and disable network access.",
    )
    args = parser.parse_args(argv)
    if args.bootstrap < 100:
        parser.error("--bootstrap must be at least 100")
    if args.permutations < 100:
        parser.error("--permutations must be at least 100")
    return args


def atomic_write_text(path: Path, content: str) -> None:
    """Atomically replace one UTF-8 text artifact."""
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(content, encoding="utf-8")
    os.replace(temporary, path)


def atomic_write_csv(path: Path, frame: pd.DataFrame) -> None:
    """Atomically replace one CSV artifact."""
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    frame.to_csv(temporary, index=False, float_format="%.10g")
    os.replace(temporary, path)


def load_catalog(path: Path, offline: bool) -> tuple[pd.DataFrame, str]:
    """Load the complete official table and return its canonical SHA-256."""
    if not path.is_file():
        if offline:
            raise FileNotFoundError(f"offline catalog not found: {path}")
        print(f"Downloading official Kraken table: {SOURCE_URL}")
        response = requests.get(
            SOURCE_URL,
            timeout=(10.0, 60.0),
            headers={"User-Agent": "StericX-study-001/1.0"},
        )
        response.raise_for_status()
        path.parent.mkdir(parents=True, exist_ok=True)
        temporary = path.with_name(path.name + ".tmp")
        temporary.write_bytes(response.content)
        os.replace(temporary, path)

    raw = path.read_bytes()
    frame = pd.read_csv(path, index_col=0)
    frame.index = pd.to_numeric(frame.index, errors="raise").astype(int)
    if len(frame) != 1_566:
        raise ValueError(f"expected 1,566 official Kraken rows, found {len(frame)}")
    required = {"smiles", "ddG_abs", *CURATED_FEATURES}
    missing = sorted(required.difference(frame.columns))
    if missing:
        raise ValueError(f"catalog lacks required columns: {', '.join(missing)}")
    return frame, hashlib.sha256(raw).hexdigest()


def calculate_metrics(actual: np.ndarray, predicted: np.ndarray) -> Metrics:
    """Return finite metrics, omitting R² for a one-point partition."""
    return Metrics(
        count=int(actual.size),
        r2=float(r2_score(actual, predicted)) if actual.size >= 2 else None,
        mae=float(mean_absolute_error(actual, predicted)),
        rmse=float(np.sqrt(mean_squared_error(actual, predicted))),
    )


def fit_published_model(
    frame: pd.DataFrame,
) -> tuple[LinearRegression, np.ndarray, np.ndarray, Metrics]:
    """Fit exactly the univariate model defined by the published notebook."""
    train = frame.loc[list(TRAIN_IDS)]
    x_train = train[[PUBLISHED_FEATURE]].to_numpy(dtype=float)
    y_train = train["ddG_abs"].to_numpy(dtype=float)
    if not np.isfinite(x_train).all() or not np.isfinite(y_train).all():
        raise ValueError("published training rows contain missing model values")
    model = LinearRegression().fit(x_train, y_train)
    prediction = model.predict(x_train)
    return model, x_train, y_train, calculate_metrics(y_train, prediction)


def fixed_feature_loo(
    x_train: np.ndarray,
    y_train: np.ndarray,
) -> tuple[np.ndarray, Metrics]:
    """Leave each reaction out while retaining the preregistered descriptor."""
    predicted = np.empty_like(y_train)
    for held_out in range(y_train.size):
        mask = np.arange(y_train.size) != held_out
        model = LinearRegression().fit(x_train[mask], y_train[mask])
        predicted[held_out] = model.predict(x_train[[held_out]])[0]
    return predicted, calculate_metrics(y_train, predicted)


def nested_regularized_loo(
    x_train: np.ndarray,
    y_train: np.ndarray,
    feature_names: tuple[str, ...],
    kind: str,
) -> BaselineResult:
    """Tune regularization entirely inside each outer LOO training fold."""
    alpha_grid = RIDGE_GRID if kind == "ridge" else LASSO_GRID
    outer_predictions = np.empty_like(y_train)
    for held_out in range(y_train.size):
        outer_mask = np.arange(y_train.size) != held_out
        x_outer = x_train[outer_mask]
        y_outer = y_train[outer_mask]
        best_alpha = choose_alpha_nested(x_outer, y_outer, alpha_grid, kind)
        scaler = StandardScaler().fit(x_outer)
        model = regularized_model(kind, best_alpha)
        model.fit(scaler.transform(x_outer), y_outer)
        outer_predictions[held_out] = model.predict(
            scaler.transform(x_train[[held_out]])
        )[0]

    final_alpha = choose_alpha_nested(x_train, y_train, alpha_grid, kind)
    scaler = StandardScaler().fit(x_train)
    final_model = regularized_model(kind, final_alpha)
    final_model.fit(scaler.transform(x_train), y_train)
    coefficients = {
        feature: float(coefficient)
        for feature, coefficient in zip(
            feature_names,
            final_model.coef_,
            strict=True,
        )
    }
    return BaselineResult(
        model=kind,
        alpha=float(final_alpha),
        metrics=calculate_metrics(y_train, outer_predictions),
        predictions=[float(value) for value in outer_predictions],
        coefficients_standardized=coefficients,
    )


def choose_alpha_nested(
    x_train: np.ndarray,
    y_train: np.ndarray,
    alpha_grid: np.ndarray,
    kind: str,
) -> float:
    """Choose alpha by inner LOO without leaking the held-out outer row."""
    scored: list[tuple[float, float]] = []
    for alpha in alpha_grid:
        predicted = np.empty_like(y_train)
        for held_out in range(y_train.size):
            mask = np.arange(y_train.size) != held_out
            scaler = StandardScaler().fit(x_train[mask])
            model = regularized_model(kind, float(alpha))
            model.fit(scaler.transform(x_train[mask]), y_train[mask])
            predicted[held_out] = model.predict(scaler.transform(x_train[[held_out]]))[
                0
            ]
        scored.append((float(mean_squared_error(y_train, predicted)), float(alpha)))
    return min(scored, key=lambda item: (item[0], item[1]))[1]


def regularized_model(kind: str, alpha: float) -> Ridge | Lasso:
    """Construct one deterministic regularized linear estimator."""
    if kind == "ridge":
        return Ridge(alpha=alpha)
    return Lasso(alpha=alpha, max_iter=500_000, tol=1.0e-7, selection="cyclic")


def bootstrap_coefficients(
    x_train: np.ndarray,
    y_train: np.ndarray,
    samples: int,
    seed: int,
) -> dict[str, list[float]]:
    """Bootstrap the preregistered univariate slope and intercept."""
    rng = np.random.default_rng(seed)
    coefficients: list[tuple[float, float]] = []
    for _ in range(samples):
        indices = rng.integers(0, y_train.size, size=y_train.size)
        x_sample = x_train[indices]
        if np.ptp(x_sample[:, 0]) <= np.finfo(float).eps:
            continue
        model = LinearRegression().fit(x_sample, y_train[indices])
        coefficients.append((float(model.intercept_), float(model.coef_[0])))
    values = np.asarray(coefficients)
    if values.shape[0] < samples // 2:
        raise RuntimeError("too few non-degenerate bootstrap samples")
    return {
        "intercept_95_ci": [
            float(value) for value in np.quantile(values[:, 0], [0.025, 0.975])
        ],
        "slope_95_ci": [
            float(value) for value in np.quantile(values[:, 1], [0.025, 0.975])
        ],
    }


def response_permutation_test(
    x_train: np.ndarray,
    y_train: np.ndarray,
    observed_r2: float,
    samples: int,
    seed: int,
) -> float:
    """Test the fixed descriptor against shuffled reaction outcomes."""
    rng = np.random.default_rng(seed)
    extreme = 0
    for _ in range(samples):
        permuted = rng.permutation(y_train)
        model = LinearRegression().fit(x_train, permuted)
        if r2_score(permuted, model.predict(x_train)) >= observed_r2:
            extreme += 1
    return float((extreme + 1) / (samples + 1))


def descriptor_diagnostics(
    x_train: np.ndarray,
    feature_names: tuple[str, ...],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Calculate descriptor correlations and classical VIF values."""
    frame = pd.DataFrame(x_train, columns=feature_names)
    correlations = frame.corr()
    vif_rows: list[dict[str, float | str]] = []
    for feature in feature_names:
        others = [column for column in feature_names if column != feature]
        model = LinearRegression().fit(frame[others], frame[feature])
        r2 = float(model.score(frame[others], frame[feature]))
        vif_rows.append(
            {
                "feature": feature,
                "vif": float(1.0 / max(1.0 - r2, np.finfo(float).eps)),
            }
        )
    return correlations, pd.DataFrame(vif_rows)


def freeze_blind_prediction(
    frame: pd.DataFrame,
    model: LinearRegression,
    output_dir: Path,
) -> tuple[pd.DataFrame, str]:
    """Write the hidden-target prediction before any scoring occurs."""
    blind = frame.loc[list(BLIND_IDS)]
    values = blind[[PUBLISHED_FEATURE]].to_numpy(dtype=float)
    train_values = frame.loc[list(TRAIN_IDS), PUBLISHED_FEATURE].to_numpy(dtype=float)
    status = np.where(
        (values[:, 0] >= train_values.min()) & (values[:, 0] <= train_values.max()),
        "inside_training_range",
        "outside_training_range",
    )
    frozen = pd.DataFrame(
        {
            "Source_ID": list(BLIND_IDS),
            "Dataset_Split": "historical_blind",
            "Feature": PUBLISHED_FEATURE,
            "Feature_Value": values[:, 0],
            "Predicted_ddG_kcal_mol": model.predict(values),
            "Applicability_Domain": status,
            "Target_Accessed_During_Fit": False,
        }
    )
    path = output_dir / "frozen_predictions.csv"
    atomic_write_csv(path, frozen)
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    return frozen, digest


def reveal_and_score(
    frame: pd.DataFrame,
    frozen: pd.DataFrame,
) -> tuple[pd.DataFrame, Metrics]:
    """Join experimental targets only after the frozen artifact exists."""
    revealed = frozen.copy()
    revealed["Experimental_ddG_kcal_mol"] = [
        float(frame.at[int(source_id), "ddG_abs"])
        for source_id in revealed["Source_ID"]
    ]
    revealed["Residual_kcal_mol"] = (
        revealed["Predicted_ddG_kcal_mol"] - revealed["Experimental_ddG_kcal_mol"]
    )
    metrics = calculate_metrics(
        revealed["Experimental_ddG_kcal_mol"].to_numpy(dtype=float),
        revealed["Predicted_ddG_kcal_mol"].to_numpy(dtype=float),
    )
    return revealed, metrics


def plot_results(
    output_dir: Path,
    x_train: np.ndarray,
    y_train: np.ndarray,
    model: LinearRegression,
    loo_predictions: np.ndarray,
    revealed: pd.DataFrame,
) -> None:
    """Render publication-style parity, relationship, and residual figures."""
    output_dir.mkdir(parents=True, exist_ok=True)
    figure, axis = plt.subplots(figsize=(6.4, 6.0))
    axis.scatter(
        y_train,
        loo_predictions,
        color="#176B87",
        s=60,
        edgecolor="white",
        linewidth=0.7,
        label="Training: fixed-feature LOO",
    )
    axis.scatter(
        revealed["Experimental_ddG_kcal_mol"],
        revealed["Predicted_ddG_kcal_mol"],
        color="#D95F02",
        marker="*",
        s=190,
        label="Historical blind",
        zorder=4,
    )
    limits = [
        min(y_train.min(), loo_predictions.min(), 0.0) - 0.1,
        max(y_train.max(), loo_predictions.max(), 2.2) + 0.1,
    ]
    axis.plot(limits, limits, linestyle="--", color="#333333", linewidth=1.2)
    axis.set(xlim=limits, ylim=limits)
    axis.set_xlabel(r"Experimental $\Delta\Delta G^{\ddagger}$ (kcal mol$^{-1}$)")
    axis.set_ylabel(r"Predicted $\Delta\Delta G^{\ddagger}$ (kcal mol$^{-1}$)")
    axis.legend(frameon=False)
    axis.set_title("StericX Study 001: Held-out Prediction")
    figure.tight_layout()
    figure.savefig(output_dir / "ni_hda_parity.png", dpi=400)
    plt.close(figure)

    figure, axis = plt.subplots(figsize=(6.4, 5.5))
    axis.scatter(x_train[:, 0], y_train, color="#176B87", s=60)
    grid = np.linspace(x_train[:, 0].min(), x_train[:, 0].max(), 300)
    axis.plot(grid, model.predict(grid[:, None]), color="#D95F02", linewidth=2)
    axis.scatter(
        revealed["Feature_Value"],
        revealed["Experimental_ddG_kcal_mol"],
        color="#7A1FA2",
        marker="*",
        s=190,
        label="Historical blind",
    )
    axis.set_xlabel(PUBLISHED_FEATURE)
    axis.set_ylabel(r"Experimental $\Delta\Delta G^{\ddagger}$ (kcal mol$^{-1}$)")
    axis.legend(frameon=False)
    axis.set_title("Published Ni-hDA Physical-Organic Relationship")
    figure.tight_layout()
    figure.savefig(output_dir / "ni_hda_relationship.png", dpi=400)
    plt.close(figure)

    residuals = loo_predictions - y_train
    figure, axis = plt.subplots(figsize=(6.4, 5.2))
    axis.axhline(0.0, color="#333333", linestyle="--", linewidth=1.1)
    axis.scatter(loo_predictions, residuals, color="#176B87", s=60)
    axis.set_xlabel(r"LOO predicted $\Delta\Delta G^{\ddagger}$ (kcal mol$^{-1}$)")
    axis.set_ylabel("Residual (predicted - experimental)")
    axis.set_title("Fixed-Feature LOO Residuals")
    figure.tight_layout()
    figure.savefig(output_dir / "ni_hda_residuals.png", dpi=400)
    plt.close(figure)


def write_model_card(
    output_dir: Path,
    results: dict[str, object],
) -> None:
    """Create a concise publication-style Markdown study report."""
    published = results["published_model"]
    assert isinstance(published, dict)
    training = published["training_metrics"]
    loo = published["fixed_feature_loo_metrics"]
    blind = results["historical_blind_evaluation"]
    native = results.get("stericx_native_ensemble_model")
    assert isinstance(training, dict)
    assert isinstance(loo, dict)
    assert isinstance(blind, dict)
    blind_table_row = (
        f"| 723 | {blind['predicted_ddg_kcal_mol']:.4f} | "
        f"{blind['experimental_ddg_kcal_mol']:.4f} | {blind['mae']:.4f} | "
        f"{blind['applicability_domain']} |"
    )
    report = f"""# StericX Study 001

## Reproduction target

This study reproduces the enantioselectivity analysis in the Sigman Group's
[Ni-catalyzed homo-Diels-Alder repository]({PUBLISHED_NOTEBOOK_URL}). The
complete official Kraken table contains {results["source_rows"]:,} ligands, but
only 11 have experimental enantioselectivity labels. The published notebook
defines ten training ligands; source ligand 723 is reserved here as a
historical blind holdout.

## Preregistered model

The descriptor `{PUBLISHED_FEATURE}` was fixed from the published notebook
before model fitting. No StericX feature search used the blind target.

| Quantity | Value |
|---|---:|
| Training observations | {published["training_count"]} |
| Training R² | {training["r2"]:.4f} |
| Training RMSE | {training["rmse"]:.4f} kcal/mol |
| Fixed-feature LOO Q² | {loo["r2"]:.4f} |
| Fixed-feature LOO RMSE | {loo["rmse"]:.4f} kcal/mol |
| Response-permutation p-value | {published["response_permutation_p_value"]:.4f} |
| Slope | {published["slope"]:.6f} |
| Intercept | {published["intercept"]:.6f} |

![Historical held-out parity](ni_hda_parity.png)

## Frozen historical holdout

The prediction was written to `frozen_predictions.csv` before the target was
joined. Its SHA-256 digest is
`{results["frozen_prediction_sha256"]}`.

| Source ID | Predicted ΔΔG‡ | Experimental ΔΔG‡ | Absolute error | Domain |
|---:|---:|---:|---:|---|
{blind_table_row}

![Physical-organic relationship](ni_hda_relationship.png)

## Native StericX descriptor ablation

The native StericX model uses only the ETKDGv3/MMFF94 ensemble Sterimol
descriptors and reported donor NBO charge. The source has no reaction-specific
IR measurement, so its constant 1650 cm⁻¹ placeholder is automatically rejected
as non-informative. This deliberately tests whether the compact StericX
descriptor set can replace the published Kraken buried-volume feature.

{native_model_markdown(native)}

## Statistical controls

- Descriptor scaling for ridge and LASSO is learned within every inner
  cross-validation fold.
- Regularization is selected by nested leave-one-out cross-validation.
- Bootstrap intervals and Y-scrambling use deterministic recorded seeds.
- Correlation and VIF tables are saved as machine-readable CSV files.
- Applicability is reported against the training descriptor range.

![LOO residuals](ni_hda_residuals.png)

## Interpretation and limitations

The primary descriptor measures conformer-sensitive variation in buried-volume
anisotropy. Its positive coefficient is consistent with increasing asymmetric
steric differentiation accompanying larger ΔΔG‡ in this reaction family.

This is a faithful **historical reproduction**, not a new prospective
experiment. Ten training points cannot establish broad catalyst
generalizability, and one holdout cannot support a population-level R².
The holdout error must be reported as-is. A publication-grade prospective claim
requires predictions recorded before new reactions are performed, preferably
across an entire ligand scaffold or mechanistic regime.

## Provenance

- Source: {SOURCE_URL}
- Source SHA-256: `{results["source_sha256"]}`
- Generated: {results["generated_at_utc"]}
- No workstation specifications are recorded.
"""
    atomic_write_text(output_dir / "STUDY_001.md", report)


def native_model_markdown(native: object) -> str:
    """Render measured native-model results without inventing missing data."""
    if not isinstance(native, dict):
        return (
            "The native Rust model artifacts were not present when this report "
            "was generated."
        )
    training = native["training"]
    loo = native["fixed_feature_loo"]
    evaluation = native["historical_blind"]
    assert isinstance(training, dict)
    assert isinstance(loo, dict)
    assert isinstance(evaluation, dict)
    return f"""| Quantity | Value |
|---|---:|
| Selected term | `{native["selected_features"][0]}` |
| Training R² | {training["r2"]:.4f} |
| Fixed-feature LOO Q² | {loo["r2"]:.4f} |
| Scaffold-group LOO Q² | {native["fixed_feature_group_loo"]["r2"]:.4f} |
| Fixed-feature LOO RMSE | {loo["rmse"]:.4f} kcal/mol |
| Historical blind MAE | {evaluation["mae_kcal_mol"]:.4f} kcal/mol |

This ablation does **not** match the published Kraken model. The negative result
is retained because it identifies the missing scientific capability:
coordination-aware buried-volume descriptors are more consequential here than
additional inference optimization."""


def main(argv: Iterable[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        args.output_dir.mkdir(parents=True, exist_ok=True)
        frame, source_sha256 = load_catalog(args.catalog, args.offline)
        labeled_count = int(frame["ddG_abs"].notna().sum())
        print(
            f"Official Kraken catalog: {len(frame):,} ligands, "
            f"{labeled_count} enantioselectivity labels"
        )

        model, x_published, y_train, training_metrics = fit_published_model(frame)
        parametric_p_value = float(stats.pearsonr(x_published[:, 0], y_train).pvalue)
        loo_predictions, loo_metrics = fixed_feature_loo(x_published, y_train)
        x_curated = frame.loc[list(TRAIN_IDS), list(CURATED_FEATURES)].to_numpy(
            dtype=float
        )
        if not np.isfinite(x_curated).all():
            raise ValueError("curated descriptor matrix contains missing values")
        ridge = nested_regularized_loo(
            x_curated,
            y_train,
            CURATED_FEATURES,
            "ridge",
        )
        lasso = nested_regularized_loo(
            x_curated,
            y_train,
            CURATED_FEATURES,
            "lasso",
        )
        bootstrap = bootstrap_coefficients(
            x_published,
            y_train,
            args.bootstrap,
            args.seed,
        )
        permutation_p = response_permutation_test(
            x_published,
            y_train,
            float(training_metrics.r2 or 0.0),
            args.permutations,
            args.seed ^ 0x5A17,
        )
        correlations, vif = descriptor_diagnostics(x_curated, CURATED_FEATURES)
        atomic_write_csv(
            args.output_dir / "descriptor_correlations.csv",
            correlations.reset_index(names="feature"),
        )
        atomic_write_csv(args.output_dir / "descriptor_vif.csv", vif)

        model_artifact = {
            "schema_version": 1,
            "model": "published_univariate_ols",
            "feature": PUBLISHED_FEATURE,
            "training_source_ids": list(TRAIN_IDS),
            "blind_source_ids": list(BLIND_IDS),
            "intercept": float(model.intercept_),
            "slope": float(model.coef_[0]),
            "training_feature_minimum": float(x_published.min()),
            "training_feature_maximum": float(x_published.max()),
            "source_sha256": source_sha256,
            "created_at_utc": datetime.now(UTC).isoformat(),
            **bootstrap,
        }
        atomic_write_text(
            args.output_dir / "published_model.json",
            json.dumps(model_artifact, indent=2, sort_keys=True) + "\n",
        )

        # This write is intentionally completed before target revelation.
        frozen, frozen_sha256 = freeze_blind_prediction(
            frame.drop(columns=["ddG_abs"]),
            model,
            args.output_dir,
        )
        revealed, blind_metrics = reveal_and_score(frame, frozen)
        atomic_write_csv(args.output_dir / "scored_blind_predictions.csv", revealed)
        plot_results(
            args.output_dir,
            x_published,
            y_train,
            model,
            loo_predictions,
            revealed,
        )

        native_summary: dict[str, object] | None = None
        native_model_path = args.output_dir / "stericx_model.json"
        native_evaluation_path = args.output_dir / "stericx_evaluation.json"
        if native_model_path.is_file() and native_evaluation_path.is_file():
            native_model = json.loads(native_model_path.read_text(encoding="utf-8"))
            native_evaluation = json.loads(
                native_evaluation_path.read_text(encoding="utf-8")
            )
            native_summary = {
                "selected_features": native_model["selected_features"],
                "training": native_model["training"],
                "fixed_feature_loo": native_model["fixed_feature_loo"],
                "fixed_feature_group_loo": native_model["fixed_feature_group_loo"],
                "ridge_baseline": native_model["ridge_baseline"],
                "lasso_baseline": native_model["lasso_baseline"],
                "response_permutation_p_value": native_model[
                    "response_permutation_p_value"
                ],
                "historical_blind": {
                    "mae_kcal_mol": native_evaluation["mae_kcal_mol"],
                    "rmse_kcal_mol": native_evaluation["rmse_kcal_mol"],
                    "applicability_warnings": native_evaluation[
                        "applicability_warnings"
                    ],
                },
            }

        results: dict[str, object] = {
            "schema_version": 1,
            "generated_at_utc": datetime.now(UTC).isoformat(),
            "source_url": SOURCE_URL,
            "source_sha256": source_sha256,
            "source_rows": len(frame),
            "labeled_rows": labeled_count,
            "published_model": {
                "training_count": len(TRAIN_IDS),
                "feature": PUBLISHED_FEATURE,
                "slope": float(model.coef_[0]),
                "intercept": float(model.intercept_),
                "training_metrics": asdict(training_metrics),
                "fixed_feature_loo_metrics": asdict(loo_metrics),
                "response_permutation_p_value": permutation_p,
                "parametric_correlation_p_value": parametric_p_value,
                **bootstrap,
            },
            "ridge_nested_loo": asdict(ridge),
            "lasso_nested_loo": asdict(lasso),
            "frozen_prediction_sha256": frozen_sha256,
            "stericx_native_ensemble_model": native_summary,
            "historical_blind_evaluation": {
                "source_id": int(revealed.iloc[0]["Source_ID"]),
                "predicted_ddg_kcal_mol": float(
                    revealed.iloc[0]["Predicted_ddG_kcal_mol"]
                ),
                "experimental_ddg_kcal_mol": float(
                    revealed.iloc[0]["Experimental_ddG_kcal_mol"]
                ),
                "residual_kcal_mol": float(revealed.iloc[0]["Residual_kcal_mol"]),
                "applicability_domain": str(revealed.iloc[0]["Applicability_Domain"]),
                **asdict(blind_metrics),
            },
            "limitations": [
                "Only ten labeled observations are used for model fitting.",
                "The single holdout is historical, not newly prospective.",
                "One holdout supports an error estimate but not a holdout R2.",
                "The complete 1,566-ligand library must not be treated as labeled.",
            ],
        }
        atomic_write_text(
            args.output_dir / "study_results.json",
            json.dumps(results, indent=2, sort_keys=True) + "\n",
        )
        write_model_card(args.output_dir, results)

        print("\nStericX Study 001 complete")
        print(
            f"  Published model: R²={training_metrics.r2:.4f}, "
            f"fixed-feature LOO Q²={loo_metrics.r2:.4f}"
        )
        print(
            f"  Nested ridge LOO: R²={ridge.metrics.r2:.4f}, "
            f"RMSE={ridge.metrics.rmse:.4f} kcal/mol"
        )
        print(
            f"  Nested LASSO LOO: R²={lasso.metrics.r2:.4f}, "
            f"RMSE={lasso.metrics.rmse:.4f} kcal/mol"
        )
        print(
            f"  Historical blind 723: predicted="
            f"{revealed.iloc[0]['Predicted_ddG_kcal_mol']:.4f}, "
            f"experimental={revealed.iloc[0]['Experimental_ddG_kcal_mol']:.4f}, "
            f"MAE={blind_metrics.mae:.4f} kcal/mol"
        )
        print(f"  Report: {args.output_dir / 'STUDY_001.md'}")
        return 0
    except (
        FileNotFoundError,
        KeyError,
        OSError,
        requests.RequestException,
        RuntimeError,
        ValueError,
    ) as exc:
        print(f"Study failed: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
