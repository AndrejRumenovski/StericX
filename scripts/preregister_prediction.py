#!/usr/bin/env python3
"""Turn the frozen Ni-hDA ligand deck into a pre-registered, falsifiable forecast.

Study 003 froze a target-free deck of ten ligands with point predictions of the
enantiodetermining free-energy difference (ddG-double-dagger) from the published
univariate model (Study 001). Point predictions are not a pre-registration: a
domain expert needs, before any measurement, (1) a quantitative uncertainty on
every prediction, (2) a statement of whether each ligand lies inside the model's
applicability domain, and (3) an exact experimental protocol with a decision rule
that will corroborate or falsify the model.

This script adds those three things WITHOUT touching the frozen deck. It reads the
deck, the published model, and the training catalog; recomputes the ordinary
least-squares prediction intervals and leverages from the ten training points;
and writes a cryptographically anchored pre-registration (`PREREGISTRATION.md`,
`preregistration.json`, an intervals table, and a forecast figure) that commits to
the frozen deck by its SHA-256. The point predictions are unchanged and are
verified byte-for-byte against the frozen deck; only the honest envelope around
them is new.

Nothing here is a synthesis or safety instruction. The candidates are
DFT-characterized Kraken ligands chosen computationally to span and challenge the
model's domain; synthetic accessibility, stability, and safety require expert
review before any experiment.
"""

from __future__ import annotations

import argparse
import hashlib
import io
import json
from datetime import UTC, datetime
from pathlib import Path
from typing import Final

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

ROOT: Final[Path] = Path(__file__).resolve().parent.parent

# Ni-hDA training ligands (Study 001 published model). Kept in sync with the
# freeze script and the model manifest; asserted against the model at run time.
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
FEATURE: Final[str] = "vbur_max_delta_qvbur_min"
GAS_CONSTANT_KCAL_MOL_K: Final[float] = 0.00198720425864083
TEMPERATURE_K: Final[float] = 298.15
CONFIDENCE: Final[float] = 0.95

# Pre-registered falsification rule (chosen before any measurement): the primary
# test is on the interpolation candidates; corroboration requires prediction-
# interval coverage AND correct ordering. Boundary ligands are scored separately.
PRIMARY_STRATUM: Final[str] = "interpolation_maximin"
COVERAGE_MIN_FRACTION: Final[float] = 0.75  # >= 6 of 8 within the 95% interval

# The reaction the model speaks to, cited (not reproduced) from its public source.
REACTION: Final[str] = "Ni-catalyzed homo-Diels-Alder (hDA), asymmetric variant"
REACTION_SOURCE: Final[str] = (
    "Sigman Group, Ni-Catalyzed-hDA repository "
    "(https://github.com/SigmanGroup/Ni-Catalyzed-hDA); the univariate "
    "enantioselectivity model reproduced in StericX Study 001."
)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--deck",
        type=Path,
        default=ROOT / "docs" / "study_003" / "prospective_ligand_deck.csv",
    )
    parser.add_argument(
        "--deck-manifest",
        type=Path,
        default=ROOT / "docs" / "study_003" / "prospective_deck_manifest.json",
    )
    parser.add_argument(
        "--model",
        type=Path,
        default=ROOT / "docs" / "study_001" / "published_model.json",
    )
    parser.add_argument(
        "--catalog",
        type=Path,
        default=ROOT / "data" / "official" / "ni_hda_kraken.csv",
    )
    parser.add_argument("--output-dir", type=Path, default=ROOT / "docs" / "study_003")
    return parser.parse_args(argv)


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _repo_relative(path: Path) -> str:
    """Path relative to the repo root, so no absolute home paths are committed."""
    try:
        return str(path.resolve().relative_to(ROOT))
    except ValueError:
        return str(path)


def ee_percent(ddg: np.ndarray | float) -> np.ndarray | float:
    """Enantiomeric excess (%) from ddG-double-dagger via the Eyring ratio."""
    ratio = np.exp(np.asarray(ddg) / (GAS_CONSTANT_KCAL_MOL_K * TEMPERATURE_K))
    return 100.0 * (ratio - 1.0) / (ratio + 1.0)


class OlsForecast:
    """Univariate OLS prediction intervals and leverages from the training fit."""

    def __init__(self, x: np.ndarray, y: np.ndarray, slope: float, intercept: float):
        self.slope = slope
        self.intercept = intercept
        self.n = len(x)
        self.x_mean = float(x.mean())
        self.sxx = float(np.sum((x - self.x_mean) ** 2))
        residuals = y - (intercept + slope * x)
        # Residual standard error of the published fit (n - 2 parameters).
        self.residual_se = float(np.sqrt(np.sum(residuals**2) / (self.n - 2)))
        self.t_value = float(stats.t.ppf(0.5 + CONFIDENCE / 2.0, self.n - 2))
        # Warning leverage h* = 3p/n (p = 2): the standard OLS applicability bound.
        self.warning_leverage = 3.0 * 2.0 / self.n

    def leverage(self, x0: np.ndarray) -> np.ndarray:
        return 1.0 / self.n + (x0 - self.x_mean) ** 2 / self.sxx

    def prediction(self, x0: np.ndarray) -> np.ndarray:
        return self.intercept + self.slope * x0

    def interval_halfwidth(self, x0: np.ndarray) -> np.ndarray:
        return self.t_value * self.residual_se * np.sqrt(1.0 + self.leverage(x0))


def build_intervals(deck: pd.DataFrame, forecast: OlsForecast) -> pd.DataFrame:
    x0 = deck[FEATURE].to_numpy(dtype=float)
    predicted = forecast.prediction(x0)
    half = forecast.interval_halfwidth(x0)
    low, high = predicted - half, predicted + half
    leverage = forecast.leverage(x0)
    table = pd.DataFrame(
        {
            "Source_ID": deck["Source_ID"].to_numpy(int),
            "Selection_Stratum": deck["Selection_Stratum"],
            "vbur_max_delta_qvbur_min": x0,
            "Predicted_ddG_kcal_mol": predicted,
            "ddG_95pi_low": low,
            "ddG_95pi_high": high,
            "Predicted_ee_percent": ee_percent(predicted),
            "ee_95pi_low": ee_percent(low),
            "ee_95pi_high": ee_percent(high),
            "leverage": leverage,
            "in_leverage_domain": leverage <= forecast.warning_leverage,
            # A prediction is only directionally falsifiable on ee when its
            # interval excludes the racemate (both bounds share ddG sign).
            "ee_interval_excludes_racemate": np.sign(low) == np.sign(high),
        }
    )
    return table


def deterministic_csv(frame: pd.DataFrame) -> bytes:
    buffer = io.StringIO()
    frame.to_csv(buffer, index=False, float_format="%.10g", lineterminator="\n")
    return buffer.getvalue().encode()


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)

    deck_bytes = args.deck.read_bytes()
    deck_sha = sha256_bytes(deck_bytes)
    deck = pd.read_csv(io.BytesIO(deck_bytes))

    manifest = json.loads(args.deck_manifest.read_text(encoding="utf-8"))
    if manifest.get("deck_sha256") != deck_sha:
        raise SystemExit(
            "deck SHA-256 does not match its manifest; the frozen deck has changed."
        )

    model_bytes = args.model.read_bytes()
    model = json.loads(model_bytes)
    if model.get("feature") != FEATURE:
        raise SystemExit("model does not use the pre-registered buried-volume feature")
    if tuple(model.get("training_source_ids", ())) != TRAIN_IDS:
        raise SystemExit("model training_source_ids do not match TRAIN_IDS")
    slope, intercept = float(model["slope"]), float(model["intercept"])

    catalog_bytes = args.catalog.read_bytes()
    catalog = pd.read_csv(io.BytesIO(catalog_bytes), index_col=0)
    catalog.index = pd.to_numeric(catalog.index).astype(int)
    train = catalog.loc[list(TRAIN_IDS)]
    x_train = train[FEATURE].to_numpy(dtype=float)
    y_train = train["ddG_abs"].to_numpy(dtype=float)
    if not np.isfinite(y_train).all():
        raise SystemExit("training ligands must all carry a measured ddG_abs")

    forecast = OlsForecast(x_train, y_train, slope, intercept)
    intervals = build_intervals(deck, forecast)

    intervals_bytes = deterministic_csv(intervals)
    intervals_sha = sha256_bytes(intervals_bytes)
    intervals_path = args.output_dir / "prospective_prediction_intervals.csv"
    intervals_path.write_bytes(intervals_bytes)

    primary = intervals[intervals["Selection_Stratum"] == PRIMARY_STRATUM]
    secondary = intervals[intervals["Selection_Stratum"] != PRIMARY_STRATUM]
    n_primary = len(primary)
    coverage_needed = int(np.ceil(COVERAGE_MIN_FRACTION * n_primary))

    prereg_at = datetime.now(UTC).isoformat()
    figure_path = args.output_dir / "prospective_prediction_forecast.png"
    write_forecast_figure(intervals, figure_path)

    preregistration = {
        "schema_version": 1,
        "preregistered_at_utc": prereg_at,
        "status": "predictions_frozen_measurements_pending",
        "reaction": REACTION,
        "reaction_source": REACTION_SOURCE,
        "prediction_target": "ddG_double_dagger_kcal_mol (enantiodetermining step)",
        "model": {
            "kind": "published univariate OLS (Study 001)",
            "feature": FEATURE,
            "slope": slope,
            "intercept": intercept,
            "training_n": forecast.n,
            "residual_standard_error_kcal_mol": forecast.residual_se,
            "training_feature_range": [float(x_train.min()), float(x_train.max())],
        },
        "uncertainty_method": {
            "kind": f"OLS {int(CONFIDENCE * 100)}% prediction interval",
            "formula": "yhat +/- t(0.975, n-2) * s * sqrt(1 + 1/n + (x0-xbar)^2/Sxx)",
            "t_value": forecast.t_value,
        },
        "applicability_domain": {
            "kind": "OLS leverage vs warning leverage h* = 3p/n",
            "warning_leverage": forecast.warning_leverage,
            "in_domain_count": int(intervals["in_leverage_domain"].sum()),
            "note": (
                "Leverage and the deck's 1-D training-range label can disagree "
                "for a ligand just outside the feature min/max; both are reported."
            ),
        },
        "falsification_rule": {
            "primary_set": f"{n_primary} {PRIMARY_STRATUM} candidates",
            "corroborated_if": (
                f"measured ddG within the 95% prediction interval for at least "
                f"{coverage_needed} of {n_primary} primary candidates AND "
                f"Spearman rho(predicted, measured) > 0"
            ),
            "falsified_otherwise": True,
            "coverage_needed": coverage_needed,
            "coverage_min_fraction": COVERAGE_MIN_FRACTION,
            "secondary_set": (
                "boundary-challenge ligands, scored as an extrapolation "
                "stress-test; reported but not gating the primary test"
            ),
            "honest_power_caveat": (
                "Near-racemate predictions have ee intervals spanning both signs "
                "and are not meaningfully falsifiable on ee; the ddG interval "
                "remains the test quantity."
            ),
        },
        "experimental_protocol": {
            "reaction": REACTION,
            "source": REACTION_SOURCE,
            "measurement": (
                "Run the asymmetric Ni-hDA reaction with each candidate ligand "
                "under the source conditions; determine product enantiomeric "
                "excess by chiral HPLC or SFC."
            ),
            "ee_to_ddG": (
                "Convert measured ee (fraction) to ddG = 2 R T atanh(ee) at "
                f"T = {TEMPERATURE_K} K, then compare to the interval above."
            ),
            "not_an_instruction": (
                "Candidates are computational recommendations; synthetic "
                "accessibility, stability, and safety require expert review."
            ),
        },
        "commits_to": {
            "deck": _repo_relative(args.deck),
            "deck_sha256": deck_sha,
            "deck_frozen_at_utc": manifest.get("frozen_at_utc"),
            "model_sha256": sha256_bytes(model_bytes),
            "catalog_sha256": sha256_bytes(catalog_bytes),
            "intervals_csv_sha256": intervals_sha,
        },
        "supersedes": (
            "adds an uncertainty/applicability-domain/protocol layer over the "
            "frozen deck; the point predictions are unchanged and verified "
            "byte-for-byte against it"
        ),
        "limitations": [
            "Ten training points; a univariate steric model of one reaction family.",
            "Candidates are DFT-characterized, not necessarily synthesizable.",
            "Must not be refit after any outcome is known.",
        ],
    }
    prereg_json_path = args.output_dir / "preregistration.json"
    prereg_json_path.write_text(
        json.dumps(preregistration, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    write_preregistration_md(
        args.output_dir / "PREREGISTRATION.md",
        preregistration,
        intervals,
        primary,
        secondary,
    )

    print("Pre-registration written")
    print(
        f"  primary set: {n_primary} candidates; corroboration needs "
        f">= {coverage_needed} within 95% PI + positive rank"
    )
    print(
        f"  in leverage domain: {int(intervals['in_leverage_domain'].sum())}/"
        f"{len(intervals)}"
    )
    print(f"  deck_sha256 committed: {deck_sha}")
    print(f"  intervals_csv_sha256: {intervals_sha}")
    return 0


def _interval_panel(
    axis, positions, labels, values, low, high, colors, ylabel, title, ylim=None
):
    """Draw one per-candidate 95%-interval error-bar panel."""
    for pos, value, lo, hi, color in zip(
        positions, values, values - low, high - values, colors, strict=True
    ):
        axis.errorbar(
            pos,
            value,
            yerr=[[lo], [hi]],
            fmt="o",
            color=color,
            capsize=3,
            markersize=6,
            linewidth=1.3,
        )
    axis.axhline(0.0, color="#999999", linewidth=0.8, linestyle=":")
    if ylim is not None:
        axis.set_ylim(*ylim)
    axis.set_xticks(positions)
    axis.set_xticklabels(labels, rotation=45, ha="right", fontsize=8)
    axis.set_xlabel("candidate (Kraken Source_ID)")
    axis.set_ylabel(ylabel)
    axis.set_title(title)


def write_forecast_figure(intervals: pd.DataFrame, output: Path) -> None:
    ordered = intervals.sort_values("Predicted_ddG_kcal_mol").reset_index(drop=True)
    labels = [str(i) for i in ordered["Source_ID"]]
    colors = [
        "#176B87" if flag else "#B4530A" for flag in ordered["in_leverage_domain"]
    ]
    positions = np.arange(len(ordered))

    figure, (ddg_ax, ee_ax) = plt.subplots(1, 2, figsize=(12, 5.2))
    _interval_panel(
        ddg_ax,
        positions,
        labels,
        ordered["Predicted_ddG_kcal_mol"].to_numpy(),
        ordered["ddG_95pi_low"].to_numpy(),
        ordered["ddG_95pi_high"].to_numpy(),
        colors,
        "predicted ddG-double-dagger (kcal/mol) +/- 95% PI",
        "Frozen forecast: enantioselectivity free energy",
    )
    _interval_panel(
        ee_ax,
        positions,
        labels,
        ordered["Predicted_ee_percent"].to_numpy(),
        ordered["ee_95pi_low"].to_numpy(),
        ordered["ee_95pi_high"].to_numpy(),
        colors,
        "predicted ee (%) +/- 95% PI",
        "Same forecast as ee: intervals cross racemic at low ee",
        ylim=(-100, 100),
    )

    handles = [
        plt.Line2D(
            [],
            [],
            marker="o",
            linestyle="",
            color="#176B87",
            label="inside leverage domain",
        ),
        plt.Line2D(
            [],
            [],
            marker="o",
            linestyle="",
            color="#B4530A",
            label="extrapolation (leverage > h*)",
        ),
    ]
    ddg_ax.legend(handles=handles, loc="upper left", fontsize=8, frameon=False)
    figure.suptitle(
        "Study 003 pre-registered Ni-hDA forecast (frozen; measurements pending)",
        fontsize=11,
    )
    figure.tight_layout()
    figure.savefig(output, dpi=300)
    plt.close(figure)


def _pi_table(rows: pd.DataFrame) -> list[str]:
    lines = [
        "| Source_ID | %Vbur feature | pred ddG | 95% PI | pred ee | ee 95% PI | "
        "leverage | in domain |",
        "|---:|---:|---:|---|---:|---|---:|:--:|",
    ]
    for _, r in rows.iterrows():
        lines.append(
            f"| {int(r['Source_ID'])} | {r['vbur_max_delta_qvbur_min']:.2f} | "
            f"{r['Predicted_ddG_kcal_mol']:.2f} | "
            f"[{r['ddG_95pi_low']:.2f}, {r['ddG_95pi_high']:.2f}] | "
            f"{r['Predicted_ee_percent']:.0f}% | "
            f"[{r['ee_95pi_low']:.0f}, {r['ee_95pi_high']:.0f}]% | "
            f"{r['leverage']:.2f} | {'yes' if r['in_leverage_domain'] else 'NO'} |"
        )
    return lines


def write_preregistration_md(
    output: Path,
    prereg: dict,
    intervals: pd.DataFrame,
    primary: pd.DataFrame,
    secondary: pd.DataFrame,
) -> None:
    model = prereg["model"]
    fr = prereg["falsification_rule"]
    commits = prereg["commits_to"]
    racemate_ambiguous = int((~intervals["ee_interval_excludes_racemate"]).sum())
    lines = [
        "# StericX Study 003 - Pre-registered Prospective Prediction",
        "",
        "**A frozen, falsifiable forecast for the "
        f"{prereg['reaction']}.** This document commits, before any measurement, "
        "to point predictions of the enantiodetermining free-energy difference "
        "(ddG-double-dagger) for ten ligands, each with a 95% prediction interval "
        "and an applicability-domain judgement, plus the exact experimental test "
        "that would corroborate or falsify the model. It is a pre-registration, "
        "not a synthesis or safety instruction.",
        "",
        f"- **Pre-registered (UTC):** `{prereg['preregistered_at_utc']}`",
        f"- **Frozen deck committed (SHA-256):** `{commits['deck_sha256']}`",
        f"  (deck frozen at `{commits['deck_frozen_at_utc']}`)",
        f"- **Model (SHA-256):** `{commits['model_sha256']}` - "
        f"{model['kind']}, feature `{model['feature']}`",
        f"- **Training catalog (SHA-256):** `{commits['catalog_sha256']}`",
        f"- **Intervals table (SHA-256):** `{commits['intervals_csv_sha256']}`",
        "",
        "The point predictions are identical to the frozen deck (verified "
        "byte-for-byte by the committed SHA-256); this document only adds the "
        "uncertainty, the applicability domain, and the falsification protocol.",
        "",
        "## The model and its uncertainty",
        "",
        f"The predictor is the {model['kind']}: "
        f"ddG = {model['intercept']:.3f} + {model['slope']:.3f} x "
        f"`{model['feature']}`, fit on **{model['training_n']} ligands** with a "
        f"residual standard error of **{model['residual_standard_error_kcal_mol']:.3f} "
        "kcal/mol**. Each prediction carries an ordinary-least-squares "
        "95% prediction "
        "interval, `yhat +/- t(0.975, n-2) * s * sqrt(1 + 1/n + (x0-xbar)^2/Sxx)`, "
        "which widens for ligands far from the training mean.",
        "",
        "## Applicability domain",
        "",
        f"Domain membership uses the OLS leverage against the warning leverage "
        f"**h\\* = 3p/n = {prereg['applicability_domain']['warning_leverage']:.2f}**. "
        f"**{prereg['applicability_domain']['in_domain_count']} of "
        f"{len(intervals)}** candidates fall inside it. "
        + prereg["applicability_domain"]["note"],
        "",
        "## The frozen forecast",
        "",
        "**Primary set (interpolation candidates - the pre-registered test):**",
        "",
        *_pi_table(primary),
        "",
        "**Secondary set (boundary challenges - extrapolation stress-test):**",
        "",
        *_pi_table(secondary),
        "",
        "![Pre-registered forecast](prospective_prediction_forecast.png)",
        "",
        "*Figure. Predicted ddG-double-dagger and ee for every candidate with 95% "
        "prediction intervals; blue = inside the leverage domain, orange = "
        "extrapolation. Generated by `scripts/preregister_prediction.py`.*",
        "",
        "### The honest limit of these predictions",
        "",
        f"For **{racemate_ambiguous} of {len(intervals)}** candidates the 95% ee "
        "interval spans both signs of enantioselectivity - i.e. the model cannot "
        "even commit to which enantiomer is favoured. Near the racemate a "
        f"{model['residual_standard_error_kcal_mol']:.2f} kcal/mol residual maps to "
        "roughly +/-60% ee, so those predictions are **not meaningfully falsifiable "
        "on ee**. The genuinely discriminating, falsifiable predictions are the "
        "high-buried-volume ligands whose intervals stay well above the racemate. "
        "The test below is therefore run on ddG, where the interval is informative "
        "for every ligand.",
        "",
        "## Pre-registered falsification rule",
        "",
        f"**Primary test** ({fr['primary_set']}):",
        "",
        f"> The model is **corroborated** iff the measured ddG-double-dagger lies "
        f"within the 95% prediction interval for at least **{fr['coverage_needed']}** "
        f"of the {len(primary)} primary candidates **and** the Spearman rank "
        "correlation between predicted and measured ddG is **positive**. "
        "Otherwise it is **falsified**.",
        "",
        f"**Secondary test:** {fr['secondary_set']}.",
        "",
        "This rule is fixed now and must not be adjusted after outcomes are known. "
        "A 95% interval is expected to miss about one candidate in twenty by "
        "chance, so the coverage bar tolerates noise while a wrong sign or more "
        "than two misses falsifies; the rank clause guards against a flat fit that "
        "happens to sit inside wide intervals.",
        "",
        "## Experimental protocol",
        "",
        f"- **Reaction:** {prereg['experimental_protocol']['reaction']}.",
        f"- **Source:** {prereg['experimental_protocol']['source']}",
        f"- **Measurement:** {prereg['experimental_protocol']['measurement']}",
        f"- **Scoring:** {prereg['experimental_protocol']['ee_to_ddG']}",
        "",
        "## Limitations (kept in view)",
        "",
        *[f"- {item}" for item in prereg["limitations"]],
        "",
        "This is the natural falsifiable hook for outreach: a frozen, dated, "
        "hash-anchored prediction with honest uncertainty, and the exact "
        "experiment that would prove it right or wrong.",
        "",
    ]
    output.write_text("\n".join(lines), encoding="utf-8")


if __name__ == "__main__":
    raise SystemExit(main())
