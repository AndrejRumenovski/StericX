#!/usr/bin/env python3
"""Study 011: retrospective, outcome-blinded ligand ranking.

The driver has an explicit design stage and an analysis stage. ``--prepare-only``
writes and hashes the outcome-independent design, exact scaffold-disjoint splits,
and seeds. A full run refuses to proceed if those files no longer match their
design lock. All StericX fits and screens finish, and their target-free predictions
are written and hashed, before ``ddG_abs`` is read for evaluation.

Reproduce the checked-in study with:

    uv run --extra science python studies/study_011_ligand_ranking.py

Use ``--prepare-only`` to regenerate only the preregistration artifacts and
``--verify-only`` to audit an existing result bundle without refitting.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import itertools
import json
import math
import os
import platform
import subprocess
import tempfile
from collections import Counter, defaultdict
from collections.abc import Iterable, Sequence
from pathlib import Path
from typing import Any, Final

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import scipy
from scipy.stats import rankdata

ROOT: Final[Path] = Path(__file__).resolve().parent.parent
STUDY_DIR: Final[Path] = ROOT / "docs" / "study_011"
METADATA: Final[Path] = ROOT / "data" / "reactions_raw.csv"
CATALOG: Final[Path] = ROOT / "data" / "official" / "ni_hda_kraken.csv"
CATALOG_PROVENANCE: Final[Path] = ROOT / "data" / "official" / "provenance.json"
SIGPACK: Final[Path] = ROOT / "data" / "reactions.sigpack"
BINARY: Final[Path] = ROOT / "target" / "release" / "stericx"
SOURCE_URL: Final[str] = (
    "https://raw.githubusercontent.com/SigmanGroup/"
    "Ni-Catalyzed-hDA/main/data/kraken.csv"
)
BASE_SEED: Final[int] = 20_260_903
BOOTSTRAP_SAMPLES: Final[int] = 1_000
PERMUTATION_SAMPLES: Final[int] = 1_000
MAX_TERMS: Final[int] = 2
CANDIDATE_COUNT: Final[int] = 3
DOMAIN_RULE: Final[str] = "max-neighbor"
DESIGN_FILES: Final[tuple[str, ...]] = (
    "DESIGN.md",
    "design.json",
    "split_definitions.json",
    "splits.csv",
    "seeds.json",
)
RESULT_FILES: Final[tuple[str, ...]] = (
    "frozen_predictions.csv",
    "frozen_predictions.sha256",
    "rankings.csv",
    "split_metrics.csv",
    "random_baseline_metrics.csv",
    "model_reports.jsonl",
    "screen_reports.jsonl",
    "split_status.csv",
    "metrics.json",
    "run_manifest.json",
    "ranking_performance.png",
    "applicability_domain.png",
    "STUDY_011.md",
)


def parse_args(argv: Iterable[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument(
        "--prepare-only",
        action="store_true",
        help="write and lock the outcome-independent design, then stop",
    )
    mode.add_argument(
        "--verify-only",
        action="store_true",
        help="verify the design lock, result hashes, and metric reconstruction",
    )
    parser.add_argument(
        "--binary",
        type=Path,
        default=BINARY,
        help="StericX release binary",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=STUDY_DIR,
        help="study artifact directory",
    )
    return parser.parse_args(argv)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def atomic_write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(text, encoding="utf-8")
    os.replace(temporary, path)


def atomic_write_json(path: Path, value: Any) -> None:
    atomic_write_text(path, json.dumps(value, indent=2, sort_keys=True) + "\n")


def atomic_write_csv(
    path: Path,
    fieldnames: Sequence[str],
    rows: Iterable[dict[str, Any]],
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    with temporary.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {
                    key: "" if row.get(key) is None else row.get(key, "")
                    for key in fieldnames
                }
            )
    os.replace(temporary, path)


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def read_allowed_metadata() -> list[dict[str, str]]:
    """Read predictors and identifiers, deliberately never indexing the target."""
    allowed = (
        "Reaction_ID",
        "Ligand_XYZ_Path",
        "NBO_Charge",
        "IR_Frequency",
        "Ligand_SMILES",
        "Ligand_Group",
        "Source_ID",
        "Source_URL",
    )
    with METADATA.open(newline="", encoding="utf-8") as handle:
        reader = csv.reader(handle)
        header = next(reader)
        indices = {name: header.index(name) for name in allowed}
        rows = [
            {name: record[index] for name, index in indices.items()}
            for record in reader
        ]
    if len(rows) != 11:
        raise ValueError(f"expected 11 labeled Ni-hDA rows, found {len(rows)}")
    if len({row["Reaction_ID"] for row in rows}) != len(rows):
        raise ValueError("Reaction_ID values are not unique")
    if any(not row["Ligand_Group"].strip() for row in rows):
        raise ValueError("every row must carry a non-empty Ligand_Group")
    return rows


def numeric_id(reaction_id: str) -> int:
    return int(reaction_id.rsplit("-", 1)[-1])


def make_splits(rows: Sequence[dict[str, str]]) -> list[dict[str, Any]]:
    """Enumerate every whole-scaffold holdout containing exactly three rows."""
    group_to_ids: dict[str, list[str]] = defaultdict(list)
    for row in rows:
        group_to_ids[row["Ligand_Group"]].append(row["Reaction_ID"])
    groups = sorted(
        group_to_ids,
        key=lambda group: min(numeric_id(value) for value in group_to_ids[group]),
    )
    combinations = [
        combination
        for size in range(1, len(groups) + 1)
        for combination in itertools.combinations(groups, size)
        if sum(len(group_to_ids[group]) for group in combination) == CANDIDATE_COUNT
    ]
    combinations.sort(
        key=lambda combination: tuple(
            sorted(
                numeric_id(value)
                for group in combination
                for value in group_to_ids[group]
            )
        )
    )
    all_ids = [row["Reaction_ID"] for row in rows]
    splits = []
    for index, held_groups in enumerate(combinations, start=1):
        held_ids = sorted(
            (value for group in held_groups for value in group_to_ids[group]),
            key=numeric_id,
        )
        held_set = set(held_ids)
        train_ids = [value for value in all_ids if value not in held_set]
        splits.append(
            {
                "split_id": f"S{index:03d}",
                "stratum": (
                    "one_three_member_scaffold"
                    if len(held_groups) == 1
                    else "three_singleton_scaffolds"
                ),
                "training_ids": train_ids,
                "held_out_ids": held_ids,
                "held_out_groups": list(held_groups),
                "training_count": len(train_ids),
                "candidate_count": len(held_ids),
                "model_seed": BASE_SEED + index,
            }
        )
    if len(splits) != 57:
        raise ValueError(f"expected 57 exhaustive splits, generated {len(splits)}")
    for split in splits:
        if split["training_count"] != 8 or split["candidate_count"] != 3:
            raise ValueError(f"invalid split size in {split['split_id']}")
        train_groups = {
            row["Ligand_Group"]
            for row in rows
            if row["Reaction_ID"] in split["training_ids"]
        }
        if train_groups.intersection(split["held_out_groups"]):
            raise ValueError(f"scaffold leakage in {split['split_id']}")
    return splits


def design_document() -> str:
    return """# Study 011 design lock: retrospective ligand ranking

Status: **outcome-independent design, frozen before Study 011 evaluation**.

## Question and estimand

Can a StericX model fitted only to observed ligand outcomes rank a small panel of
previously unseen ligands so that the more enantioselective Ni-hDA ligands are
prioritized? The response is the published absolute transition-state free-energy
difference, `ddG_abs` (kcal/mol); larger is defined as better. The primary estimand
is top-1 recovery across all planned screens, including failed screens.

## Dataset and eligibility

Use the 11 experimentally labeled Ni-catalyzed homo-Diels-Alder ligands already
cached in `data/official/ni_hda_kraken.csv` and aligned to StericX records in
`data/reactions.sigpack`. The authors publicly distribute the source table at the
URL recorded in `data/official/provenance.json`. Unlabeled Kraken rows are ineligible
because they have no outcome against which a ranking can be evaluated.

## Splits fixed without outcomes

`Ligand_Group` is the grouping variable. Enumerate every combination of complete
groups whose combined size is exactly three ligands. This produces 57 exhaustive,
non-random splits: one holdout containing the sole three-member scaffold and 56
holdouts containing three singleton scaffolds. Each split has eight training
observations and a three-ligand candidate panel. No group may cross train/test.
Exact membership is in `split_definitions.json` and row-expanded `splits.csv`.

The exhaustive construction prevents choosing a favorable split. It does not make
the overlapping splits statistically independent. The three-member scaffold enters
one split, whereas each singleton enters 21 splits, so the two strata must also be
reported separately.

## Model and screening protocol

For each split, in this order:

1. Write row-aligned split metadata. Only the eight `train` records may contribute
   outcomes to fitting.
2. Run `stericx fit` with the native physical-organic feature space, BIC forward
   selection, at most two non-intercept terms, the built-in `|r| <= 0.95` guard,
   1,000 bootstrap samples, 1,000 response permutations, and the split seed in
   `seeds.json`. Scaling and selection occur on training rows only.
3. Create a three-row, target-free candidate CSV using an explicit predictor
   allow-list. Run `stericx screen` with `--optimize maximize`, bond-axis Sterimol,
   and the training-derived `max-neighbor` applicability rule. Do not filter
   extrapolations. Rank descending; StericX's identifier tie-break is retained.
4. Finish all 57 fits/screens, write `frozen_predictions.csv`, and hash it.
5. Only then read `ddG_abs`, join outcomes for scoring, and write revealed results.

A split that cannot fit or rank all three candidates remains a failure. It is not
dropped from primary intention-to-screen metrics.

## Metrics fixed in advance

Primary (all 57 planned splits):

- top-1 recovery: whether predicted rank 1 is the observed best ligand; a failed
  screen scores 0;
- top-2 recovery: recall of the observed top two among the predicted top two; a
  failed screen scores 0;
- rank of the observed best ligand (1--3); a failed screen scores 4.

Secondary (successfully ranked splits unless stated otherwise):

- Spearman correlation between predicted and observed response when neither vector
  is constant;
- nDCG over all three non-negative outcomes;
- MAE and RMSE in kcal/mol, both per-split and over all repeated predictions;
- observed outcome of the predicted top 1/top 2 and its gain and ratio relative to
  the candidate-panel mean;
- 95% prediction-interval empirical coverage;
- absolute error by StericX domain verdict and range status, plus descriptive
  association with nearest-training-distance ratio and leverage ratio.

Observed ranks use descending outcome then numeric ligand ID to resolve exact ties.
Pooled errors are descriptive because the same ligand appears in multiple splits.
No hypothesis-test threshold is used to select or suppress a result.

## Random baseline

Enumerate all `3! = 6` possible candidate rankings in every panel. This exact
baseline needs no random seed. Chance expectations are top-1 recovery `1/3`, top-2
recall `2/3`, mean best-ligand rank `2`, and mean Spearman `0`. Enrichment is the
StericX mean divided by the corresponding exact random mean. Outcome-selection
gain is compared with the exact panel-wise random expectation.

## Interpretation constraints

This is a small retrospective test of ranking three measured ligands, not a claim
about recovery from the full 1,566-ligand Kraken library, yield, synthetic
accessibility, prospective performance, or causal ligand effects. Failed and poor
results must remain visible. Because splits overlap, no naive confidence interval
or p-value treating 57 folds as independent will be reported.
"""


def design_spec(splits: Sequence[dict[str, Any]]) -> dict[str, Any]:
    strata = Counter(split["stratum"] for split in splits)
    return {
        "study": "study_011",
        "title": "Retrospective ligand ranking",
        "design_status": "fixed_before_study_011_outcome_revelation",
        "dataset": {
            "name": "Sigman Group Ni-catalyzed homo-Diels-Alder",
            "catalog": str(CATALOG.relative_to(ROOT)),
            "metadata": str(METADATA.relative_to(ROOT)),
            "sigpack": str(SIGPACK.relative_to(ROOT)),
            "source_url": SOURCE_URL,
            "eligible_labeled_ligands": 11,
            "response": "ddG_abs",
            "response_units": "kcal/mol",
            "optimization": "maximize",
        },
        "splitting": {
            "group_column": "Ligand_Group",
            "scheme": "exhaustive whole-group combinations totaling 3 observations",
            "planned_splits": len(splits),
            "training_observations_per_split": 8,
            "candidate_observations_per_split": 3,
            "strata": dict(sorted(strata.items())),
            "randomized": False,
        },
        "model": {
            "engine": "stericx fit",
            "family": "mechanistically_constrained_ols",
            "max_non_intercept_terms": MAX_TERMS,
            "selection": "training-only BIC forward selection",
            "correlation_guard_absolute": 0.95,
            "bootstrap_samples": BOOTSTRAP_SAMPLES,
            "permutation_samples": PERMUTATION_SAMPLES,
        },
        "screen": {
            "engine": "stericx screen",
            "candidate_set": "three held-out ligands only",
            "ranking": "predicted ddG descending; ligand identifier tie-break",
            "sterimol_axis": "bond",
            "domain_rule": DOMAIN_RULE,
            "filter_out_of_domain": False,
        },
        "primary_metrics": [
            "intention_to_screen_top1_recovery",
            "intention_to_screen_top2_recall",
            "intention_to_screen_rank_of_best",
        ],
        "secondary_metrics": [
            "spearman_rho",
            "ndcg",
            "mae_kcal_mol",
            "rmse_kcal_mol",
            "selection_gain_over_panel_mean",
            "selection_enrichment_ratio",
            "prediction_interval_coverage",
            "applicability_domain_error_strata",
        ],
        "failed_split_policy": {
            "top1_recovery": 0.0,
            "top2_recall": 0.0,
            "rank_of_best": CANDIDATE_COUNT + 1,
            "other_metrics": None,
        },
        "random_baseline": {
            "method": "exact enumeration of all 6 rankings per candidate panel",
            "random_seed": None,
        },
    }


def prepare_design(output_dir: Path) -> dict[str, Any]:
    rows = read_allowed_metadata()
    splits = make_splits(rows)
    output_dir.mkdir(parents=True, exist_ok=True)
    atomic_write_text(output_dir / "DESIGN.md", design_document())
    atomic_write_json(output_dir / "design.json", design_spec(splits))
    atomic_write_json(output_dir / "split_definitions.json", splits)

    split_rows = []
    row_by_id = {row["Reaction_ID"]: row for row in rows}
    for split in splits:
        held = set(split["held_out_ids"])
        for row in rows:
            split_rows.append(
                {
                    "split_id": split["split_id"],
                    "stratum": split["stratum"],
                    "reaction_id": row["Reaction_ID"],
                    "source_id": row["Source_ID"],
                    "ligand_group": row["Ligand_Group"],
                    "role": "held_out" if row["Reaction_ID"] in held else "train",
                }
            )
    atomic_write_csv(
        output_dir / "splits.csv",
        ("split_id", "stratum", "reaction_id", "source_id", "ligand_group", "role"),
        split_rows,
    )
    atomic_write_json(
        output_dir / "seeds.json",
        {
            "base_seed": BASE_SEED,
            "purpose": "StericX bootstrap and response-permutation resampling only",
            "split_seeds": {split["split_id"]: split["model_seed"] for split in splits},
            "split_generation_seed": None,
            "split_generation_note": "splits are exhaustive and non-random",
            "random_baseline_seed": None,
            "random_baseline_note": "all six rankings are enumerated exactly",
            "ranking_tie_break": "numeric ligand identifier ascending",
        },
    )
    lock = {
        "status": "locked_before_study_011_outcome_revelation",
        "files": {name: sha256_file(output_dir / name) for name in DESIGN_FILES},
        "source_inputs": {
            str(METADATA.relative_to(ROOT)): sha256_file(METADATA),
            str(CATALOG.relative_to(ROOT)): sha256_file(CATALOG),
            str(CATALOG_PROVENANCE.relative_to(ROOT)): sha256_file(CATALOG_PROVENANCE),
            str(SIGPACK.relative_to(ROOT)): sha256_file(SIGPACK),
            str(Path(__file__).resolve().relative_to(ROOT)): sha256_file(
                Path(__file__).resolve()
            ),
        },
    }
    atomic_write_json(output_dir / "design_lock.json", lock)
    # An accidental caller-side mutation cannot affect what was serialized.
    if any(
        row_by_id[value]["Ligand_Group"] in split["held_out_groups"]
        for split in splits
        for value in split["training_ids"]
    ):
        raise AssertionError("design contains group leakage")
    return lock


def verify_design_lock(output_dir: Path) -> dict[str, Any]:
    lock_path = output_dir / "design_lock.json"
    if not lock_path.is_file():
        raise FileNotFoundError(
            f"missing {lock_path}; run this driver with --prepare-only first"
        )
    lock = json.loads(lock_path.read_text(encoding="utf-8"))
    errors = []
    for name, expected in lock["files"].items():
        path = output_dir / name
        actual = sha256_file(path) if path.is_file() else "missing"
        if actual != expected:
            errors.append(f"{name}: expected {expected}, found {actual}")
    script_key = str(Path(__file__).resolve().relative_to(ROOT))
    expected_script = lock["source_inputs"].get(script_key)
    actual_script = sha256_file(Path(__file__).resolve())
    if expected_script != actual_script:
        amendment_path = output_dir / "implementation_amendments.json"
        amendment = (
            json.loads(amendment_path.read_text(encoding="utf-8"))
            if amendment_path.is_file()
            else {}
        )
        accepted = (
            amendment.get("original_design_lock_sha256") == sha256_file(lock_path)
            and amendment.get("original_driver_sha256") == expected_script
            and amendment.get("amended_driver_sha256") == actual_script
            and amendment.get("changes_analysis_design") is False
        )
        if not accepted:
            errors.append(
                f"{script_key}: design was locked against {expected_script}, "
                f"found {actual_script}"
            )
    for source in (METADATA, CATALOG, CATALOG_PROVENANCE, SIGPACK):
        key = str(source.relative_to(ROOT))
        expected = lock["source_inputs"].get(key)
        actual = sha256_file(source)
        if expected != actual:
            errors.append(f"{key}: expected {expected}, found {actual}")
    if errors:
        raise RuntimeError("design lock verification failed:\n  " + "\n  ".join(errors))
    return lock


def run_command(command: Sequence[str]) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        command,
        cwd=ROOT,
        check=True,
        capture_output=True,
        text=True,
    )


def write_split_metadata(
    path: Path,
    rows: Sequence[dict[str, str]],
    held_out: set[str],
) -> None:
    atomic_write_csv(
        path,
        ("Reaction_ID", "Dataset_Split", "Ligand_Group"),
        (
            {
                "Reaction_ID": row["Reaction_ID"],
                "Dataset_Split": (
                    "blind" if row["Reaction_ID"] in held_out else "train"
                ),
                "Ligand_Group": row["Ligand_Group"],
            }
            for row in rows
        ),
    )


def write_candidate_library(
    path: Path,
    rows: Sequence[dict[str, str]],
    held_out: set[str],
) -> None:
    fields = (
        "Reaction_ID",
        "Ligand_SMILES",
        "Ligand_XYZ_Path",
        "NBO_Charge",
        "IR_Frequency",
    )
    selected = []
    for row in rows:
        if row["Reaction_ID"] not in held_out:
            continue
        geometry = Path(row["Ligand_XYZ_Path"])
        if not geometry.is_absolute():
            geometry = ROOT / "data" / geometry
        selected.append(
            {
                **{field: row[field] for field in fields if field != "Ligand_XYZ_Path"},
                "Ligand_XYZ_Path": str(geometry.resolve()),
            }
        )
    atomic_write_csv(path, fields, selected)


def compact_screen_report(report: dict[str, Any], split_id: str) -> dict[str, Any]:
    """Retain complete screening decisions without volatile temporary paths."""
    keys = (
        "model",
        "selected_features",
        "required_inputs",
        "training_count",
        "training_r2",
        "training_rmse_kcal_mol",
        "library_size",
        "screened",
        "returned",
        "ranking_order",
        "model_optimization",
        "ranking_overridden",
        "domain_rule",
        "domain_rule_description",
        "domain_threshold",
        "domain_filter_applied",
        "skipped",
        "excluded",
        "inside_domain",
        "warning_leverage",
        "high_leverage",
        "trust_summary",
        "domain_summary",
        "neighbor_rule",
        "uncertainty_note",
        "loo_q2",
        "loo_rmse_kcal_mol",
        "group_loo_q2",
        "hits",
    )
    return {"split_id": split_id, **{key: report.get(key) for key in keys}}


def screen_all_splits(
    binary: Path,
    output_dir: Path,
    rows: Sequence[dict[str, str]],
    splits: Sequence[dict[str, Any]],
) -> tuple[
    list[dict[str, Any]],
    list[dict[str, Any]],
    list[dict[str, Any]],
    list[dict[str, Any]],
]:
    """Fit and screen every split without reading held-out target values."""
    frozen_rows: list[dict[str, Any]] = []
    model_reports: list[dict[str, Any]] = []
    screen_reports: list[dict[str, Any]] = []
    statuses: list[dict[str, Any]] = []
    for split in splits:
        split_id = split["split_id"]
        held_out = set(split["held_out_ids"])
        status = {
            "split_id": split_id,
            "stratum": split["stratum"],
            "status": "failed",
            "phase": "setup",
            "error": "",
            "training_count": 8,
            "candidate_count": 3,
        }
        try:
            with tempfile.TemporaryDirectory(prefix=f"stericx-{split_id}-") as raw_work:
                work = Path(raw_work)
                metadata_path = work / "metadata.csv"
                candidates_path = work / "candidates.csv"
                model_path = work / "model.json"
                frozen_path = work / "fit_frozen_predictions.csv"
                portable_path = work / "portable_model.json"
                write_split_metadata(metadata_path, rows, held_out)
                write_candidate_library(candidates_path, rows, held_out)

                status["phase"] = "fit"
                run_command(
                    (
                        str(binary),
                        "fit",
                        "--data",
                        str(SIGPACK),
                        "--metadata",
                        str(metadata_path),
                        "--output",
                        str(model_path),
                        "--predictions",
                        str(frozen_path),
                        "--max-terms",
                        str(MAX_TERMS),
                        "--bootstrap",
                        str(BOOTSTRAP_SAMPLES),
                        "--permutations",
                        str(PERMUTATION_SAMPLES),
                        "--seed",
                        str(split["model_seed"]),
                        "--portable-model",
                        str(portable_path),
                        "--model-id",
                        f"study-011-{split_id.lower()}",
                        "--reaction-family",
                        "Ni-catalyzed homo-Diels-Alder",
                        "--catalyst-metal",
                        "Ni",
                        "--ligand-class",
                        "monodentate phosphorus(III)",
                        "--source-url",
                        SOURCE_URL,
                        "--response-temp-k",
                        "298.15",
                        "--optimize",
                        "maximize",
                    )
                )
                model_report = json.loads(model_path.read_text(encoding="utf-8"))
                fit_predictions = {
                    row["Reaction_ID"]: float(row["Predicted_ddG_kcal_mol"])
                    for row in read_csv(frozen_path)
                }
                model_reports.append(
                    {
                        "split_id": split_id,
                        "seed": split["model_seed"],
                        "report": model_report,
                    }
                )

                status["phase"] = "screen"
                screened = run_command(
                    (
                        str(binary),
                        "screen",
                        str(portable_path),
                        str(candidates_path),
                        "--format",
                        "json",
                        "--domain-rule",
                        DOMAIN_RULE,
                        "--sterimol-axis",
                        "bond",
                    )
                )
                report = json.loads(screened.stdout)
                if (
                    report.get("screened") != CANDIDATE_COUNT
                    or len(report["hits"]) != CANDIDATE_COUNT
                ):
                    raise RuntimeError(
                        f"expected {CANDIDATE_COUNT} ranked candidates, got "
                        f"screened={report.get('screened')} "
                        f"hits={len(report.get('hits', []))}"
                    )
                if report.get("skipped") != 0:
                    raise RuntimeError(
                        f"screen skipped {report['skipped']} candidate(s)"
                    )
                hit_ids = {hit["ligand"] for hit in report["hits"]}
                if hit_ids != held_out:
                    raise RuntimeError(
                        f"screened candidate IDs {sorted(hit_ids)} do not match "
                        f"holdout {sorted(held_out)}"
                    )
                screen_reports.append(compact_screen_report(report, split_id))
                for hit in report["hits"]:
                    fit_prediction = fit_predictions[hit["ligand"]]
                    frozen_rows.append(
                        {
                            "split_id": split_id,
                            "stratum": split["stratum"],
                            "reaction_id": hit["ligand"],
                            "predicted_rank": hit["rank"],
                            "predicted_ddg_kcal_mol": hit["predicted_ddg_kcal_mol"],
                            "fit_frozen_prediction_kcal_mol": fit_prediction,
                            "screen_minus_fit_prediction_kcal_mol": (
                                hit["predicted_ddg_kcal_mol"] - fit_prediction
                            ),
                            "domain_verdict": hit["domain_verdict"],
                            "applicability": hit["applicability"],
                            "trust": hit["trust"],
                            "nearest_training_distance": hit[
                                "nearest_training_distance"
                            ],
                            "nearest_training_ligand": hit["nearest_training_ligand"],
                            "nearest_training_threshold": hit[
                                "nearest_training_threshold"
                            ],
                            "nearest_training_ratio": hit["nearest_training_ratio"],
                            "mahalanobis_distance": hit["mahalanobis_distance"],
                            "maximum_extrapolation": hit["maximum_extrapolation"],
                            "leverage": hit["leverage"],
                            "leverage_ratio": hit["leverage_ratio"],
                            "prediction_interval_low": hit["prediction_interval_low"],
                            "prediction_interval_high": hit["prediction_interval_high"],
                            "coefficient_band_low": hit["coefficient_band_low"],
                            "coefficient_band_high": hit["coefficient_band_high"],
                            "bootstrap_mean_low": (
                                hit["uncertainty"]["lower"]
                                if hit["uncertainty"]
                                else None
                            ),
                            "bootstrap_mean_high": (
                                hit["uncertainty"]["upper"]
                                if hit["uncertainty"]
                                else None
                            ),
                            "descriptors_json": json.dumps(
                                hit["descriptors"],
                                sort_keys=True,
                                separators=(",", ":"),
                            ),
                            "selected_features": "|".join(report["selected_features"]),
                            "model_seed": split["model_seed"],
                        }
                    )
                status.update(status="success", phase="complete")
        except (
            OSError,
            ValueError,
            KeyError,
            RuntimeError,
            subprocess.CalledProcessError,
        ) as error:
            if isinstance(error, subprocess.CalledProcessError):
                detail = "\n".join(
                    part.strip()
                    for part in (error.stdout or "", error.stderr or "")
                    if part.strip()
                )
                status["error"] = f"exit {error.returncode}: {detail}"
            else:
                status["error"] = str(error)
        statuses.append(status)
    return frozen_rows, model_reports, screen_reports, statuses


FROZEN_FIELDS: Final[tuple[str, ...]] = (
    "split_id",
    "stratum",
    "reaction_id",
    "predicted_rank",
    "predicted_ddg_kcal_mol",
    "fit_frozen_prediction_kcal_mol",
    "screen_minus_fit_prediction_kcal_mol",
    "domain_verdict",
    "applicability",
    "trust",
    "nearest_training_distance",
    "nearest_training_ligand",
    "nearest_training_threshold",
    "nearest_training_ratio",
    "mahalanobis_distance",
    "maximum_extrapolation",
    "leverage",
    "leverage_ratio",
    "prediction_interval_low",
    "prediction_interval_high",
    "coefficient_band_low",
    "coefficient_band_high",
    "bootstrap_mean_low",
    "bootstrap_mean_high",
    "descriptors_json",
    "selected_features",
    "model_seed",
)


def write_frozen_artifacts(
    output_dir: Path,
    frozen_rows: Sequence[dict[str, Any]],
    model_reports: Sequence[dict[str, Any]],
    screen_reports: Sequence[dict[str, Any]],
    statuses: Sequence[dict[str, Any]],
) -> str:
    frozen_path = output_dir / "frozen_predictions.csv"
    atomic_write_csv(frozen_path, FROZEN_FIELDS, frozen_rows)
    frozen_digest = sha256_file(frozen_path)
    atomic_write_text(
        output_dir / "frozen_predictions.sha256",
        f"{frozen_digest}  frozen_predictions.csv\n",
    )
    atomic_write_text(
        output_dir / "model_reports.jsonl",
        "".join(json.dumps(row, sort_keys=True) + "\n" for row in model_reports),
    )
    atomic_write_text(
        output_dir / "screen_reports.jsonl",
        "".join(json.dumps(row, sort_keys=True) + "\n" for row in screen_reports),
    )
    atomic_write_csv(
        output_dir / "split_status.csv",
        (
            "split_id",
            "stratum",
            "status",
            "phase",
            "error",
            "training_count",
            "candidate_count",
        ),
        statuses,
    )
    return frozen_digest


def reveal_outcomes(rows: Sequence[dict[str, str]]) -> dict[str, float]:
    """First study function that accesses the held-out outcome column."""
    wanted = {int(row["Source_ID"]): row["Reaction_ID"] for row in rows}
    outcomes: dict[str, float] = {}
    with CATALOG.open(newline="", encoding="utf-8") as handle:
        reader = csv.reader(handle)
        header = next(reader)
        outcome_column = header.index("ddG_abs")
        for record in reader:
            try:
                source_id = int(record[0])
            except (ValueError, IndexError):
                continue
            if source_id in wanted:
                value = record[outcome_column].strip()
                if not value:
                    raise ValueError(
                        f"source ligand {source_id} has no ddG_abs outcome"
                    )
                outcomes[wanted[source_id]] = float(value)
    if set(outcomes) != set(wanted.values()):
        missing = sorted(set(wanted.values()).difference(outcomes))
        raise ValueError(f"catalog outcomes missing for {missing}")
    if any(not math.isfinite(value) or value < 0.0 for value in outcomes.values()):
        raise ValueError("ddG_abs outcomes must be finite and non-negative")
    return outcomes


def safe_spearman(
    predicted: Sequence[float], observed: Sequence[float]
) -> float | None:
    left = np.asarray(predicted, dtype=float)
    right = np.asarray(observed, dtype=float)
    if left.size < 3 or np.ptp(left) == 0.0 or np.ptp(right) == 0.0:
        return None
    left_rank = rankdata(left, method="average")
    right_rank = rankdata(right, method="average")
    value = float(np.corrcoef(left_rank, right_rank)[0, 1])
    return value if math.isfinite(value) else None


def ndcg(predicted_order: Sequence[str], outcomes: dict[str, float]) -> float | None:
    discounts = np.log2(np.arange(2, len(predicted_order) + 2, dtype=float))
    observed = np.asarray([outcomes[value] for value in predicted_order], dtype=float)
    ideal = np.asarray(sorted(observed, reverse=True), dtype=float)
    denominator = float(np.sum(ideal / discounts))
    return (
        float(np.sum(observed / discounts) / denominator) if denominator > 0.0 else None
    )


def ranking_metrics(
    ordered_ids: Sequence[str],
    predicted_values: Sequence[float] | None,
    outcomes: dict[str, float],
) -> dict[str, Any]:
    actual_order = sorted(
        ordered_ids, key=lambda value: (-outcomes[value], numeric_id(value))
    )
    positions = {value: index + 1 for index, value in enumerate(ordered_ids)}
    top1 = float(ordered_ids[0] == actual_order[0])
    top2 = len(set(ordered_ids[:2]).intersection(actual_order[:2])) / 2.0
    candidate_mean = float(np.mean([outcomes[value] for value in ordered_ids]))
    top1_outcome = outcomes[ordered_ids[0]]
    top2_outcome = float(np.mean([outcomes[value] for value in ordered_ids[:2]]))
    score_values = (
        list(predicted_values)
        if predicted_values is not None
        else [float(len(ordered_ids) - index) for index in range(len(ordered_ids))]
    )
    observed_values = [outcomes[value] for value in ordered_ids]
    return {
        "top1_recovery": top1,
        "top2_recall": top2,
        "rank_of_best": positions[actual_order[0]],
        "spearman_rho": safe_spearman(score_values, observed_values),
        "ndcg": ndcg(ordered_ids, outcomes),
        "candidate_mean_ddg": candidate_mean,
        "selected_top1_ddg": top1_outcome,
        "selected_top2_mean_ddg": top2_outcome,
        "top1_gain_over_panel_mean": top1_outcome - candidate_mean,
        "top2_gain_over_panel_mean": top2_outcome - candidate_mean,
        "top1_outcome_enrichment": (
            top1_outcome / candidate_mean if candidate_mean != 0.0 else None
        ),
        "top2_outcome_enrichment": (
            top2_outcome / candidate_mean if candidate_mean != 0.0 else None
        ),
    }


def evaluate(
    frozen_rows: Sequence[dict[str, Any]],
    statuses: Sequence[dict[str, Any]],
    splits: Sequence[dict[str, Any]],
    outcomes: dict[str, float],
) -> tuple[
    list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]], dict[str, Any]
]:
    frozen_by_split: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in frozen_rows:
        frozen_by_split[row["split_id"]].append(dict(row))
    status_by_split = {row["split_id"]: row for row in statuses}
    revealed_rows: list[dict[str, Any]] = []
    split_metrics: list[dict[str, Any]] = []
    random_rows: list[dict[str, Any]] = []

    for split in splits:
        split_id = split["split_id"]
        status = status_by_split[split_id]
        successful = status["status"] == "success"
        screen_rows = sorted(
            frozen_by_split.get(split_id, []),
            key=lambda row: int(row["predicted_rank"]),
        )
        if successful and len(screen_rows) != CANDIDATE_COUNT:
            raise ValueError(
                f"{split_id} is marked successful without three predictions"
            )
        if successful:
            ordered_ids = [row["reaction_id"] for row in screen_rows]
            predicted = [float(row["predicted_ddg_kcal_mol"]) for row in screen_rows]
            fold = ranking_metrics(ordered_ids, predicted, outcomes)
            residuals = [
                pred - outcomes[value]
                for value, pred in zip(ordered_ids, predicted, strict=True)
            ]
            fold["mae_kcal_mol"] = float(np.mean(np.abs(residuals)))
            fold["rmse_kcal_mol"] = float(np.sqrt(np.mean(np.square(residuals))))
            actual_order = sorted(
                ordered_ids, key=lambda value: (-outcomes[value], numeric_id(value))
            )
            actual_ranks = {
                value: index + 1 for index, value in enumerate(actual_order)
            }
            for row, residual in zip(screen_rows, residuals, strict=True):
                low = optional_float(row.get("prediction_interval_low"))
                high = optional_float(row.get("prediction_interval_high"))
                observed = outcomes[row["reaction_id"]]
                revealed_rows.append(
                    {
                        **row,
                        "observed_ddg_kcal_mol": observed,
                        "observed_rank": actual_ranks[row["reaction_id"]],
                        "residual_kcal_mol": residual,
                        "absolute_error_kcal_mol": abs(residual),
                        "prediction_interval_covers": (
                            None
                            if low is None or high is None
                            else int(low <= observed <= high)
                        ),
                        "is_observed_best": int(actual_ranks[row["reaction_id"]] == 1),
                    }
                )
        else:
            fold = {
                "top1_recovery": 0.0,
                "top2_recall": 0.0,
                "rank_of_best": CANDIDATE_COUNT + 1,
                "spearman_rho": None,
                "ndcg": None,
                "candidate_mean_ddg": float(
                    np.mean([outcomes[value] for value in split["held_out_ids"]])
                ),
                "selected_top1_ddg": None,
                "selected_top2_mean_ddg": None,
                "top1_gain_over_panel_mean": None,
                "top2_gain_over_panel_mean": None,
                "top1_outcome_enrichment": None,
                "top2_outcome_enrichment": None,
                "mae_kcal_mol": None,
                "rmse_kcal_mol": None,
            }
        split_metrics.append(
            {
                "split_id": split_id,
                "stratum": split["stratum"],
                "status": status["status"],
                **fold,
            }
        )

        candidate_ids = sorted(split["held_out_ids"], key=numeric_id)
        for permutation_index, ordering in enumerate(
            itertools.permutations(candidate_ids), start=1
        ):
            baseline = ranking_metrics(ordering, None, outcomes)
            random_rows.append(
                {
                    "split_id": split_id,
                    "stratum": split["stratum"],
                    "permutation": permutation_index,
                    "ranking": "|".join(ordering),
                    **baseline,
                }
            )

    metrics = aggregate_metrics(
        split_metrics,
        revealed_rows,
        random_rows,
        statuses,
        splits,
    )
    return revealed_rows, split_metrics, random_rows, metrics


def optional_float(value: Any) -> float | None:
    if value is None or value == "":
        return None
    result = float(value)
    return result if math.isfinite(result) else None


def mean_defined(rows: Sequence[dict[str, Any]], key: str) -> float | None:
    values = [optional_float(row.get(key)) for row in rows]
    finite = [value for value in values if value is not None]
    return float(np.mean(finite)) if finite else None


def divide(left: float | None, right: float | None) -> float | None:
    if left is None or right is None or right == 0.0:
        return None
    return left / right


def error_summary(rows: Sequence[dict[str, Any]]) -> dict[str, Any]:
    residuals = np.asarray(
        [float(row["residual_kcal_mol"]) for row in rows], dtype=float
    )
    coverage_values = [
        int(row["prediction_interval_covers"])
        for row in rows
        if row.get("prediction_interval_covers") is not None
    ]
    return {
        "count": len(rows),
        "mae_kcal_mol": (float(np.mean(np.abs(residuals))) if residuals.size else None),
        "rmse_kcal_mol": (
            float(np.sqrt(np.mean(np.square(residuals)))) if residuals.size else None
        ),
        "mean_error_kcal_mol": float(np.mean(residuals)) if residuals.size else None,
        "prediction_interval_count": len(coverage_values),
        "prediction_interval_coverage": (
            float(np.mean(coverage_values)) if coverage_values else None
        ),
    }


def aggregate_metrics(
    split_metrics: Sequence[dict[str, Any]],
    revealed_rows: Sequence[dict[str, Any]],
    random_rows: Sequence[dict[str, Any]],
    statuses: Sequence[dict[str, Any]],
    splits: Sequence[dict[str, Any]],
) -> dict[str, Any]:
    successful = [row for row in split_metrics if row["status"] == "success"]
    planned = len(splits)
    random_top1 = mean_defined(random_rows, "top1_recovery")
    random_top2 = mean_defined(random_rows, "top2_recall")
    primary_top1 = mean_defined(split_metrics, "top1_recovery")
    primary_top2 = mean_defined(split_metrics, "top2_recall")

    by_verdict: dict[str, list[dict[str, Any]]] = defaultdict(list)
    by_range: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in revealed_rows:
        by_verdict[row["domain_verdict"]].append(row)
        by_range[row["applicability"]].append(row)
    distance_rows = [
        row
        for row in revealed_rows
        if optional_float(row.get("nearest_training_ratio")) is not None
    ]
    leverage_rows = [
        row
        for row in revealed_rows
        if optional_float(row.get("leverage_ratio")) is not None
    ]
    domain = {
        "domain_rule": DOMAIN_RULE,
        "by_domain_verdict": {
            key: error_summary(value) for key, value in sorted(by_verdict.items())
        },
        "by_feature_range_status": {
            key: error_summary(value) for key, value in sorted(by_range.items())
        },
        "spearman_absolute_error_vs_nearest_training_ratio": safe_spearman(
            [float(row["absolute_error_kcal_mol"]) for row in distance_rows],
            [float(row["nearest_training_ratio"]) for row in distance_rows],
        ),
        "spearman_absolute_error_vs_leverage_ratio": safe_spearman(
            [float(row["absolute_error_kcal_mol"]) for row in leverage_rows],
            [float(row["leverage_ratio"]) for row in leverage_rows],
        ),
    }

    stratum_metrics = {}
    for stratum in sorted({row["stratum"] for row in split_metrics}):
        subset = [row for row in split_metrics if row["stratum"] == stratum]
        subset_success = [row for row in subset if row["status"] == "success"]
        stratum_metrics[stratum] = {
            "planned_splits": len(subset),
            "successful_splits": len(subset_success),
            "intention_to_screen_top1_recovery": mean_defined(subset, "top1_recovery"),
            "intention_to_screen_top2_recall": mean_defined(subset, "top2_recall"),
            "intention_to_screen_mean_rank_of_best": mean_defined(
                subset, "rank_of_best"
            ),
            "successful_only_mean_spearman_rho": mean_defined(
                subset_success, "spearman_rho"
            ),
        }

    max_screen_fit_delta = max(
        (
            abs(float(row["screen_minus_fit_prediction_kcal_mol"]))
            for row in revealed_rows
        ),
        default=None,
    )
    return {
        "study": "study_011",
        "result_status": (
            "complete_all_splits"
            if len(successful) == planned
            else "complete_with_failures"
        ),
        "execution": {
            "planned_splits": planned,
            "successful_splits": len(successful),
            "failed_splits": planned - len(successful),
            "failure_details": [
                {
                    "split_id": row["split_id"],
                    "phase": row["phase"],
                    "error": row["error"],
                }
                for row in statuses
                if row["status"] != "success"
            ],
            "ranked_predictions": len(revealed_rows),
            "maximum_absolute_screen_minus_fit_prediction_kcal_mol": (
                max_screen_fit_delta
            ),
        },
        "primary_intention_to_screen": {
            "top1_recovery": primary_top1,
            "top2_recall": primary_top2,
            "mean_rank_of_best": mean_defined(split_metrics, "rank_of_best"),
            "top1_enrichment_over_exact_random": divide(primary_top1, random_top1),
            "top2_enrichment_over_exact_random": divide(primary_top2, random_top2),
            "failed_split_scoring": {
                "top1_recovery": 0.0,
                "top2_recall": 0.0,
                "rank_of_best": CANDIDATE_COUNT + 1,
            },
        },
        "successful_ranking_only": {
            "splits": len(successful),
            "top1_recovery": mean_defined(successful, "top1_recovery"),
            "top2_recall": mean_defined(successful, "top2_recall"),
            "mean_rank_of_best": mean_defined(successful, "rank_of_best"),
            "mean_spearman_rho": mean_defined(successful, "spearman_rho"),
            "mean_ndcg": mean_defined(successful, "ndcg"),
            "mean_top1_gain_over_panel_mean_kcal_mol": mean_defined(
                successful, "top1_gain_over_panel_mean"
            ),
            "mean_top2_gain_over_panel_mean_kcal_mol": mean_defined(
                successful, "top2_gain_over_panel_mean"
            ),
            "mean_top1_outcome_enrichment": mean_defined(
                successful, "top1_outcome_enrichment"
            ),
            "mean_top2_outcome_enrichment": mean_defined(
                successful, "top2_outcome_enrichment"
            ),
        },
        "prediction_error": {
            **error_summary(revealed_rows),
            "macro_mean_split_mae_kcal_mol": mean_defined(successful, "mae_kcal_mol"),
            "macro_mean_split_rmse_kcal_mol": mean_defined(successful, "rmse_kcal_mol"),
            "note": "repeated ligand predictions are not statistically independent",
        },
        "exact_random_baseline": {
            "rankings_per_split": math.factorial(CANDIDATE_COUNT),
            "evaluated_rankings": len(random_rows),
            "top1_recovery": random_top1,
            "top2_recall": random_top2,
            "mean_rank_of_best": mean_defined(random_rows, "rank_of_best"),
            "mean_spearman_rho": mean_defined(random_rows, "spearman_rho"),
            "mean_ndcg": mean_defined(random_rows, "ndcg"),
            "mean_top1_gain_over_panel_mean_kcal_mol": mean_defined(
                random_rows, "top1_gain_over_panel_mean"
            ),
            "mean_top2_gain_over_panel_mean_kcal_mol": mean_defined(
                random_rows, "top2_gain_over_panel_mean"
            ),
        },
        "by_split_stratum": stratum_metrics,
        "applicability_domain": domain,
        "interpretation_guardrails": [
            "The 57 exhaustive folds overlap and are not independent replicates.",
            "The candidate panels contain three ligands, not the full Kraken library.",
            "The one three-member scaffold contributes one panel; singleton "
            "ligands contribute 21 panels each.",
            "Pooled error metrics repeatedly count the same experimental ligand.",
        ],
    }


RANKING_FIELDS: Final[tuple[str, ...]] = (
    *FROZEN_FIELDS,
    "observed_ddg_kcal_mol",
    "observed_rank",
    "residual_kcal_mol",
    "absolute_error_kcal_mol",
    "prediction_interval_covers",
    "is_observed_best",
)

SPLIT_METRIC_FIELDS: Final[tuple[str, ...]] = (
    "split_id",
    "stratum",
    "status",
    "top1_recovery",
    "top2_recall",
    "rank_of_best",
    "spearman_rho",
    "ndcg",
    "mae_kcal_mol",
    "rmse_kcal_mol",
    "candidate_mean_ddg",
    "selected_top1_ddg",
    "selected_top2_mean_ddg",
    "top1_gain_over_panel_mean",
    "top2_gain_over_panel_mean",
    "top1_outcome_enrichment",
    "top2_outcome_enrichment",
)

RANDOM_FIELDS: Final[tuple[str, ...]] = (
    "split_id",
    "stratum",
    "permutation",
    "ranking",
    "top1_recovery",
    "top2_recall",
    "rank_of_best",
    "spearman_rho",
    "ndcg",
    "candidate_mean_ddg",
    "selected_top1_ddg",
    "selected_top2_mean_ddg",
    "top1_gain_over_panel_mean",
    "top2_gain_over_panel_mean",
    "top1_outcome_enrichment",
    "top2_outcome_enrichment",
)


def plot_ranking_performance(
    path: Path,
    metrics: dict[str, Any],
    split_metrics: Sequence[dict[str, Any]],
    rankings: Sequence[dict[str, Any]],
) -> None:
    plt.style.use("seaborn-v0_8-whitegrid")
    figure, axes = plt.subplots(2, 2, figsize=(11, 8.5), constrained_layout=True)
    primary = metrics["primary_intention_to_screen"]
    random = metrics["exact_random_baseline"]

    names = ["Top-1", "Top-2"]
    model_values = [primary["top1_recovery"], primary["top2_recall"]]
    random_values = [random["top1_recovery"], random["top2_recall"]]
    x = np.arange(2)
    width = 0.36
    axes[0, 0].bar(x - width / 2, model_values, width, label="StericX", color="#176B87")
    axes[0, 0].bar(
        x + width / 2, random_values, width, label="Exact random", color="#D8A03C"
    )
    axes[0, 0].set_xticks(x, names)
    axes[0, 0].set_ylim(0, 1.05)
    axes[0, 0].set_ylabel("Recovery fraction")
    axes[0, 0].set_title("Intention-to-screen recovery")
    axes[0, 0].legend(frameon=False)

    ranks = [int(row["rank_of_best"]) for row in split_metrics]
    bins = np.arange(0.5, 5.5, 1.0)
    axes[0, 1].hist(ranks, bins=bins, density=True, color="#176B87", alpha=0.85)
    axes[0, 1].plot(
        [1, 2, 3], [1 / 3] * 3, "o--", color="#D8A03C", label="Exact random"
    )
    axes[0, 1].axvline(3.5, color="#7A7A7A", linewidth=1)
    axes[0, 1].set_xticks([1, 2, 3, 4], ["1", "2", "3", "4 (failed)"])
    axes[0, 1].set_ylabel("Fraction of planned splits")
    axes[0, 1].set_title("Rank of observed-best ligand")
    axes[0, 1].legend(frameon=False)

    rho = [
        float(row["spearman_rho"])
        for row in split_metrics
        if row.get("spearman_rho") is not None
    ]
    if rho:
        jitter = np.linspace(-0.06, 0.06, len(rho))
        axes[1, 0].scatter(jitter, rho, color="#176B87", alpha=0.7, s=28)
        axes[1, 0].boxplot(
            rho,
            positions=[0],
            widths=[0.32],
            showfliers=False,
            patch_artist=False,
        )
    axes[1, 0].axhline(0.0, color="#D8A03C", linestyle="--", label="Random mean")
    axes[1, 0].set_xlim(-0.5, 0.5)
    axes[1, 0].set_xticks([0], ["Successful splits"])
    axes[1, 0].set_ylim(-1.08, 1.08)
    axes[1, 0].set_ylabel("Spearman rho")
    axes[1, 0].set_title("Within-panel rank correlation")
    axes[1, 0].legend(frameon=False)

    if rankings:
        predicted = np.asarray(
            [float(row["predicted_ddg_kcal_mol"]) for row in rankings]
        )
        observed = np.asarray([float(row["observed_ddg_kcal_mol"]) for row in rankings])
        colors = [
            "#176B87" if row["domain_verdict"] == "interpolation" else "#C44E52"
            for row in rankings
        ]
        axes[1, 1].scatter(observed, predicted, c=colors, alpha=0.55, s=25)
        low = min(float(observed.min()), float(predicted.min()))
        high = max(float(observed.max()), float(predicted.max()))
        pad = max((high - low) * 0.05, 0.02)
        axes[1, 1].plot(
            [low - pad, high + pad], [low - pad, high + pad], "--", color="#555555"
        )
        axes[1, 1].set_xlim(low - pad, high + pad)
        axes[1, 1].set_ylim(low - pad, high + pad)
    axes[1, 1].set_xlabel(r"Observed $|\Delta\Delta G^\ddagger|$ (kcal mol$^{-1}$)")
    axes[1, 1].set_ylabel(r"Predicted $|\Delta\Delta G^\ddagger|$ (kcal mol$^{-1}$)")
    axes[1, 1].set_title("All repeated held-out predictions")
    figure.suptitle("Study 011 — retrospective ligand ranking", fontsize=15)
    figure.savefig(path, dpi=180)
    plt.close(figure)


def plot_applicability(path: Path, rankings: Sequence[dict[str, Any]]) -> None:
    plt.style.use("seaborn-v0_8-whitegrid")
    figure, axes = plt.subplots(1, 2, figsize=(11, 4.6), constrained_layout=True)
    with_distance = [
        row
        for row in rankings
        if optional_float(row.get("nearest_training_ratio")) is not None
    ]
    verdicts = sorted({row["domain_verdict"] for row in rankings})
    palette = {
        "interpolation": "#176B87",
        "sparse_interpolation": "#D8A03C",
        "extrapolation": "#C44E52",
        "unknown": "#7A7A7A",
    }
    for verdict in verdicts:
        subset = [row for row in with_distance if row["domain_verdict"] == verdict]
        axes[0].scatter(
            [float(row["nearest_training_ratio"]) for row in subset],
            [float(row["absolute_error_kcal_mol"]) for row in subset],
            label=f"{verdict} (n={len(subset)})",
            color=palette.get(verdict, "#7A7A7A"),
            alpha=0.65,
            s=28,
        )
    axes[0].axvline(1.0, color="#555555", linestyle="--", linewidth=1)
    axes[0].set_xlabel("Nearest-training distance / domain threshold")
    axes[0].set_ylabel(r"Absolute error (kcal mol$^{-1}$)")
    axes[0].set_title("Error versus training-set proximity")
    if verdicts:
        axes[0].legend(frameon=False, fontsize=8)

    box_data = [
        [
            float(row["absolute_error_kcal_mol"])
            for row in rankings
            if row["domain_verdict"] == verdict
        ]
        for verdict in verdicts
    ]
    if box_data:
        axes[1].boxplot(box_data, tick_labels=verdicts, showfliers=True)
    axes[1].tick_params(axis="x", rotation=20)
    axes[1].set_ylabel(r"Absolute error (kcal mol$^{-1}$)")
    axes[1].set_title("Error by applicability verdict")
    figure.suptitle("Study 011 — applicability-domain behavior", fontsize=15)
    figure.savefig(path, dpi=180)
    plt.close(figure)


def fmt(value: Any, digits: int = 3) -> str:
    if value is None:
        return "not defined"
    return f"{float(value):.{digits}f}"


def write_study_report(
    output_dir: Path,
    metrics: dict[str, Any],
    frozen_digest: str,
    design_lock_digest: str,
) -> None:
    primary = metrics["primary_intention_to_screen"]
    successful = metrics["successful_ranking_only"]
    error = metrics["prediction_error"]
    random = metrics["exact_random_baseline"]
    execution = metrics["execution"]
    top1_model = fmt(primary["top1_recovery"])
    top1_random = fmt(random["top1_recovery"])
    top1_enrichment = fmt(primary["top1_enrichment_over_exact_random"])
    top2_model = fmt(primary["top2_recall"])
    top2_random = fmt(random["top2_recall"])
    top2_enrichment = fmt(primary["top2_enrichment_over_exact_random"])
    mean_best_rank = fmt(primary["mean_rank_of_best"])
    random_best_rank = fmt(random["mean_rank_of_best"])
    mean_rho = fmt(successful["mean_spearman_rho"])
    distance_rho = fmt(
        metrics["applicability_domain"][
            "spearman_absolute_error_vs_nearest_training_ratio"
        ]
    )
    leverage_rho = fmt(
        metrics["applicability_domain"]["spearman_absolute_error_vs_leverage_ratio"]
    )
    screen_fit_delta = fmt(
        execution["maximum_absolute_screen_minus_fit_prediction_kcal_mol"], 6
    )
    failures = execution["failure_details"]
    failure_text = (
        "No split failed to fit or screen."
        if not failures
        else "\n".join(
            f"- `{row['split_id']}` failed during `{row['phase']}`: {row['error']}"
            for row in failures
        )
    )
    strata = metrics["by_split_stratum"]
    stratum_rows = "\n".join(
        "| {name} | {planned} | {success} | {top1} | {top2} | {rho} |".format(
            name=name,
            planned=value["planned_splits"],
            success=value["successful_splits"],
            top1=fmt(value["intention_to_screen_top1_recovery"]),
            top2=fmt(value["intention_to_screen_top2_recall"]),
            rho=fmt(value["successful_only_mean_spearman_rho"]),
        )
        for name, value in strata.items()
    )
    domain_rows = "\n".join(
        "| {name} | {count} | {mae} | {rmse} | {coverage} |".format(
            name=name,
            count=value["count"],
            mae=fmt(value["mae_kcal_mol"]),
            rmse=fmt(value["rmse_kcal_mol"]),
            coverage=fmt(value["prediction_interval_coverage"]),
        )
        for name, value in metrics["applicability_domain"]["by_domain_verdict"].items()
    )
    conclusion = (
        "StericX exceeded the exact random baseline for top-1 recovery."
        if primary["top1_recovery"] > random["top1_recovery"]
        else "StericX did not exceed the exact random baseline for top-1 recovery."
    )
    report = f"""# Study 011: retrospective ligand ranking

## Result

{conclusion} Across all {execution["planned_splits"]} preregistered screens, top-1
recovery was **{top1_model}** versus **{top1_random}** for exact random selection
(enrichment **{top1_enrichment}x**). Top-2 recovery was **{top2_model}** versus
**{top2_random}** (enrichment **{top2_enrichment}x**), and the mean rank of the
best held-out ligand was **{mean_best_rank}** versus **{random_best_rank}** at
random. These are intention-to-screen
metrics: any failed screen receives the penalties fixed in `DESIGN.md`.

Of {execution["planned_splits"]} planned splits, **{execution["successful_splits"]}**
produced a complete three-ligand ranking and **{execution["failed_splits"]}** failed.
Among successful rankings, mean Spearman rho was **{mean_rho}**,
mean nDCG was **{fmt(successful["mean_ndcg"])}**, pooled MAE was
**{fmt(error["mae_kcal_mol"])} kcal/mol**, and pooled RMSE was
**{fmt(error["rmse_kcal_mol"])} kcal/mol**. The same ligand is predicted repeatedly,
so pooled errors and fold averages are descriptive rather than independent estimates.

![Ranking performance](ranking_performance.png)

## Predefined strata

| Holdout stratum | Planned | Successful | Top-1 recovery | Top-2 recall | Mean rho |
|---|---:|---:|---:|---:|---:|
{stratum_rows}

The single three-member-scaffold fold is shown separately because the exhaustive
design includes it once, while each singleton ligand appears in 21 distinct panels.
No stratum or split was selected after outcomes were revealed.

## Applicability-domain behavior

StericX predictions were ranked without filtering extrapolations. The continuous
nearest-training distance, leverage, range status, and categorical verdict for every
candidate are retained in `rankings.csv`.

| Domain verdict | Repeated predictions | MAE | RMSE | 95% PI coverage |
|---|---:|---:|---:|---:|
{domain_rows}

The descriptive Spearman association between absolute error and normalized nearest-
training distance was **{distance_rho}**; the association with leverage ratio was
**{leverage_rho}**.
Small and repeated samples prevent a calibrated claim about these domain strata.

![Applicability-domain behavior](applicability_domain.png)

## Method

The design was written and checksum-locked before this study accessed held-out
`ddG_abs` values. The design-lock SHA-256 is `{design_lock_digest}`. Exact split
membership is in `split_definitions.json` and `splits.csv`; seeds are in `seeds.json`.

For each split, the driver changed only row-aligned split labels, ran `stericx fit`
on the eight allowed observations, built a target-free three-ligand candidate CSV,
then called `stericx screen`. All predictions were completed before target reveal.
The frozen-prediction SHA-256 is `{frozen_digest}`. The greatest absolute difference
between `fit`'s frozen prediction and the actual `screen` prediction was
**{screen_fit_delta} kcal/mol**; screen predictions are the values evaluated.
The difference is expected from the locked practical-screening path: `fit` consumes
Boltzmann-averaged conformer Sterimol descriptors in the sigpack, whereas `screen`
recomputes Sterimol from each candidate's representative `Ligand_XYZ_Path`. It is a
real limitation of this experiment, not silently reconciled after target reveal.

The response is published `ddG_abs` at 298.15 K and larger is treated as more useful
enantioselectivity. Each candidate panel contains exactly three fully held-out
ligands. Native StericX descriptors and its constrained OLS/BIC training workflow
were used without post-result model changes. The exact random comparator enumerates
all six possible rankings for every panel; `random_baseline_metrics.csv` keeps every
one.

## Post-lock implementation QA

After the first outcome reveal, quality assurance changed only plot typography,
removed volatile timing text from a model audit file, and expanded the explanation
of the already-recorded screen-versus-fit descriptor mismatch. The original design
lock is preserved. `implementation_amendments.json` records both driver hashes and
the hashes proving that predictions, rankings, metrics, and baselines did not change.

## Failed and unfavorable results

{failure_text}

Poor correlations, bottom-ranked best ligands, negative selection gains, out-of-
domain predictions, and large errors are retained row-by-row in the machine-readable
artifacts. Primary metrics include failures rather than silently restricting the
denominator to successful screens.

## Reproducibility and artifacts

Run from the repository root:

```bash
cargo build --release
uv run --extra science python studies/study_011_ligand_ranking.py
uv run --extra science python studies/study_011_ligand_ranking.py --verify-only
```

- `DESIGN.md`, `design.json`, `design_lock.json`: outcome-independent protocol and lock
- `split_definitions.json`, `splits.csv`, `seeds.json`: exact partitions and seeds
- `frozen_predictions.csv`: target-free screen output written before reveal
- `rankings.csv`: predictions joined to revealed outcomes and rank annotations
- `split_metrics.csv`: every planned fold, including explicit failed-fold penalties
- `random_baseline_metrics.csv`: all 342 exact random rankings
- `model_reports.jsonl`, `screen_reports.jsonl`: per-fold model and screen audit trails
- `split_status.csv`: execution success/failure without suppressed folds
- `metrics.json`, `run_manifest.json`: aggregate results and provenance
- `implementation_amendments.json`: post-lock non-analytic QA disclosure
- `ranking_performance.png`, `applicability_domain.png`: generated figures

The source dataset is the authors' public [Ni-Catalyzed-hDA repository]
(https://github.com/SigmanGroup/Ni-Catalyzed-hDA), also linked by the
[CaltechAUTHORS article record](https://authors.library.caltech.edu/records/ef9d0-y2c15).
StericX pins the exact cached CSV by SHA-256 in `data/official/provenance.json`.

## Limitations

- Only 11 ligands have the required measured outcome; eight train and three rank in
  each fold. Spearman values for a three-item panel are coarse and unstable.
- The 57 folds overlap heavily and cannot be treated as 57 independent experiments.
- Candidate panels are deliberately small. This does not show recovery from all
  1,566 Kraken ligands, whose outcomes are mostly unknown.
- `ddG_abs` captures magnitude of enantioselectivity, not yield, cost, availability,
  stability, or synthetic practicality; “useful” is limited to that response.
- This is retrospective validation. Even strict computational hiding cannot replace
  a prospectively measured, untouched ligand screen.
- Scaffold labels come from the existing StericX preparation pipeline; alternate
  structure standardization could change grouping.
- The model is trained on ensemble-averaged Sterimol features but the current
  screening CSV path featurizes one representative conformer. The resulting
  screen-versus-fit prediction differences confound model-ranking quality with the
  descriptor aggregation mismatch; both predictions and their difference are saved.
- Applicability strata contain repeated predictions and sometimes very small counts;
  their error differences are descriptive only.
- The upstream data are publicly distributed by the authors, but this study does not
  make a new claim about upstream licensing or confer redistribution rights.
"""
    atomic_write_text(output_dir / "STUDY_011.md", report)


def result_manifest(
    binary: Path,
    output_dir: Path,
    design_lock_digest: str,
    frozen_digest: str,
) -> dict[str, Any]:
    version = run_command((str(binary), "--version")).stdout.strip()
    amendment_path = output_dir / "implementation_amendments.json"
    return {
        "study": "study_011",
        "driver": str(Path(__file__).resolve().relative_to(ROOT)),
        "single_command": (
            "uv run --extra science python studies/study_011_ligand_ranking.py"
        ),
        "stericx_binary": str(binary.resolve()),
        "stericx_version": version,
        "python": platform.python_version(),
        "numpy": np.__version__,
        "scipy": scipy.__version__,
        "matplotlib": matplotlib.__version__,
        "design_lock_sha256": design_lock_digest,
        "frozen_predictions_sha256": frozen_digest,
        "inputs": {
            str(path.relative_to(ROOT)): sha256_file(path)
            for path in (METADATA, CATALOG, CATALOG_PROVENANCE, SIGPACK, binary)
        },
        "artifacts": {
            name: sha256_file(output_dir / name)
            for name in RESULT_FILES
            if name not in {"run_manifest.json", "STUDY_011.md"}
            and (output_dir / name).is_file()
        },
        "post_lock_implementation_amendment": (
            {
                "artifact": amendment_path.name,
                "sha256": sha256_file(amendment_path),
            }
            if amendment_path.is_file()
            else None
        ),
        "note": (
            "No timestamps are stored so identical inputs can reproduce stable "
            "artifacts."
        ),
    }


def write_results(
    output_dir: Path,
    rankings: Sequence[dict[str, Any]],
    split_metrics: Sequence[dict[str, Any]],
    random_rows: Sequence[dict[str, Any]],
    metrics: dict[str, Any],
    binary: Path,
    frozen_digest: str,
) -> None:
    atomic_write_csv(output_dir / "rankings.csv", RANKING_FIELDS, rankings)
    atomic_write_csv(
        output_dir / "split_metrics.csv", SPLIT_METRIC_FIELDS, split_metrics
    )
    atomic_write_csv(
        output_dir / "random_baseline_metrics.csv", RANDOM_FIELDS, random_rows
    )
    design_lock_digest = sha256_file(output_dir / "design_lock.json")
    metrics = {
        **metrics,
        "design_lock_sha256": design_lock_digest,
        "frozen_predictions_sha256": frozen_digest,
    }
    atomic_write_json(output_dir / "metrics.json", metrics)
    plot_ranking_performance(
        output_dir / "ranking_performance.png", metrics, split_metrics, rankings
    )
    plot_applicability(output_dir / "applicability_domain.png", rankings)
    write_study_report(output_dir, metrics, frozen_digest, design_lock_digest)
    manifest = result_manifest(binary, output_dir, design_lock_digest, frozen_digest)
    atomic_write_json(output_dir / "run_manifest.json", manifest)


def verify_results(output_dir: Path) -> None:
    verify_design_lock(output_dir)
    missing = [name for name in RESULT_FILES if not (output_dir / name).is_file()]
    if missing:
        raise FileNotFoundError(f"missing result artifacts: {', '.join(missing)}")
    frozen_line = (
        (output_dir / "frozen_predictions.sha256")
        .read_text(encoding="utf-8")
        .split()[0]
    )
    actual_frozen = sha256_file(output_dir / "frozen_predictions.csv")
    if frozen_line != actual_frozen:
        raise RuntimeError(
            f"frozen prediction digest mismatch: {frozen_line} != {actual_frozen}"
        )
    manifest = json.loads(
        (output_dir / "run_manifest.json").read_text(encoding="utf-8")
    )
    errors = []
    for name, expected in manifest["artifacts"].items():
        actual = sha256_file(output_dir / name)
        if actual != expected:
            errors.append(f"{name}: expected {expected}, found {actual}")
    metrics = json.loads((output_dir / "metrics.json").read_text(encoding="utf-8"))
    split_rows = read_csv(output_dir / "split_metrics.csv")
    random_rows = read_csv(output_dir / "random_baseline_metrics.csv")
    checks = (
        (
            "primary top1",
            mean_defined(split_rows, "top1_recovery"),
            metrics["primary_intention_to_screen"]["top1_recovery"],
        ),
        (
            "primary top2",
            mean_defined(split_rows, "top2_recall"),
            metrics["primary_intention_to_screen"]["top2_recall"],
        ),
        (
            "random top1",
            mean_defined(random_rows, "top1_recovery"),
            metrics["exact_random_baseline"]["top1_recovery"],
        ),
        (
            "random top2",
            mean_defined(random_rows, "top2_recall"),
            metrics["exact_random_baseline"]["top2_recall"],
        ),
    )
    for name, reconstructed, stored in checks:
        if reconstructed is None or not math.isclose(
            reconstructed, float(stored), rel_tol=1e-12, abs_tol=1e-12
        ):
            errors.append(f"{name}: reconstructed {reconstructed}, stored {stored}")
    if len(split_rows) != 57:
        errors.append(f"split_metrics.csv has {len(split_rows)} rows, expected 57")
    if len(random_rows) != 57 * 6:
        errors.append(
            f"random_baseline_metrics.csv has {len(random_rows)} rows, expected 342"
        )
    if errors:
        raise RuntimeError("result verification failed:\n  " + "\n  ".join(errors))
    print("Study 011 verification passed")
    print(f"  design lock: {sha256_file(output_dir / 'design_lock.json')}")
    print(f"  frozen predictions: {actual_frozen}")
    print(f"  planned splits: {len(split_rows)}")


def main(argv: Iterable[str] | None = None) -> int:
    args = parse_args(argv)
    output_dir = args.output_dir.resolve()
    if args.prepare_only:
        lock = prepare_design(output_dir)
        print(f"Prepared and locked Study 011 design in {output_dir}")
        print(f"  design files: {len(lock['files'])}")
        print(f"  lock sha256: {sha256_file(output_dir / 'design_lock.json')}")
        return 0
    if args.verify_only:
        verify_results(output_dir)
        return 0
    if not output_dir.joinpath("design_lock.json").is_file():
        prepare_design(output_dir)
    lock = verify_design_lock(output_dir)
    binary = args.binary.resolve()
    if not binary.is_file():
        raise FileNotFoundError(
            f"StericX binary not found: {binary}; run `cargo build --release`"
        )
    rows = read_allowed_metadata()
    splits = json.loads(
        (output_dir / "split_definitions.json").read_text(encoding="utf-8")
    )
    frozen, models, screens, statuses = screen_all_splits(
        binary, output_dir, rows, splits
    )
    frozen_digest = write_frozen_artifacts(
        output_dir, frozen, models, screens, statuses
    )
    # Target revelation begins only after every planned screen has terminated and
    # the combined target-free prediction artifact has been hashed.
    outcomes = reveal_outcomes(rows)
    rankings, fold_metrics, random_rows, metrics = evaluate(
        frozen, statuses, splits, outcomes
    )
    write_results(
        output_dir,
        rankings,
        fold_metrics,
        random_rows,
        metrics,
        binary,
        frozen_digest,
    )
    # Re-verify the lock after analysis so the result records the same design that
    # authorized the run.
    if lock != verify_design_lock(output_dir):
        raise RuntimeError("design lock changed during study execution")
    verify_results(output_dir)
    success = sum(row["status"] == "success" for row in statuses)
    print(f"Study 011 complete: {success}/{len(splits)} splits ranked")
    print(f"  report: {output_dir / 'STUDY_011.md'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
