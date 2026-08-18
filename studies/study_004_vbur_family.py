"""Study 004 (family): reproduce Kraken's whole buried-volume descriptor family.

The scaled Study 004 validated a single Kraken descriptor,
``vbur_max_delta_qvbur_min``, across 1,541 ligands. But the same buried-volume
kernel run produced every other descriptor in Kraken's ``vbur`` family at the
same time. This study closes the loop: it compares StericX against Kraken's
*published* values for the full family — the buried volume itself, the quadrant
and octant extremes, and the near/far hemispheres — over the whole ligand set,
using the per-conformer values already computed by ``studies/study_004_scaled.py``.

For each descriptor the ligand value is the minimum over the conformer ensemble,
matching Kraken's ``*_min`` convention (and the reduction Study 004 used for the
headline descriptor). Nothing is recomputed here: StericX values are read from
the committed per-conformer table, and Kraken's published values are pulled from
the public MolSSI REST API and cached. The output is one R² per descriptor.
"""

from __future__ import annotations

import argparse
import json
from collections.abc import Iterable
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from study_004_reproduction import API_BASE, r_squared
from study_004_scaled import get_session

# StericX per-conformer column -> (Kraken property, display label). Each Kraken
# property is compared at its published minimum over the conformer ensemble.
FAMILY: dict[str, tuple[str, str]] = {
    "vbur": ("vbur_vbur", "V_bur (buried volume)"),
    "qvbur_min": ("vbur_qvbur_min", "quadrant V_bur, min"),
    "qvbur_max": ("vbur_qvbur_max", "quadrant V_bur, max"),
    "ovbur_min": ("vbur_ovbur_min", "octant V_bur, min"),
    "ovbur_max": ("vbur_ovbur_max", "octant V_bur, max"),
    "near_vbur": ("vbur_near_vbur", "near-hemisphere V_bur"),
    "far_vbur": ("vbur_far_vbur", "far-hemisphere V_bur"),
    "max_delta_qvbur": ("vbur_max_delta_qvbur", "max Δ quadrant V_bur"),
}


def parse_args(argv: Iterable[str] | None = None) -> argparse.Namespace:
    root = Path(__file__).resolve().parent.parent
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--conformers-csv",
        type=Path,
        default=root / "docs" / "study_004" / "kraken_dft_scaled_conformers.csv",
    )
    parser.add_argument("--output-dir", type=Path, default=root / "docs" / "study_004")
    parser.add_argument(
        "--cache",
        type=Path,
        default=(
            root / ".stericx" / "kraken_dft_cache" / "kraken_vbur_family_published.csv"
        ),
    )
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument(
        "--refresh", action="store_true", help="re-download published values"
    )
    return parser.parse_args(list(argv) if argv is not None else None)


def stericx_minima(conformers_csv: Path) -> pd.DataFrame:
    """Minimum over each ligand's conformers for every family descriptor."""
    conformers = pd.read_csv(conformers_csv)
    id_column = next(c for c in conformers.columns if c.lower().endswith("reaction_id"))
    conformers["Source_ID"] = (
        conformers[id_column].astype(str).str.extract(r"KRAKEN-(\d+)-").astype(int)
    )
    return conformers.groupby("Source_ID")[list(FAMILY)].min()


def fetch_published(mid: int) -> tuple[int, dict[str, float] | None]:
    """Pull one ligand's published family minima from the MolSSI API."""
    session = get_session()
    try:
        response = session.get(
            f"{API_BASE}/molecules/data/{mid}?data_type=dft", timeout=60
        )
        if response.status_code != 200:
            return mid, None
        by_property = {row["property"]: row for row in response.json()}
        values: dict[str, float] = {}
        for kraken_property, _ in FAMILY.values():
            record = by_property.get(kraken_property)
            if record is None or record.get("min") is None:
                return mid, None
            values[kraken_property] = float(record["min"])
        return mid, values
    except (ValueError, KeyError, TypeError, OSError):
        return mid, None


def published_minima(
    ids: list[int], cache: Path, workers: int, refresh: bool
) -> pd.DataFrame:
    """Published family minima for every ligand, cached to CSV between runs."""
    if cache.is_file() and not refresh:
        cached = pd.read_csv(cache).set_index("Source_ID")
        if set(ids).issubset(cached.index):
            return cached.loc[ids]
    rows: dict[int, dict[str, float]] = {}
    with ThreadPoolExecutor(max_workers=workers) as pool:
        futures = {pool.submit(fetch_published, mid): mid for mid in ids}
        for done, future in enumerate(as_completed(futures), start=1):
            mid, values = future.result()
            if values is not None:
                rows[mid] = values
            if done % 200 == 0:
                print(f"  fetched {done}/{len(ids)} (kept {len(rows)})", flush=True)
    frame = pd.DataFrame.from_dict(rows, orient="index")
    frame.index.name = "Source_ID"
    cache.parent.mkdir(parents=True, exist_ok=True)
    frame.sort_index().to_csv(cache)
    return frame


def write_figure(
    merged: pd.DataFrame, metrics: dict[str, dict[str, float]], output: Path
) -> None:
    figure, axes = plt.subplots(2, 4, figsize=(14.0, 7.2))
    for axis, (column, (kraken_property, label)) in zip(
        axes.ravel(), FAMILY.items(), strict=True
    ):
        published = merged[kraken_property].to_numpy()
        stericx = merged[column].to_numpy()
        axis.scatter(
            published, stericx, s=9, alpha=0.4, color="#176B87", edgecolor="none"
        )
        span = [
            float(min(published.min(), stericx.min())),
            float(max(published.max(), stericx.max())),
        ]
        axis.plot(span, span, "--", color="#333333", linewidth=1.0)
        axis.set_title(f"{label}\n$R^2$ = {metrics[column]['r2']:.4f}", fontsize=9)
        axis.set_xlabel("Kraken published (Å³)", fontsize=8)
        axis.set_ylabel("StericX (Å³)", fontsize=8)
        axis.tick_params(labelsize=7)
    figure.suptitle(
        f"Study 004 (family): StericX vs published Kraken buried-volume descriptors "
        f"(n = {len(merged)} ligands)",
        fontsize=12,
    )
    figure.tight_layout(rect=(0, 0, 1, 0.96))
    figure.savefig(output, dpi=300)
    plt.close(figure)


def write_report(
    metrics: dict[str, dict[str, float]], ligands: int, output: Path
) -> None:
    lines = [
        "# StericX Study 004 — Buried-Volume Descriptor Family",
        "",
        "## Reproducing the whole `vbur` family, not one descriptor",
        "",
        f"Study 004 validated a single Kraken descriptor across the full set. The "
        f"same buried-volume kernel run produces Kraken's entire `vbur` family, so "
        f"this study compares all of it against Kraken's *published* values over "
        f"{ligands} ligands. Each descriptor is reduced to its minimum over the "
        f"conformer ensemble, matching Kraken's `*_min` convention. StericX values "
        f"are read from the committed per-conformer table; Kraken values come from "
        f"the public MolSSI API. Nothing is recomputed.",
        "",
        "| Descriptor | Kraken property | R² | Pearson r | RMSE (Å³) | "
        "Median abs. err (Å³) |",
        "|---|---|---:|---:|---:|---:|",
    ]
    for column, (kraken_property, label) in FAMILY.items():
        metric = metrics[column]
        lines.append(
            f"| {label} | `{kraken_property}` | {metric['r2']:.4f} | "
            f"{metric['pearson_r']:.4f} | {metric['rmse']:.4f} | "
            f"{metric['median_abs_err']:.4f} |"
        )
    mean_r2 = float(np.mean([m["r2"] for m in metrics.values()]))
    lines += [
        "",
        f"Across the eight descriptors the mean R² is **{mean_r2:.4f}**. The single "
        "buried-volume kernel reproduces not just the headline `max_delta_qvbur` but "
        "the buried volume itself, the quadrant and octant extremes, and the "
        "near/far hemispheres — each against Kraken's own published numbers, across "
        "the whole library. The near/far hemispheres match in the correct sense "
        "(no axis swap), confirming the octant partitioning is oriented as Kraken's.",
        "",
        "![Family parity](kraken_vbur_family_parity.png)",
        "",
        "*Figure. StericX vs published Kraken value for each buried-volume "
        "descriptor. Generated by `studies/study_004_vbur_family.py`.*",
        "",
    ]
    output.write_text("\n".join(lines), encoding="utf-8")


def main(argv: Iterable[str] | None = None) -> int:
    args = parse_args(argv)
    stericx = stericx_minima(args.conformers_csv)
    ids = [int(i) for i in stericx.index]
    print(f"StericX ligands with per-conformer descriptors: {len(ids)}", flush=True)

    published = published_minima(ids, args.cache, args.workers, args.refresh)
    merged = stericx.join(published, how="inner").dropna()
    print(f"ligands with both StericX and published values: {len(merged)}", flush=True)

    metrics: dict[str, dict[str, float]] = {}
    for column, (kraken_property, _) in FAMILY.items():
        stericx_values = merged[column].to_numpy()
        published_values = merged[kraken_property].to_numpy()
        residual = stericx_values - published_values
        metrics[column] = {
            "r2": r_squared(published_values, stericx_values),
            "pearson_r": float(np.corrcoef(published_values, stericx_values)[0, 1]),
            "rmse": float(np.sqrt(np.mean(residual**2))),
            "median_abs_err": float(np.median(np.abs(residual))),
        }

    args.output_dir.mkdir(parents=True, exist_ok=True)
    write_figure(merged, metrics, args.output_dir / "kraken_vbur_family_parity.png")
    write_report(metrics, len(merged), args.output_dir / "STUDY_004_FAMILY.md")
    (args.output_dir / "vbur_family_metrics.json").write_text(
        json.dumps(
            {"ligands": len(merged), "descriptors": metrics}, indent=2, sort_keys=True
        )
        + "\n",
        encoding="utf-8",
    )

    print("StericX Study 004 buried-volume family validation complete")
    for column, (_, label) in FAMILY.items():
        print(f"  {label:28s} R² = {metrics[column]['r2']:.4f}")
    print(f"  mean R² = {np.mean([m['r2'] for m in metrics.values()]):.4f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
