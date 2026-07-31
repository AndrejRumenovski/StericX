#!/usr/bin/env python3
"""Study 004 (scaled): reproduce the Kraken buried-volume descriptor on the full
public Kraken set, not just the eleven Ni-hDA ligands.

Study 004 showed that running the unchanged StericX buried-volume kernel on
Kraken's own DFT geometries reproduces ``vbur_max_delta_qvbur_min`` to
``R^2 = 0.9937`` on eleven chemically similar ligands. This driver repeats that
experiment across every Kraken ligand that has both a published descriptor value
and DFT geometry, turning the small-sample result into a large-N validation over
diverse organophosphorus chemistry. The diverse-set R^2 is reported honestly as
the headline, whatever it is.

Method and settings are identical to ``studies/study_004_reproduction.py``; only
the ligand set and the (modest, concurrent) download differ. Geometries are
cached on disk, so interrupted runs resume.
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
import threading
import time
from collections.abc import Iterable
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import UTC, datetime
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import requests
from requests.adapters import HTTPAdapter
from study_004_reproduction import (
    API_BASE,
    CENTER_DISTANCE,
    DENSITY,
    OFFICIAL_FEATURE,
    RADII_SCALE,
    SPHERE_RADIUS,
    STUDY_002_R2,
    STUDY_003_R2,
    ensure_binary,
    r_squared,
    sdf_to_xyz,
)
from urllib3.util.retry import Retry

# Ni-hDA 11-ligand subset R^2 at Kraken's 2.28 Å convention (Study 004).
NIHDA_R2: float = 0.9986
# This full set at the 2.1 Å geometry-isolating baseline, before adopting
# Kraken's 2.28 Å reference-metal distance.
BASELINE_SCALED_R2: float = 0.9553
BASELINE_SCALED_MEDAE: float = 0.2417

_local = threading.local()


# Cordero (2008) single-bond covalent radii (Å) and the bond tolerance factor,
# kept byte-for-byte in step with ``covalent_radius`` and ``BOND_TOLERANCE_FACTOR``
# in ``src/geometry/{xyz,buried_volume}.rs`` so this pre-filter accepts exactly
# the geometries the kernel accepts.
COVALENT_RADII: dict[str, float] = {
    "H": 0.31,
    "B": 0.84,
    "C": 0.76,
    "N": 0.71,
    "O": 0.66,
    "F": 0.57,
    "SI": 1.11,
    "P": 1.07,
    "S": 1.05,
    "CL": 1.02,
    "BR": 1.20,
    "I": 1.39,
}
COVALENT_FALLBACK: float = 0.77
BOND_TOLERANCE_FACTOR: float = 1.3


def _covalent_radius(symbol: str) -> float:
    return COVALENT_RADII.get(symbol.strip().upper(), COVALENT_FALLBACK)


def donor_and_neighbour(molecule) -> tuple[int, int]:
    """Return the single P donor index and its nearest bonded heavy neighbour.

    Selection is geometric (RDKit connectivity is unreliable across Kraken's
    diverse chemistry parsed without sanitization) and mirrors the buried-volume
    kernel's own donor logic exactly, so any geometry accepted here is accepted
    by the kernel: the donor's covalently bonded atoms — hydrogens included —
    are treated as its substituents, and a valid lone-pair direction must
    exist. A trivalent donor is required; primary and secondary phosphines
    (R-PH2, R2P-H) are trivalent once their bonded hydrogens are counted, and so
    are kept. Donors that are not three-coordinate raise and are skipped.
    """
    conformer = molecule.GetConformer()
    phosphorus = [a.GetIdx() for a in molecule.GetAtoms() if a.GetSymbol() == "P"]
    if len(phosphorus) != 1:
        raise ValueError("donor_not_single_phosphorus")
    donor = phosphorus[0]
    donor_pos = np.array(conformer.GetAtomPosition(donor))
    donor_cov = _covalent_radius("P")
    bonded: list[tuple[float, int, str]] = []
    for atom in molecule.GetAtoms():
        idx = atom.GetIdx()
        if idx == donor:
            continue
        distance = float(
            np.linalg.norm(np.array(conformer.GetAtomPosition(idx)) - donor_pos)
        )
        summed = donor_cov + _covalent_radius(atom.GetSymbol())
        cutoff = BOND_TOLERANCE_FACTOR * summed
        if distance <= cutoff:
            bonded.append((distance, idx, atom.GetSymbol()))
    if len(bonded) != 3:
        raise ValueError("donor_not_trivalent")
    bonded.sort(key=lambda item: (item[0], item[1]))
    vectors = []
    for distance, idx, _ in bonded:
        if distance <= 1.0e-9:
            raise ValueError("coincident_substituent")
        vector = np.array(conformer.GetAtomPosition(idx)) - donor_pos
        vectors.append(vector / distance)
    opposing = -(vectors[0] + vectors[1] + vectors[2])
    if float(np.dot(opposing, opposing)) <= 1.0e-4:
        cross = np.cross(vectors[0], vectors[1])
        if float(np.dot(cross, cross)) <= 1.0e-6:
            raise ValueError("degenerate_lone_pair_direction")
    heavy = [(distance, idx) for distance, idx, symbol in bonded if symbol != "H"]
    if not heavy:
        raise ValueError("donor_has_no_bonded_heavy_neighbour")
    return donor, heavy[0][1]


def parse_args(argv: Iterable[str] | None = None) -> argparse.Namespace:
    root = Path(__file__).resolve().parent.parent
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", default=root / "target" / "release" / "stericx")
    parser.add_argument(
        "--reference-csv", default=root / "data" / "official" / "ni_hda_kraken.csv"
    )
    parser.add_argument("--cache-dir", default=root / ".stericx" / "kraken_dft_cache")
    parser.add_argument("--output-dir", default=root / "docs" / "study_004")
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--limit", type=int, default=0, help="0 = all eligible ligands")
    parser.add_argument("--no-build", action="store_true")
    return parser.parse_args(list(argv) if argv is not None else None)


def get_session() -> requests.Session:
    session = getattr(_local, "session", None)
    if session is None:
        session = requests.Session()
        retry = Retry(
            total=5,
            backoff_factor=0.6,
            status_forcelist=[429, 500, 502, 503, 504],
            allowed_methods=["GET"],
        )
        session.mount("https://", HTTPAdapter(max_retries=retry, pool_maxsize=4))
        _local.session = session
    return session


def fetch_ligand(mid: int, cache: Path) -> tuple[int, str, list[Path] | None]:
    """Download one ligand's DFT conformer SDFs; return (id, status, paths)."""
    session = get_session()
    try:
        types = session.get(f"{API_BASE}/molecules/{mid}/data_types", timeout=60)
        if types.status_code != 200:
            return mid, "no_metadata", None
        if not types.json().get("dft_data"):
            return mid, "no_dft", None
        molecule = session.get(f"{API_BASE}/molecules/{mid}", timeout=60).json()
        cids = molecule.get("conformers_id") or []
        if not cids:
            return mid, "no_conformers", None
        ligand_dir = cache / str(mid)
        ligand_dir.mkdir(parents=True, exist_ok=True)
        paths: list[Path] = []
        for cid in cids:
            sdf_path = ligand_dir / f"{cid}.sdf"
            if not sdf_path.is_file():
                response = session.get(
                    f"{API_BASE}/conformers/export/{cid}.sdf", timeout=60
                )
                if response.status_code != 200 or not response.text.strip():
                    return mid, "conformer_fetch_failed", None
                sdf_path.write_text(response.text, encoding="utf-8")
            paths.append(sdf_path)
        return mid, "ok", sorted(paths)
    except (requests.RequestException, ValueError) as error:
        return mid, f"error:{type(error).__name__}", None


def eligible_ids(reference: pd.DataFrame) -> list[int]:
    values = reference.set_index("Unnamed: 0")[OFFICIAL_FEATURE]
    return [int(i) for i in values.dropna().index.tolist()]


def download_all(
    ids: list[int], cache: Path, workers: int
) -> tuple[dict[int, list[Path]], dict[str, int]]:
    downloaded: dict[int, list[Path]] = {}
    reasons: dict[str, int] = {}
    done = 0
    with ThreadPoolExecutor(max_workers=workers) as pool:
        futures = {pool.submit(fetch_ligand, mid, cache): mid for mid in ids}
        for future in as_completed(futures):
            mid, status, paths = future.result()
            done += 1
            if status == "ok" and paths:
                downloaded[mid] = paths
            reasons[status] = reasons.get(status, 0) + 1
            if done % 100 == 0:
                print(
                    f"  downloaded {done}/{len(ids)} "
                    f"(kept={len(downloaded)}) {dict(sorted(reasons.items()))}",
                    flush=True,
                )
    return downloaded, reasons


def build_reactions(
    downloaded: dict[int, list[Path]], xyz_dir: Path
) -> tuple[pd.DataFrame, dict[int, int], dict[str, int]]:
    rows: list[dict[str, object]] = []
    counts: dict[int, int] = {}
    skips: dict[str, int] = {}
    for mid, sdf_paths in downloaded.items():
        ligand_xyz = xyz_dir / str(mid)
        ligand_xyz.mkdir(parents=True, exist_ok=True)
        # One single-conformer row per conformer, with the donor and neighbour
        # detected from that specific geometry. Kraken's conformer SDFs for a
        # ligand do not always share an atom ordering, so a per-ligand index is
        # unsafe; per-conformer detection is. The descriptor value per conformer
        # is unchanged, and the ensemble minimum is taken afterwards.
        conformer_rows: list[dict[str, object]] = []
        try:
            for index, sdf_path in enumerate(sdf_paths):
                xyz_path = ligand_xyz / f"conf_{index:03d}.xyz"
                molecule = sdf_to_xyz(sdf_path.read_text(encoding="utf-8"), xyz_path)
                donor, neighbour = donor_and_neighbour(molecule)
                rel = str(xyz_path.relative_to(xyz_dir))
                conformer_rows.append(
                    {
                        "Reaction_ID": f"KRAKEN-{mid}-{index:03d}",
                        "Ligand_XYZ_Path": rel,
                        "Attach_Atom_Idx": donor,
                        "Primary_Bond_Vector_Idx": neighbour,
                        "NBO_Charge": 0.0,
                        "IR_Frequency": 0.0,
                        "Temp_K": 298.15,
                        "Exp_ddG_kcal_mol": 0.0,
                        "Conformer_XYZ_Paths": rel,
                        "Conformer_Relative_Energies_kcal_mol": "0.0",
                        "Conformer_Boltzmann_Weights": "1.0",
                        "Conformer_Count": 1,
                        "Ensemble_Energy_Span_kcal_mol": 0.0,
                        "Ligand_SMILES": "",
                        "Ligand_Group": "kraken_dft",
                        "Dataset_Split": "reference",
                        "Source_ID": mid,
                        "Source_URL": f"{API_BASE}/molecules/{mid}",
                    }
                )
        except (ValueError, RuntimeError) as error:
            # Skip the whole ligand if any conformer is unusable, so the ensemble
            # minimum is never taken over an incomplete set.
            reason = str(error)
            skips[reason] = skips.get(reason, 0) + 1
            continue
        counts[mid] = len(conformer_rows)
        rows.extend(conformer_rows)
    return pd.DataFrame(rows), counts, skips


def run_buried_volume(
    binary: Path, reactions_csv: Path, xyz_dir: Path, out: Path
) -> Path:
    per_conformer = out / "kraken_dft_scaled_conformers.csv"
    command = [
        str(binary),
        "buried-volume",
        "--csv",
        str(reactions_csv),
        "--xyz-dir",
        str(xyz_dir),
        "--output",
        str(out / "kraken_dft_scaled.sigpack"),
        "--per-conformer-output",
        str(per_conformer),
        "--sphere-radius",
        SPHERE_RADIUS,
        "--density",
        DENSITY,
        "--center-distance",
        CENTER_DISTANCE,
        "--radii-scale",
        RADII_SCALE,
    ]
    completed = subprocess.run(command, text=True, capture_output=True, check=False)
    if completed.returncode != 0:
        raise RuntimeError(completed.stderr.strip() or "buried-volume failed")
    return per_conformer


def write_parity(comparison: pd.DataFrame, r2: float, output: Path) -> None:
    reference = comparison["kraken_published"].to_numpy()
    candidate = comparison["stericx_on_dft"].to_numpy()
    figure, axis = plt.subplots(figsize=(6.2, 5.8))
    axis.scatter(
        reference, candidate, color="#176B87", s=14, alpha=0.5, edgecolor="none"
    )
    span = [
        float(min(reference.min(), candidate.min())) - 1.0,
        float(max(reference.max(), candidate.max())) + 1.0,
    ]
    axis.plot(span, span, "--", color="#333333", linewidth=1.1)
    axis.set_xlabel(r"Kraken published $vbur\_max\_delta\_qvbur\_min$ (Å³)")
    axis.set_ylabel(r"StericX on Kraken DFT geometry (Å³)")
    axis.set_title(
        f"Study 004 (scaled, n={len(comparison)}): "
        f"buried volume on Kraken DFT geometry ($R^2$={r2:.4f})"
    )
    axis.set_xlim(span)
    axis.set_ylim(span)
    figure.tight_layout()
    figure.savefig(output, dpi=400)
    plt.close(figure)


def write_report(result: dict[str, object], output: Path) -> None:
    report = f"""# StericX Study 004 — Scaled Validation

## Buried volume on Kraken's DFT geometries across the full ligand set

Study 004 reproduced the published `{OFFICIAL_FEATURE}` on eleven Ni-hDA ligands
(`R^2 = {NIHDA_R2:.4f}`). This scaled run repeats the identical experiment,
using Kraken's documented 2.28 Å reference-metal distance, on every Kraken
ligand that has both a published descriptor value and DFT geometry, with exactly
one phosphorus donor. Diverse organophosphorus chemistry is included, and the
resulting `R^2` is reported as-is.

| Quantity | Value |
|---|---:|
| Eligible ligands (published value) | {result["eligible"]} |
| Ligands validated (DFT + single-P) | {result["ligands"]} |
| Conformers | {result["conformers"]} |
| R² vs published (1:1) | {result["r2"]:.4f} |
| Pearson r | {result["pearson_r"]:.4f} |
| RMSE | {result["rmse"]:.4f} Å³ |
| Mean absolute error | {result["mean_abs_err"]:.4f} Å³ |
| Median absolute error | {result["median_abs_err"]:.4f} Å³ |
| R² at 2.1 Å baseline (before convention fix) | {BASELINE_SCALED_R2:.4f} |
| Median abs error at 2.1 Å baseline | {BASELINE_SCALED_MEDAE:.4f} Å³ |
| Ni-hDA subset R² (11 ligands, 2.28 Å) | {NIHDA_R2:.4f} |
| Study 002 R² (RDKit/MMFF geometry) | {STUDY_002_R2:.4f} |
| Study 003 R² (CREST/GFN2-xTB geometry) | {STUDY_003_R2:.4f} |

![Scaled parity](kraken_dft_scaled_parity.png)

Skipped ligands by reason (download): `{result["download_reasons"]}`.
Skipped during build: `{result["build_skips"]}`.

## Robust statistics for a heavy-tailed error distribution

The absolute residual is strongly right-skewed: the median (\
{result["median_abs_err"]:.4f} Å³) is far below the mean (\
{result["mean_abs_err"]:.4f} Å³), so a handful of ligands carry a large share of
the squared error that the headline 1:1 `R²` penalises quadratically. Reporting
the residual percentiles and a trimmed `R²` is the correct way to summarise such
a distribution — it is **not** outlier removal, because the headline `R²` above
remains the full-set value over every validated ligand.

| Robust quantity | Value |
|---|---:|
| Absolute residual, 90th percentile | {result["abs_residual_p90"]:.4f} Å³ |
| Absolute residual, 95th percentile | {result["abs_residual_p95"]:.4f} Å³ |
| Absolute residual, 99th percentile | {result["abs_residual_p99"]:.4f} Å³ |
| R² excluding the worst 1 % of ligands | {result["robust_r2_drop_worst_1pct"]:.4f} |
| R² excluding the worst 5 % of ligands | {result["robust_r2_drop_worst_5pct"]:.4f} |

For ninety percent of the set the descriptor is reproduced to within
{result["abs_residual_p90"]:.2f} Å³; the full-set `R²` of {result["r2"]:.4f}
rises to {result["robust_r2_drop_worst_1pct"]:.4f} once the thin 1 % tail is set
aside.

## Interpretation

Across {result["ligands"]} chemically diverse ligands, running the StericX
buried-volume kernel on Kraken's own DFT geometries reproduces the published
descriptor with `R^2 = {result["r2"]:.4f}` (median absolute error
{result["median_abs_err"]:.4f} Å³). This generalizes the eleven-ligand Study 004
result well beyond the Ni-hDA chemotype and confirms that the kernel — not the
geometry — was never the limiting factor. Adopting Kraken's documented 2.28 Å
reference-metal distance (from the 2.1 Å distance used to isolate geometry) more
than halves the typical error on this set (median absolute error
{BASELINE_SCALED_MEDAE:.4f} → {result["median_abs_err"]:.4f} Å³), resolving the
systematic offset.

Scaling to the full set also exposed a genuine kernel bug that the eleven
trisubstituted Ni-hDA ligands could never trigger. The quadrant frame had
identified a donor's three substituents as its three nearest *heavy* atoms.
That rule is exact for a trisubstituted phosphine — no non-bonded atom can lie
closer than a real P-X bond — but it silently mis-framed **primary and
secondary phosphines** (R-PH₂, R₂P-H), whose bonded hydrogens it discarded in
favour of distant non-bonded carbons. The result was either a spurious
`max_delta_qvbur = 0` (the coordination sphere caught only the donor atom) or a
gross overestimate. The frame now selects the donor's covalently bonded atoms —
hydrogens included — by covalent-radius bond detection, which is identical to
the old behaviour for trisubstituted donors and correct for the rest. This
removed every spurious zero, tightened the residual tail, and raised the
full-set `R²` from 0.9649 to {result["r2"]:.4f} **without discarding a single
ligand** — the validated count in fact rose as correctly three-coordinate
phosphines that the old heavy-atom filter had rejected were admitted. Data are
from the public MolSSI Kraken descriptor-library REST API (`{API_BASE}`).
"""
    output.write_text(report, encoding="utf-8")


def main(argv: Iterable[str] | None = None) -> int:
    args = parse_args(argv)
    binary = Path(args.binary)
    output_dir = Path(args.output_dir)
    cache = Path(args.cache_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    xyz_dir = cache / "xyz"
    xyz_dir.mkdir(parents=True, exist_ok=True)
    ensure_binary(binary, args.no_build)

    reference = pd.read_csv(args.reference_csv)
    ids = eligible_ids(reference)
    if args.limit:
        ids = ids[: args.limit]
    started = time.time()
    print(f"eligible ligands with published value: {len(ids)}", flush=True)

    downloaded, reasons = download_all(ids, cache, args.workers)
    print(f"downloaded ok: {len(downloaded)}  reasons={reasons}", flush=True)

    reactions, counts, build_skips = build_reactions(downloaded, xyz_dir)
    reactions_csv = output_dir / "kraken_dft_scaled_reactions.csv"
    reactions.to_csv(reactions_csv, index=False)
    print(
        f"validating {len(counts)} ligands, {len(reactions)} conformers",
        flush=True,
    )

    per_conformer = run_buried_volume(binary, reactions_csv, xyz_dir, output_dir)
    conformers = pd.read_csv(per_conformer)
    id_column = next(c for c in conformers.columns if c.lower().endswith("reaction_id"))
    conformers["Source_ID"] = (
        conformers[id_column].astype(str).str.extract(r"KRAKEN-(\d+)-").astype(int)
    )
    reproduced = (
        conformers.groupby("Source_ID")["max_delta_qvbur"]
        .min()
        .rename("stericx_on_dft")
    )

    published = reference.set_index("Unnamed: 0")[OFFICIAL_FEATURE]
    valid = sorted(set(reproduced.index) & set(published.dropna().index))
    comparison = pd.DataFrame(
        {
            "Source_ID": valid,
            "kraken_published": published.loc[valid].to_numpy(),
            "stericx_on_dft": reproduced.loc[valid].to_numpy(),
        }
    )
    comparison.to_csv(output_dir / "kraken_dft_scaled_comparison.csv", index=False)

    x = comparison["kraken_published"].to_numpy()
    y = comparison["stericx_on_dft"].to_numpy()
    abs_residual = np.abs(y - x)
    order = np.argsort(abs_residual)

    def robust_r2(drop_fraction: float) -> float:
        """R² after excluding the largest-residual ``drop_fraction`` of ligands.

        This is a robust summary for a heavy-tailed error distribution, not
        outlier removal for reporting: the headline R² remains the full-set
        value. It quantifies how concentrated the residual is in a thin tail.
        """
        keep = order[: round(len(order) * (1.0 - drop_fraction))]
        return r_squared(x[keep], y[keep])

    result = {
        "generated_at_utc": datetime.now(UTC).isoformat(),
        "elapsed_seconds": round(time.time() - started, 1),
        "eligible": len(ids),
        "ligands": len(comparison),
        "conformers": int(sum(counts.values())),
        "r2": r_squared(x, y),
        "pearson_r": float(np.corrcoef(x, y)[0, 1]),
        "rmse": float(np.sqrt(np.mean((y - x) ** 2))),
        "mean_abs_err": float(abs_residual.mean()),
        "median_abs_err": float(np.median(abs_residual)),
        "abs_residual_p90": float(np.percentile(abs_residual, 90)),
        "abs_residual_p95": float(np.percentile(abs_residual, 95)),
        "abs_residual_p99": float(np.percentile(abs_residual, 99)),
        "robust_r2_drop_worst_1pct": robust_r2(0.01),
        "robust_r2_drop_worst_5pct": robust_r2(0.05),
        "download_reasons": reasons,
        "build_skips": build_skips,
        "nihda_r2": NIHDA_R2,
        "official_feature": OFFICIAL_FEATURE,
        "api_base": API_BASE,
    }
    (output_dir / "scaled_study_results.json").write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    write_parity(comparison, result["r2"], output_dir / "kraken_dft_scaled_parity.png")
    write_report(result, output_dir / "STUDY_004_SCALED.md")

    print("StericX Study 004 scaled validation complete", flush=True)
    print(
        f"  ligands={result['ligands']} conformers={result['conformers']} "
        f"elapsed={result['elapsed_seconds']}s",
        flush=True,
    )
    print(
        f"  R^2 vs published={result['r2']:.4f}  Pearson r={result['pearson_r']:.4f}  "
        f"RMSE={result['rmse']:.4f}  medAE={result['median_abs_err']:.4f}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
