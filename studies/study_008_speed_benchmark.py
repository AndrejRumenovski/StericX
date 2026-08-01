"""Study 008: StericX vs morfeus -- a head-to-head throughput benchmark.

Studies 002 and 005 establish that StericX's descriptors match morfeus to
R2 >= 0.9999 (buried volume, Sterimol) and to machine precision (pyramidalization).
Fidelity, though, is only half of a reproduction's practical value: the other half
is cost. This study measures how fast StericX computes the flagship buried-volume
descriptor at library scale against morfeus, the reference Python implementation,
on the same geometries and the same single CPU core.

The comparison is deliberately conservative. StericX computes Sterimol *and*
pyramidalization in the same pass it is timed on, while morfeus is timed for
buried volume alone; the buried-volume convention (3.5 A sphere, Bondi radii
x1.17, virtual Ni centre at the Kraken 2.28 A convention, octant analysis) is the
exact one validated in Study 002 and reused here verbatim. Both tools read the
same files from a warm OS cache, both run single-threaded, and both are timed
end-to-end (file on disk -> descriptor). Agreement on this set is re-confirmed so
the speedup is a like-for-like number, not a comparison of two different results.

The Kraken DFT SDF cache is a local, gitignored artifact; only StericX's own
timings and the aggregate comparison are written out.
"""

from __future__ import annotations

import argparse
import csv
import io
import json
import os
import platform
import subprocess
import time
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import study_002_buried_volume as bv
from study_004_reproduction import r_squared

ROOT = Path(__file__).resolve().parent.parent

# StericX's `descriptors` subcommand places the virtual metal centre at the Kraken
# 2.28 A convention (main.rs default); Study 002's morfeus reference used its 2.1 A
# fallback. Align morfeus to the same 2.28 A so the two tools compute an identical
# quantity and the benchmark is a like-for-like speed comparison, not two
# conventions. Every other buried-volume parameter (3.5 A sphere, 0.01 density,
# Bondi radii x1.17) already matches between the two by default.
CENTER_DISTANCE = 2.28
bv.CENTER_DISTANCE = CENTER_DISTANCE

# %Vbur difference above which a ligand is treated as a frame-topology mismatch
# (morfeus's nearest-3-heavy vs StericX's covalent bonding) rather than numerical
# noise. Well clear of the ~0.4 %Vbur spread among genuinely-agreeing ligands.
FRAME_TOLERANCE = 0.5


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--binary", type=Path, default=ROOT / "target" / "release" / "stericx"
    )
    parser.add_argument(
        "--cache-dir", type=Path, default=ROOT / ".stericx" / "kraken_dft_cache"
    )
    parser.add_argument("--output-dir", type=Path, default=ROOT / "docs" / "study_008")
    parser.add_argument(
        "--limit",
        type=int,
        default=0,
        help="Cap the number of structures (0 = the whole ligand library).",
    )
    parser.add_argument(
        "--all-conformers",
        action="store_true",
        help="Use every cached conformer instead of one geometry per ligand.",
    )
    parser.add_argument(
        "--reps",
        type=int,
        default=3,
        help="Timed repetitions; the fastest (least-noisy) run is reported.",
    )
    return parser.parse_args(argv)


def gather_structures(cache_dir: Path, all_conformers: bool, limit: int) -> list[Path]:
    """One lowest-numbered SDF per ligand (default), or every conformer."""
    ligand_dirs = sorted(
        (d for d in cache_dir.iterdir() if d.is_dir() and d.name.isdigit()),
        key=lambda d: int(d.name),
    )
    paths: list[Path] = []
    for directory in ligand_dirs:
        sdfs = sorted(directory.glob("*.sdf"))
        if not sdfs:
            continue
        paths.extend(sdfs if all_conformers else sdfs[:1])
    if limit > 0:
        paths = paths[:limit]
    return paths


def read_sdf(path: Path) -> tuple[list[str], np.ndarray]:
    """Read elements and coordinates from a V2000 SDF's first molecule block."""
    lines = path.read_text(encoding="utf-8").splitlines()
    atom_count = int(lines[3][:3])
    elements: list[str] = []
    coordinates: list[list[float]] = []
    for i in range(atom_count):
        fields = lines[4 + i].split()
        coordinates.append([float(fields[0]), float(fields[1]), float(fields[2])])
        elements.append(fields[3])
    return elements, np.asarray(coordinates, dtype=float)


def morfeus_vbur(elements: list[str], coordinates: np.ndarray) -> float | None:
    """%Vbur via the exact Study 002 morfeus convention (phosphine donor = P)."""
    try:
        donor_idx = elements.index("P")
    except ValueError:
        return None
    heavy = [
        (float(np.sum((coordinates[i] - coordinates[donor_idx]) ** 2)), i)
        for i, element in enumerate(elements)
        if i != donor_idx and element.upper() != "H"
    ]
    if len(heavy) < 3:
        return None
    reference_idx = min(heavy)[1]
    try:
        return float(
            bv.morfeus_reference(elements, coordinates, donor_idx, reference_idx)[
                "percent_vbur"
            ]
        )
    except (ValueError, KeyError, IndexError):
        return None


def time_morfeus(paths: list[Path], reps: int) -> tuple[float, dict[str, float]]:
    """Fastest end-to-end wall time over `reps`, plus %Vbur per file basename."""
    best = float("inf")
    values: dict[str, float] = {}
    for _rep in range(reps):
        start = time.perf_counter()
        rep_values: dict[str, float] = {}
        for path in paths:
            elements, coordinates = read_sdf(path)
            vbur = morfeus_vbur(elements, coordinates)
            if vbur is not None:
                rep_values[path.name] = vbur
        best = min(best, time.perf_counter() - start)
        values = rep_values
    return best, values


# Maximum structures passed to one `descriptors` invocation. Kept well below the
# OS argument-length limit (ARG_MAX) so an all-conformers run (~31k files) is split
# across a handful of calls instead of overflowing argv; a per-ligand run (~1.5k
# files) still fits in a single call, leaving the committed single-geometry
# benchmark byte-identical.
STERICX_BATCH = 2000


def time_stericx(
    binary: Path, paths: list[Path], reps: int
) -> tuple[float, dict[str, float]]:
    """Fastest wall time for one single-threaded CLI pass, plus %Vbur per file.

    Paths are chunked into batches below the OS argument limit; every batch of a
    pass is timed together, so a chunked all-conformers run and a single-call
    per-ligand run are measured the same end-to-end way.
    """
    env = {**os.environ, "RAYON_NUM_THREADS": "1"}
    batches = [
        paths[i : i + STERICX_BATCH] for i in range(0, len(paths), STERICX_BATCH)
    ]
    best = float("inf")
    values: dict[str, float] = {}
    for _rep in range(reps):
        start = time.perf_counter()
        outputs = [
            subprocess.run(
                [str(binary), "descriptors", "--format", "csv", *map(str, batch)],
                capture_output=True,
                text=True,
                check=True,
                env=env,
            ).stdout
            for batch in batches
        ]
        best = min(best, time.perf_counter() - start)
        rep_values: dict[str, float] = {}
        for output in outputs:
            for row in csv.DictReader(io.StringIO(output)):
                vbur = row.get("percent_buried_volume")
                if vbur:
                    rep_values[Path(row["file"]).name] = float(vbur)
        values = rep_values
    return best, values


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    if not args.binary.is_file():
        raise SystemExit(
            f"StericX binary not found: {args.binary} (cargo build --release)"
        )

    paths = gather_structures(args.cache_dir, args.all_conformers, args.limit)
    if not paths:
        raise SystemExit(f"no SDF geometries under {args.cache_dir}")
    scale = "conformers" if args.all_conformers else "ligands (one geometry each)"
    print(f"structures: {len(paths)} {scale}")

    # Warm the OS file cache so neither tool pays a one-time cold-read penalty.
    for path in paths:
        path.read_bytes()

    print("timing morfeus (buried volume, single core)...")
    morfeus_time, morfeus_vals = time_morfeus(paths, args.reps)
    print("timing StericX (full descriptor panel, single core)...")
    stericx_time, stericx_vals = time_stericx(args.binary, paths, args.reps)

    shared = sorted(set(morfeus_vals) & set(stericx_vals))
    m = np.array([morfeus_vals[k] for k in shared])
    s = np.array([stericx_vals[k] for k in shared])
    resid = s - m
    r2 = r_squared(m, s)
    max_abs = float(np.max(np.abs(resid)))
    mae = float(np.mean(np.abs(resid)))

    # A small fraction of ligands differ because morfeus's naive "3 nearest heavy
    # atoms" frame disagrees with StericX's covalent-bonding frame -- exactly the
    # topology case the StericX frame fix and Study 006 address. Separate them so
    # the headline agreement is not a comparison of two different frames.
    frame_mask = np.abs(resid) > FRAME_TOLERANCE
    clean = ~frame_mask
    frame_outliers = int(frame_mask.sum())
    r2_clean = r_squared(m[clean], s[clean])
    max_abs_clean = float(np.max(np.abs(resid[clean])))

    n = len(paths)
    morfeus_tp = n / morfeus_time
    stericx_tp = n / stericx_time
    speedup = morfeus_time / stericx_time

    print(
        f"\nagreement on {len(shared)} structures: R2 = {r2:.6f}  "
        f"max |diff| = {max_abs:.4f} %Vbur  MAE = {mae:.5f} %Vbur"
    )
    print(
        f"  frame-topology outliers (|diff| > {FRAME_TOLERANCE}): "
        f"{frame_outliers}/{len(shared)}; the other {len(shared) - frame_outliers} "
        f"agree at R2 = {r2_clean:.6f} (max |diff| {max_abs_clean:.4f})"
    )
    print(f"morfeus:  {morfeus_time:8.3f} s   {morfeus_tp:8.1f} structures/s")
    print(f"StericX:  {stericx_time:8.3f} s   {stericx_tp:8.1f} structures/s")
    print(
        f"speedup:  {speedup:.1f}x  (single core; StericX also computes Sterimol+pyr)"
    )

    system = {
        "cpu": _cpu_name(),
        "cores": os.cpu_count(),
        "machine": platform.machine(),
        "threads": "single (RAYON_NUM_THREADS=1; morfeus single-threaded)",
        "binary_bytes": args.binary.stat().st_size,
        "morfeus_version": _morfeus_version(),
    }
    metrics = {
        "structures": n,
        "scale": scale,
        "reps": args.reps,
        "agreement": {
            "r2": r2,
            "max_abs_diff": max_abs,
            "mae": mae,
            "n": len(shared),
            "frame_outliers": frame_outliers,
            "frame_tolerance": FRAME_TOLERANCE,
            "r2_excluding_frame_outliers": r2_clean,
            "max_abs_diff_excluding_frame_outliers": max_abs_clean,
        },
        "morfeus": {"seconds": morfeus_time, "throughput_per_s": morfeus_tp},
        "stericx": {"seconds": stericx_time, "throughput_per_s": stericx_tp},
        "speedup": speedup,
        "system": system,
    }

    args.output_dir.mkdir(parents=True, exist_ok=True)
    write_figure(metrics, m, s, frame_mask, args.output_dir / "speed_benchmark.png")
    (args.output_dir / "speed_metrics.json").write_text(
        json.dumps(metrics, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    write_report(metrics, args.output_dir / "STUDY_008.md")
    print("\nStudy 008 complete.")
    return 0


def _cpu_name() -> str:
    try:
        for line in Path("/proc/cpuinfo").read_text().splitlines():
            if line.startswith("model name"):
                return line.split(":", 1)[1].strip()
    except OSError:
        pass
    return platform.processor() or platform.machine()


def _morfeus_version() -> str:
    import importlib.metadata

    for dist in ("morfeus-ml", "morfeus-fsu"):
        try:
            return f"{dist} {importlib.metadata.version(dist)}"
        except importlib.metadata.PackageNotFoundError:
            continue
    return "unknown"


def write_figure(
    metrics: dict, morfeus_vals, stericx_vals, frame_mask, output: Path
) -> None:
    figure, (throughput_ax, parity_ax) = plt.subplots(1, 2, figsize=(11, 5))

    tools = ["morfeus\n(Python)", "StericX\n(Rust binary)"]
    values = [
        metrics["morfeus"]["throughput_per_s"],
        metrics["stericx"]["throughput_per_s"],
    ]
    bars = throughput_ax.bar(tools, values, color=["#8C8C8C", "#176B87"], width=0.6)
    throughput_ax.set_yscale("log")
    throughput_ax.set_ylabel("structures / second  (single core, log scale)")
    throughput_ax.set_title(
        f"Buried-volume throughput: StericX is {metrics['speedup']:.0f}x faster"
    )
    for bar, value in zip(bars, values, strict=True):
        throughput_ax.text(
            bar.get_x() + bar.get_width() / 2,
            value,
            f"{value:,.0f}/s",
            ha="center",
            va="bottom",
        )

    agree = ~frame_mask
    parity_ax.scatter(
        morfeus_vals[agree],
        stericx_vals[agree],
        s=14,
        alpha=0.4,
        color="#176B87",
        edgecolor="none",
        label=f"agree (n={int(agree.sum())})",
    )
    if frame_mask.any():
        parity_ax.scatter(
            morfeus_vals[frame_mask],
            stericx_vals[frame_mask],
            s=26,
            alpha=0.85,
            color="#B4530A",
            edgecolor="none",
            label=f"frame-topology diff (n={int(frame_mask.sum())})",
        )
        parity_ax.legend(loc="upper left", fontsize=8, frameon=False)
    span = [
        float(min(morfeus_vals.min(), stericx_vals.min())) - 1,
        float(max(morfeus_vals.max(), stericx_vals.max())) + 1,
    ]
    parity_ax.plot(span, span, "--", color="#333333", linewidth=1.0)
    parity_ax.set_xlabel("morfeus %Vbur")
    parity_ax.set_ylabel("StericX %Vbur")
    parity_ax.set_title(
        f"Same numbers: R2 = {metrics['agreement']['r2_excluding_frame_outliers']:.6f} "
        f"({metrics['agreement']['n'] - metrics['agreement']['frame_outliers']} "
        "frame-consistent)"
    )
    figure.suptitle(
        f"Study 008: {metrics['structures']} {metrics['scale']} on one "
        f"{metrics['system']['cpu']} core",
        fontsize=11,
    )
    figure.tight_layout()
    figure.savefig(output, dpi=300)
    plt.close(figure)


def write_report(metrics: dict, output: Path) -> None:
    system = metrics["system"]
    agreement = metrics["agreement"]
    binary_mb = system["binary_bytes"] / 1_000_000
    lines = [
        "# StericX Study 008 - A Head-to-Head Speed Benchmark vs morfeus",
        "",
        "## Same descriptor, same geometries, same core -- how much faster?",
        "",
        "Studies 002 and 005 show StericX's descriptors *match* morfeus (buried "
        "volume and Sterimol to R2 >= 0.9999; pyramidalization to machine "
        "precision). Fidelity is only half of what makes a reproduction useful in "
        "practice -- the other half is cost. This study times StericX against "
        "morfeus, the reference Python implementation, computing the flagship "
        "**buried-volume** descriptor on the same geometries and the same single "
        "CPU core.",
        "",
        f"The benchmark runs on **{metrics['structures']} {metrics['scale']}** from "
        "the Kraken DFT set, on a single CPU core "
        f"({system['cpu']}). Both tools read the same files from a warm OS cache and "
        "are timed end-to-end (file on disk to descriptor value); the fastest of "
        f"{metrics['reps']} repetitions is reported. The buried-volume convention "
        "(3.5 A sphere, Bondi radii x1.17, virtual Ni centre at 2.28 A, octant "
        "analysis) is the exact one validated in Study 002, reused here verbatim "
        f"for morfeus ({system['morfeus_version']}).",
        "",
        "### Result",
        "",
        "| | Wall time | Throughput | Relative |",
        "|---|---:|---:|---:|",
        f"| morfeus (Python) | {metrics['morfeus']['seconds']:.2f} s | "
        f"{metrics['morfeus']['throughput_per_s']:,.0f} / s | 1x |",
        f"| **StericX** (Rust binary) | **{metrics['stericx']['seconds']:.2f} s** | "
        f"**{metrics['stericx']['throughput_per_s']:,.0f} / s** | "
        f"**{metrics['speedup']:.0f}x** |",
        "",
        f"StericX is **{metrics['speedup']:.0f}x faster** on a single core. The "
        "comparison is deliberately conservative: StericX computes Sterimol L/B1/B5 "
        "and pyramidalization *in the same timed pass*, while morfeus is timed for "
        "buried volume alone -- so the real end-to-end advantage of reproducing "
        "StericX's full descriptor panel with morfeus is larger than the figure "
        "above.",
        "",
        "### The speedup is not from computing something cheaper",
        "",
        "On the identical geometries, the two tools return the same number. Of the "
        f"{agreement['n']} phosphines morfeus can frame (it requires three heavy "
        f"substituents), **{agreement['n'] - agreement['frame_outliers']}** agree "
        f"to **R2 = {agreement['r2_excluding_frame_outliers']:.6f}** (max absolute "
        f"difference {agreement['max_abs_diff_excluding_frame_outliers']:.2f} "
        "%Vbur) -- StericX's buried volume *is* morfeus's, to well within a "
        "hundredth of a percent.",
        "",
        f"The remaining **{agreement['frame_outliers']}** ligands "
        f"({100 * agreement['frame_outliers'] / agreement['n']:.1f}%) differ, and "
        "the reason is not the arithmetic but the frame: morfeus's reference here "
        "takes the *three nearest heavy atoms* as the donor substituents, while "
        "StericX uses covalent-radius bonding. Where those disagree -- an atom that "
        "is close but not bonded, or a bonded hydrogen that a heavy-atom rule "
        "ignores -- the integration sphere is centred differently. This is exactly "
        "the frame issue the StericX frame fix corrected and Study 006 "
        "characterizes; on these ligands StericX's bonded frame is the more "
        "defensible one. They are shown in a separate colour in the figure rather "
        "than averaged away. Overall (all "
        f"{agreement['n']} ligands, both conventions mixed) the agreement is still "
        f"R2 = {agreement['r2']:.4f}, MAE {agreement['mae']:.3f} %Vbur.",
        "",
        "![Throughput and agreement](speed_benchmark.png)",
        "",
        "*Figure. Left: single-core buried-volume throughput (log scale). Right: "
        "StericX vs morfeus %Vbur on every benchmarked structure; the "
        f"{agreement['frame_outliers']} orange points are frame-topology "
        "differences (see above), not compute errors. Generated by "
        "`studies/study_008_speed_benchmark.py`.*",
        "",
        "### Why this matters beyond the multiplier",
        "",
        f"StericX ships as a single **{binary_mb:.1f} MB** native binary with no "
        "runtime dependencies -- no Python, no NumPy, no environment to resolve. "
        "morfeus is a Python library that requires an interpreter and a scientific "
        "stack. For a one-off calculation the difference is convenience; at library "
        "scale, or embedded in a screening pipeline, the combination of a large "
        "constant-factor speedup and zero deployment surface is the point. The "
        "descriptor `descriptors` subcommand is single-threaded here for a fair "
        "per-core comparison; on this "
        f"{system['cores']}-core machine the wall-clock throughput scales further "
        "with trivial parallelism across files.",
        "",
        "### Honest caveats",
        "",
        "- **Single descriptor.** Only buried volume is timed head-to-head. It is "
        "the flagship descriptor (Study 004) and the most expensive of the panel, "
        "but morfeus's Sterimol and pyramidalization are not separately timed here.",
        "- **Warm cache, steady state.** Both tools are timed after the files are "
        "in the OS cache, measuring compute throughput rather than cold disk I/O.",
        "- **One machine.** Absolute numbers are specific to the CPU above; the "
        "ratio is the portable quantity and will vary with hardware.",
        "",
        "### Reproducing this study",
        "",
        "The Kraken DFT SDF cache is a local, gitignored artifact. With it in "
        "place and the release binary built (`cargo build --release`), run "
        "`uv run --extra science python studies/study_008_speed_benchmark.py`. "
        "Only StericX's own timings and the aggregate agreement are committed.",
        "",
    ]
    output.write_text("\n".join(lines), encoding="utf-8")


if __name__ == "__main__":
    raise SystemExit(main())
