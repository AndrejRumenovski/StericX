#!/usr/bin/env python3
"""Render source-bound fixed-resolution evidence without scoring release gates."""

# Long generated Markdown paragraphs/table rows are single string literals.
# ruff: noqa: E501

import argparse
import csv
import hashlib
import itertools
import json
import math
import shutil
from fractions import Fraction
from pathlib import Path


def record(path):
    path = Path(path).resolve()
    with path.open("rb") as stream:
        digest = hashlib.file_digest(stream, "sha256").hexdigest()
    return {"path": str(path), "bytes": path.stat().st_size, "sha256": digest}


def save(path, obj):
    path.write_text(json.dumps(obj, indent=2) + "\n")


def upper_mean(values):
    exact = sum(map(Fraction, values)) / len(values)
    result = float(exact)
    return math.nextafter(result, math.inf) if Fraction(result) < exact else result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    run = args.run.resolve()
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    manifest = json.loads((run / "manifest.json").read_text())
    for item in [*manifest["files"], manifest["summary"]]:
        if record(item["path"]) != item:
            raise ValueError(f"Evidence changed: {item['path']}")
    cases = json.loads((run / "cases.json").read_text())
    summary = json.loads((run / "summary.json").read_text())
    references = {}
    rows = []
    for case in cases:
        directory = run / case["id"]
        references[case["id"]] = json.loads(
            (directory / "continuous_reference.json").read_text()
        )
        rows.extend(json.loads((directory / "results.json").read_text()))
    complete = [row for row in rows if "status" not in row]
    save(
        output / "convergence_results.json",
        {"cases": cases, "references": references, "rows": rows},
    )
    shutil.copyfile(run / "summary.json", output / "convergence_summary.json")
    columns = [
        "case",
        "method",
        "level",
        "control",
        "sphere_nodes",
        "occupied_nodes",
        "volume_angstrom3",
        "percent",
        "reference_lower",
        "reference_upper",
        "absolute_error_lower",
        "absolute_error_upper",
        "relative_error_lower",
        "relative_error_upper",
        "successive_volume_difference",
        "elapsed_seconds",
        "cost_scope",
        "cell_error_bound_angstrom3",
        "cell_interval_lower",
        "cell_interval_upper",
        "reference_interval_width",
        "sphere_mask_disagreements",
        "numerator_mask_disagreements",
        "native_disputed_nodes",
        "native_unexplained_nodes",
        "analytical_truth_in_cell_bound",
        "clamped",
    ]
    with (output / "convergence_residuals.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=columns)
        writer.writeheader()
        for row in complete:
            ref = references[row["case_id"]]
            interval = ref["analytical_interval"] or [ref["lower"], ref["upper"]]
            enclosure = row["cell_enclosure"]
            relative = row["relative_error_interval"] or [None, None]
            writer.writerow(
                dict(
                    zip(
                        columns,
                        [
                            row["case_id"],
                            row["method"],
                            row["level"],
                            row["control"],
                            row["sphere_nodes"],
                            row["occupied_nodes"],
                            row["volume"],
                            row["percent"],
                            *interval,
                            *row["absolute_error_interval"],
                            *relative,
                            row["successive_volume_difference"],
                            row["elapsed_seconds"],
                            row["cost_scope"],
                            enclosure["total_absolute_error_upper"],
                            enclosure["certified_volume_lower"],
                            enclosure["certified_volume_upper"],
                            ref["width"],
                            enclosure["denominator_reference_mask_disagreements"],
                            enclosure["numerator_reference_mask_disagreements"],
                            row.get("native_disputed_node_count"),
                            row.get("native_unexplained_disputed_node_count"),
                            row["analytical_truth_inside_cell_enclosure"],
                            row.get("clamped"),
                        ],
                        strict=True,
                    )
                )
            )
    methods = sorted({row["method"] for row in complete})
    analytics = [
        case["id"]
        for case in cases
        if references[case["id"]]["analytical_interval"] is not None
    ]
    stats = []
    for method in methods:
        for level in range(4):
            selected = [
                row
                for row in complete
                if row["method"] == method
                and row["level"] == level
                and row["case_id"] in analytics
            ]
            stats.append(
                {
                    "method": method,
                    "level": level,
                    "analytical_case_count": len(selected),
                    "maximum_absolute_error_upper": max(
                        row["absolute_error_interval"][1] for row in selected
                    ),
                    "mean_absolute_error_upper": upper_mean(
                        [row["absolute_error_interval"][1] for row in selected]
                    ),
                    "maximum_cell_error_bound": max(
                        row["cell_enclosure"]["total_absolute_error_upper"]
                        for row in selected
                    ),
                }
            )
    regressions = []
    for method in methods:
        for case in analytics:
            selected = sorted(
                [r for r in complete if r["method"] == method and r["case_id"] == case],
                key=lambda r: r["level"],
            )
            for before, after in itertools.pairwise(selected):
                if (
                    after["absolute_error_interval"][0]
                    > before["absolute_error_interval"][1]
                ):
                    regressions.append(
                        {
                            "case": case,
                            "method": method,
                            "from_level": before["level"],
                            "to_level": after["level"],
                            "previous_absolute_error": before[
                                "absolute_error_interval"
                            ],
                            "next_absolute_error": after["absolute_error_interval"],
                        }
                    )
    save(
        output / "convergence_trends.json",
        {
            "analytical_aggregate_statistics": stats,
            "nonmonotonic_error_steps": regressions,
        },
    )
    native_disputes = sum(row.get("native_disputed_node_count", 0) for row in complete)
    text = [
        "# Fixed-resolution buried-volume convergence\n",
        "This study preserves every scheduled comparison and does not itself score the approved release policy. The Rust production engine and all reference implementations remain unchanged.\n",
        f"The accepted campaign contains **{len(cases)} geometries, {len(rows)} comparisons**, and {len(summary['unsuccessful_comparisons'])} unsuccessful comparisons. It includes 17 mathematical primitives and three captured native frames of original ligand 401/0 at R=2.5 Å.\n",
        "## Methods and provenance\n",
        "Native and Morfeus density controls are 0.08, 0.01, 0.00125 and 0.00015625 Å³. DBSTEP spacings are 0.4, 0.2, 0.1 and 0.05 Å. These controls generate different lattices. All comparisons use identical accepted binary32 atomic coordinates and scaled radii promoted exactly to the reference precision; the native kernel retains its actual rounded squared radii.\n",
        "Morfeus is called through public `BuriedVolume`, with an excluded zero-radius ghost center at the origin, explicit radii, `radii_scale=1`, and hydrogen inclusion after shared explicit filtering. DBSTEP has no equivalent public custom-radii input: its unmodified `occupied`/`buried_vol` and `occupied_direct`/`buried_vol` routines are therefore identified as support-function lanes. Their actual rounded axes cover an explicit cube [-3,3]³ Å³. Separate unchanged public DBSTEP calls and their reference defect are documented in [DBSTEP_SOURCE_REVIEW.md](../DBSTEP_SOURCE_REVIEW.md).\n",
        "The [actual-cell derivation](../CONVERGENCE_CELL_BOUNDS.md) includes finite coverage, nonuniform axis widths, direct-reference roundoff, observed mask differences, count/weight conversion, ratio normalization and clamp arithmetic. DBSTEP fast denominator masks are explicitly reconstructed and count-checked, not independently observed point decisions. The bounds are conservative and often broad; they do not establish fine numerical precision.\n",
        "All raw native requests/responses, reference-returned points, axes, full masks, disputed-node witnesses, source snapshots, dependency hashes and elapsed costs are bound by the campaign manifest and [receipt](convergence_receipt.json). Native timings include observation calls and JSON transfer and are **not production latency benchmarks**.\n",
        "## Independent reference checks\n",
        f"Adaptive boxes use the prospectively fixed depths 8 (primitives) and 11 (original molecular frames). Analytical box-enclosure failures: {len(summary['analytical_box_failures'])}; non-nested interval failures: {len(summary['non_nested_box_failures'])}; analytical cell-bound failures: {len(summary['analytical_cell_bound_failures'])}.\n",
        f"Across convergence grids, {native_disputes} native disputed-node observations retain exact rational margins and explicit binary32 operation traces; unexplained native disputed nodes: {summary['native_unexplained_disputed_node_count']}. The separate 10-point boundary set contains {summary['exact_boundary_native_disagreements']} native/exact disagreements, with {len(summary['exact_boundary_unexplained'])} unexplained. These witnesses distinguish rounded point inclusion from the continuous integral.\n",
        "### Original ligand 401 continuous bounds\n",
        "| Frame | Lower volume (Å³) | Upper volume (Å³) | Interval width (Å³) | Reference cost (s) |",
        "|---|---:|---:|---:|---:|",
    ]
    for case, ref in references.items():
        if case.startswith("original401"):
            text.append(
                f"| {case} | {ref['lower']:.9f} | {ref['upper']:.9f} | {ref['width']:.9f} | {ref['elapsed_seconds']:.3f} |"
            )
    text += [
        "",
        "These intervals concern the continuous volume of the accepted balls, including each frame's stored coordinate rounding. Their widths must be respected when interpreting the historical ~0.024 Å³ difference. Full convergence rows and signed/absolute/relative error intervals are in [convergence_residuals.csv](convergence_residuals.csv); [convergence_results.json](convergence_results.json) retains settings, counts and bound components.\n",
        "| Frame | Native at density0.01 (Å³) | Proven absolute-error lower bound (Å³) | Native at finest density (Å³) |",
        "|---|---:|---:|---:|",
    ]
    for case in references:
        if case.startswith("original401"):
            default = next(
                row
                for row in complete
                if row["case_id"] == case
                and row["method"] == "stericx_native"
                and row["level"] == 1
            )
            fine = next(
                row
                for row in complete
                if row["case_id"] == case
                and row["method"] == "stericx_native"
                and row["level"] == 3
            )
            text.append(
                f"| {case} | {default['volume']:.9f} | {default['absolute_error_interval'][0]:.9f} | {fine['volume']:.9f} |"
            )
    text += [
        "",
        "A positive error lower bound proves that the corresponding default-density estimate differs from the continuous integral. Finite-grid correctness does not imply continuum exactness or the historical arithmetic-only 0.00025 Å³ accuracy. The finite-width molecular intervals do not support ranking differences substantially smaller than their widths.\n",
        "### Maximum analytical error by level\n",
        "These maxima summarize the fixed synthetic analytical cases; they are not a universal accuracy claim or a molecule benchmark.\n",
        "| Method | Level | Maximum absolute error upper bound (Å³) | Maximum conservative cell error bound (Å³) |",
        "|---|---:|---:|---:|",
    ]
    for row in stats:
        text.append(
            f"| {row['method']} | {row['level']} | {row['maximum_absolute_error_upper']:.9f} | {row['maximum_cell_error_bound']:.6f} |"
        )
    text += [
        "",
        f"There are {len(regressions)} rigorously increasing analytical-error steps among the finer successive grids. Every such step is retained in [convergence_trends.json](convergence_trends.json). Center-point integration is not monotonic; finer grids are not automatically closer for every geometry.\n",
        "## Boundary/convention findings\n",
        "For exact external tangency (evaluation sphere R=2.5 Å, atom r=1 Å at x=3.5 Å), the continuous intersection is zero. Some native/DBSTEP grids include one closed-boundary contact node and consequently report a small nonzero finite-grid estimate. Morfeus's public prefilter uses `atom.radius + sphere.radius > norm(atom.coordinates)` and removes the exactly tangent atom, reporting zero. This difference first occurs in Morfeus's prefilter, before point occupancy; it is a finite-grid boundary convention with the same zero continuous-volume limit. It does not justify changing validated Rust defaults.\n",
        "The original strict 401 comparison failure remains a separate finite-domain reclip issue. This convergence campaign neither erases its failed historical criterion nor fits a tolerance to that residual. No production correction was introduced.\n",
        "## Remaining limits\n",
        "The shape tests are mathematical inputs, not characterized experimental molecules. Density and spacing are distinct numerical controls. Continuous box intervals remain finite-width; shell bounds overcount overlapping boundaries and can be much larger than observed analytical errors. Region summaries and molecular rigid motions are evaluated in the separate point-domain/adversarial campaigns. Cross-platform installation and matched performance are separate release requirements.\n",
    ]
    (output / "CONVERGENCE_REPORT.md").write_text("\n".join(text))
    copies = [record(path) for path in sorted(output.iterdir()) if path.is_file()]
    save(
        output / "convergence_receipt.json",
        {
            "kind": "convergence_evidence_summary",
            "release_gate_scored": False,
            "controller": record(__file__),
            "campaign_manifest": record(run / "manifest.json"),
            "campaign_summary": record(run / "summary.json"),
            "committed_outputs": copies,
            "raw_evidence_count": len(manifest["files"]),
        },
    )
    print(
        json.dumps(
            {
                "output": str(output),
                "rows": len(rows),
                "native_disputed_nodes": native_disputes,
                "nonmonotonic_steps": len(regressions),
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
