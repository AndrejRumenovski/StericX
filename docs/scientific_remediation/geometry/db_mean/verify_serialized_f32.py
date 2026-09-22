#!/usr/bin/env python3
"""Supplement the retained v4 witness check with its declared f32 output type."""

import json
from pathlib import Path

from recheck_v4 import f32, record

HERE = Path(__file__).resolve().parent
RUN = HERE / "final_v4"


def main():
    receipt = json.loads((RUN / "receipt.json").read_text())
    for item in receipt["bindings"] + receipt["files"]:
        if record(Path(item["path"])) != item:
            raise ValueError("retained witness input/output changed")
    failed = [key for key, passed in receipt["checks"].items() if not passed]
    if failed != ["distance_exact_high", "distance_exact_low"]:
        raise ValueError("unexpected original check failure")
    values = json.loads((RUN / "search_varying.stdout").read_text())["hits"]
    comparisons = [
        {
            "ligand": hit["ligand"],
            "serialized_decimal": hit["distance"],
            "observed_f32": f32(hit["distance"]),
            "analytic_expected_f32": receipt["expected_distances_f32"][hit["ligand"]],
            "exact_f32_equality": f32(hit["distance"])
            == receipt["expected_distances_f32"][hit["ligand"]],
        }
        for hit in values
    ]
    result = {
        "original_receipt": record(RUN / "receipt.json"),
        "helpers": [record(Path(__file__)), record(HERE / "recheck_v4.py")],
        "reason": (
            "Original checker compared shortest-roundtrip f32 JSON decimal parsed "
            "as Python f64 directly to its full f32 numeric value. Decode at the "
            "declared f32 API boundary, then require exact equality. No tolerance."
        ),
        "comparisons": comparisons,
        "other_original_checks_passed": all(
            passed for key, passed in receipt["checks"].items() if key not in failed
        ),
        "all_arithmetic_checks_passed": all(
            row["exact_f32_equality"] for row in comparisons
        ),
        "original_failed_receipt_preserved": True,
    }
    with (RUN / "serialized_f32_comparison.json").open("x") as stream:
        json.dump(result, stream, indent=2)
        stream.write("\n")
    if not result["all_arithmetic_checks_passed"]:
        raise SystemExit(2)


if __name__ == "__main__":
    main()
