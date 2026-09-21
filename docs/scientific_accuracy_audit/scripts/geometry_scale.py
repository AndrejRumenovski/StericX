# ruff: noqa: E501
"""Descriptor scale comparison only: repository library values are NOT accuracy truth."""

import csv
import hashlib
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
source = ROOT / "frozen/repository/data/ligand_db/kraken_phosphines.csv"
if not source.exists():
    raise RuntimeError("must use frozen library")
rows = list(csv.DictReader(source.open()))
out = {}
changes = {
    "sterimol_l": 1.0013580322265625e-5,
    "sterimol_b1": 0.019944190979003906,
    "sterimol_b5": 4.76837158203125e-6,
    "percent_buried_volume": 0.14926569734644,
    "buried_volume": 0.2680721841369511,
    "max_delta_qvbur": 0.1491403579711914,
    "pyr_p": 3.874301910400391e-6,
    "pyr_alpha": 0.0002593994140625,
}
for field, deviation in changes.items():
    x = np.array([float(r[field]) for r in rows])
    ordered = np.sort(x)
    gaps = np.diff(ordered)
    nearest = np.minimum(np.r_[np.inf, gaps], np.r_[gaps, np.inf])
    iqr = float(np.quantile(x, 0.75) - np.quantile(x, 0.25))
    nn = float(np.median(nearest))
    distinct = float(np.median(gaps[gaps > 0]))
    out[field] = {
        "N": len(x),
        "iqr": iqr,
        "median_univariate_nearest_neighbor_spacing": nn,
        "median_positive_adjacent_spacing": distinct,
        "measured_deviation": deviation,
        "deviation_over_iqr": deviation / iqr,
        "deviation_over_median_nearest_spacing": deviation / nn if nn else None,
        "deviation_type": "rigid-transformation observed maximum"
        if field in ["sterimol_l", "sterimol_b1", "sterimol_b5", "pyr_p", "pyr_alpha"]
        else "default lens error"
        if field in ["percent_buried_volume", "buried_volume"]
        else "default-versus-very-fine maximum in4molecule campaign",
    }
alignment = json.loads((ROOT / "geometry/alignment/reference_results.json").read_text())
alignment_scale = []
for case in alignment["comparison"]:
    if not case["id"].endswith("__controlled"):
        continue
    for short, field in [("l", "sterimol_l"), ("b1", "sterimol_b1"), ("b5", "sterimol_b5")]:
        deviation = abs(case["errors"][short])
        scale = out[field]
        alignment_scale.append({
            "id": case["id"],
            "descriptor": field,
            "absolute_error": deviation,
            "error_over_iqr": deviation / scale["iqr"],
            "error_over_median_nearest_spacing": deviation / scale["median_univariate_nearest_neighbor_spacing"],
            "meaning": "Controlled real Kraken case; B1 combines scan and alignment effects. Not an experimental effect size.",
        })
(ROOT / "geometry/scientific_scale.json").write_text(
    json.dumps(
        {
            "source": str(source.relative_to(ROOT)),
            "source_sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
            "interpretation": "These rounded model-library outputs are only an empirical descriptor scale, never a correctness target. Univariate nearest spacings in a dense1541row database are not chemical effect sizes or experimental resolution.",
            "fields": out,
            "additional_targeted_real_geometry_observations": alignment_scale,
        },
        indent=2,
    )
    + "\n"
)
print(json.dumps(out, indent=2))
