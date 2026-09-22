"""Compare every captured topology BV scalar/frame/bin with fresh references.

All residuals retain full Python float precision. B, Br and other substantive
radius-convention differences are separate strata. Metrics are observations,
not fitted tolerances, and every per-ID/per-plane/bin residual is retained.
"""

import argparse
import heapq
import importlib.metadata
import itertools
import json
import math
import struct
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
from sealed_reference import (
    MANIFEST,
    REFERENCE,
    load_reference,
    refuse_audit_output,
    sha,
)

reference = load_reference()


def verify(record):
    path = Path(record["path"])
    if sha(path) != record["sha256"] or path.stat().st_size != record["bytes"]:
        raise ValueError(f"Changed captured file: {path}")
    return path


def rows(path):
    with path.open() as handle:
        for line in handle:
            yield json.loads(line)


def main(args):
    refuse_audit_output(args.output)
    args.output.mkdir(parents=True, exist_ok=False)
    obs_manifest = args.observations / "manifest.json"
    if (
        sha(obs_manifest)
        != obs_manifest.with_suffix(".json.sha256").read_text().strip()
    ):
        raise ValueError("Observation manifest seal mismatch")
    observed = json.loads(obs_manifest.read_text())
    lane = observed["lanes"]["kraken_topology_all_bins"]
    requests, outputs = verify(lane["requests"]), verify(lane["stdout"])
    ref_receipt = args.references / "complete.json"
    receipt = json.loads(ref_receipt.read_text())
    references = args.references / "reference_results.jsonl"
    if (
        not receipt["complete"]
        or not receipt["complete_source_stream_validated"]
        or receipt["limit"] is not None
        or receipt["requests"] != 31721
        or receipt["reference_errors"]
    ):
        raise ValueError(
            "Independent reference scope is not the complete successful corpus"
        )
    if sha(references) != receipt["reference_results_sha256"]:
        raise ValueError("Changed independent reference stream")
    for path, expected in receipt["inputs"].items():
        if sha(Path(path)) != expected:
            raise ValueError(f"Changed reference input/helper: {path}")
    inputs = {
        str(p): sha(p)
        for p in [
            obs_manifest,
            requests,
            outputs,
            ref_receipt,
            references,
            MANIFEST,
            REFERENCE,
            Path(__file__),
        ]
    }
    stats, top = defaultdict(list), defaultdict(list)
    strata, mismatch_elements = Counter(), Counter()
    count = frame_count = bin_count = 0
    seen = set()

    def record(stratum, field, error, id_, actual, expected, plane=None, index=None):
        if not all(math.isfinite(value) for value in [error, actual, expected]):
            raise ValueError(f"Nonfinite comparison: {id_} {field}")
        for group in ["all", stratum]:
            key = group + "/" + field
            stats[key].append(error)
            item = (
                abs(error),
                id_,
                -1 if plane is None else plane,
                -1 if index is None else index,
                actual,
                expected,
            )
            if len(top[key]) < 20:
                heapq.heappush(top[key], item)
            elif item > top[key][0]:
                heapq.heapreplace(top[key], item)

    with (args.output / "per_id_residuals.jsonl").open("x") as destination:
        for request, observation, independent in itertools.zip_longest(
            rows(requests), rows(outputs), rows(references)
        ):
            if any(value is None for value in [request, observation, independent]):
                raise ValueError("Reference/private observation lengths differ")
            id_ = request["id"]
            if observation["id"] != id_ or independent["id"] != id_ or id_ in seen:
                raise ValueError("Reference/private observation IDs differ or repeat")
            seen.add(id_)
            native_radii = [atom["vdw_radius"] for atom in observation["atoms"]]
            expected_radii = reference.radii_for(request)
            mismatches = []
            for index, (actual, expected) in enumerate(
                zip(native_radii, expected_radii, strict=True)
            ):
                # Same f32 rounding is arithmetic precision, not a distinct table.
                if struct.pack(">f", actual) != struct.pack(">f", expected):
                    element = request["atoms"][index]["element"].upper()
                    mismatches.append(
                        {
                            "index": index,
                            "element": element,
                            "sut_radius": actual,
                            "independent_radius": float(expected),
                        }
                    )
                    mismatch_elements[element] += 1
            elements = {atom["element"].upper() for atom in request["atoms"]}
            stratum = (
                ("B" if "B" in elements else "no_B")
                + "+"
                + ("Br" if "BR" in elements else "no_Br")
            )
            other = sorted({item["element"] for item in mismatches} - {"B", "BR"})
            if other:
                stratum += "+other_radius:" + ",".join(other)
            strata[stratum] += 1
            public = observation["buried_volume"]
            expected = independent["buried_volume"]["independent"]
            if "error" in public or "reference_error" in independent:
                raise ValueError("Unexpected unsuccessful full-corpus comparison")
            residual = {
                field: public[field] - value for field, value in expected.items()
            }
            if residual != independent["residual"]:
                raise ValueError(
                    "Private/public captures differ for the referenced public values"
                )
            for field, value in residual.items():
                record(
                    stratum,
                    "public." + field,
                    value,
                    id_,
                    public[field],
                    expected[field],
                )
            dump = observation["dump"]
            sut_frames = {frame["plane"]: frame for frame in dump["orientations"]}
            ref_frames = {
                frame["plane"]: frame["independent"]
                for frame in independent["buried_volume"]["orientations"]
            }
            if len(sut_frames) != 3 or sut_frames.keys() != ref_frames.keys():
                raise ValueError("Private/reference plane identities differ")
            frames = []
            for plane, native in sut_frames.items():
                expected_frame = ref_frames[plane]
                if dump["grid_count"] != expected_frame["grid_count"]:
                    raise ValueError("Default-grid sample counts differ")
                frame = {"plane": plane, "grid_count": dump["grid_count"]}
                for name, size in [("quadrants", 4), ("octants", 8)]:
                    actual_values, expected_values = native[name], expected_frame[name]
                    if len(actual_values) != size or len(expected_values) != size:
                        raise ValueError("Missing private/reference bin")
                    frame[name] = []
                    for index, (actual, expected_value) in enumerate(
                        zip(actual_values, expected_values, strict=True)
                    ):
                        delta = actual - expected_value
                        frame[name].append(delta)
                        record(
                            stratum,
                            "private." + name,
                            delta,
                            id_,
                            actual,
                            expected_value,
                            plane,
                            index,
                        )
                        bin_count += 1
                for name in ["buried_volume", "near_vbur", "far_vbur"]:
                    delta = native[name] - expected_frame[name]
                    frame[name] = delta
                    record(
                        stratum,
                        "private." + name,
                        delta,
                        id_,
                        native[name],
                        expected_frame[name],
                        plane,
                    )
                frames.append(frame)
                frame_count += 1
            destination.write(
                json.dumps(
                    {
                        "id": id_,
                        "radius_stratum": stratum,
                        "radius_mismatches": mismatches,
                        "public": residual,
                        "frames": frames,
                    },
                    allow_nan=False,
                )
                + "\n"
            )
            count += 1
    if (count, frame_count, bin_count) != (31721, 95163, 1141956):
        raise ValueError("Full public/frame/bin coverage failed")
    metrics = {}
    for key, values in stats.items():
        array = np.array(values, dtype=np.float64)
        absolute = np.abs(array)
        metrics[key] = {
            "count": len(values),
            "max_abs": float(absolute.max()),
            "mean_signed": float(array.mean()),
            "mean_abs": float(absolute.mean()),
            "rms": float(np.sqrt(np.mean(array * array))),
            "quantiles_abs": {
                str(q): float(np.quantile(absolute, q)) for q in [0.5, 0.9, 0.99, 0.999]
            },
            "zero_count": int(np.count_nonzero(array == 0)),
        }
    top_rows = {
        key: [
            {
                "absolute_error": value[0],
                "id": value[1],
                "plane": None if value[2] == -1 else value[2],
                "index": None if value[3] == -1 else value[3],
                "sut": value[4],
                "independent": value[5],
                "residual": value[4] - value[5],
            }
            for value in sorted(values, reverse=True)
        ]
        for key, values in top.items()
    }
    (args.output / "metrics.json").write_text(
        json.dumps(metrics, indent=2, allow_nan=False) + "\n"
    )
    (args.output / "top20.json").write_text(
        json.dumps(top_rows, indent=2, allow_nan=False) + "\n"
    )
    manifest = {
        "kind": "complete_private_bin_reference_comparison",
        "global_scientific_pass": False,
        "thresholds_inferred_from_observed_errors": False,
        "inputs": inputs,
        "rows": count,
        "frames": frame_count,
        "bins": bin_count,
        "radius_strata": strata,
        "radius_mismatch_atom_counts": mismatch_elements,
        "packages": {
            name: importlib.metadata.version(name)
            for name in ["numpy", "scipy", "morfeus-ml"]
        },
        "artifacts": {
            name: sha(args.output / name)
            for name in ["per_id_residuals.jsonl", "metrics.json", "top20.json"]
        },
    }
    (args.output / "complete.json").write_text(json.dumps(manifest, indent=2) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ["references", "observations", "output"]:
        parser.add_argument("--" + name, type=Path, required=True)
    main(parser.parse_args())
