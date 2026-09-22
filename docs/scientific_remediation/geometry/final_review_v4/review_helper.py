"""Review corrected captures against the sealed original; no new tolerance gate.

Reports all status changes, successful nonfinite outputs, all atom-order/rigid
groups, and full private coverage. Numeric differences are measurements, never
automatically permitted residuals. Independent reference review remains separate.
"""

import argparse
import itertools
import json
import math
import struct
from collections import Counter
from pathlib import Path

from sealed_reference import refuse_audit_output, sha

BASELINE_SHA256 = "1a4463460daebdd37a17412b56ebc959e75df08e000a03f62156b8f8f5178806"
LANES = [
    "geometry",
    "geometry_focused",
    "geometry_alignment",
    "geometry_alignment_rotation",
    "geometry_numerical",
    "geometry_historical",
    "kraken",
    "kraken_all_bins",
    "kraken_topology",
    "kraken_topology_all_bins",
]
METHODS = [
    "center",
    "sterimol_bond",
    "sterimol_dummy_raw",
    "sterimol_coordination",
    "pyramidalization",
    "buried_volume",
    "dump",
]


def verify(record):
    path = Path(record["path"])
    if sha(path) != record["sha256"] or path.stat().st_size != record["bytes"]:
        raise ValueError(f"Changed captured file: {path}")
    return path


def rows(path):
    with path.open() as handle:
        for line in handle:
            yield json.loads(line)


def error(value):
    return value.get("error") if isinstance(value, dict) else None


def classification(old, new, lane):
    old_error, new_error = error(old), error(new)
    if new_error:
        if "ambiguous at input coordinate precision" in new_error:
            return "fixed_probe_coordinate_precision_ambiguity"
        if any(
            text in new_error
            for text in [
                "invalid",
                "non-finite",
                "not finite",
                "coincid",
                "distinct",
                "out of bounds",
                "finite f32",
                "empty",
                "degenerate pyramidalization",
                "two distinct axis",
            ]
        ):
            return "checked_invalid_or_unrepresentable_geometry"
    if old_error and not new_error:
        if "symmetric zero" in old_error:
            return "valid_symmetric_volume"
        if "collinear with the coordination axis" in old_error:
            return "resolved_nonzero_small_plane"
        if "topology" in lane and "found 4 bonded substituents" in old_error:
            return "authoritative_topology_resolves_extra_contact"
    return "requires_review"


def scalar_changes(before, after, maxima, id_, key):
    if not isinstance(before, dict) or not isinstance(after, dict):
        return
    for field, value in after.items():
        old = before.get(field)
        if isinstance(value, (int, float)) and isinstance(old, (int, float)):
            delta = abs(value - old)
            label = key + "." + field
            if delta > maxima.get(label, {}).get("absolute_change", -1):
                maxima[label] = {"id": id_, "absolute_change": delta}


def finite_and_bits(row, method, value):
    if method in ["center", "dump"] or not isinstance(value, dict) or error(value):
        return []
    problems = []
    if value.get("_nonfinite"):
        problems.append(
            {"id": row["id"], "method": method, "nonfinite": value["_nonfinite"]}
        )
    for field, bits in value.get("_bits", {}).items():
        scalar = value.get(field)
        if not isinstance(scalar, (int, float)) or not math.isfinite(scalar):
            problems.append(
                {"id": row["id"], "method": method, "field": field, "value": scalar}
            )
        elif struct.pack(">f", scalar).hex() != bits:
            raise ValueError(f"Public float/bit mismatch: {row['id']} {method}.{field}")
    return problems


def compare_lane(name, before, after):
    source = name.replace("kraken_topology", "kraken")
    old_lane, new_lane = before["lanes"][source], after["lanes"][name]
    if old_lane["returncode"] or new_lane["returncode"]:
        raise ValueError("Nonzero observer process status")
    old_output, new_output = verify(old_lane["stdout"]), verify(new_lane["stdout"])
    old_requests, new_requests = (
        verify(old_lane["requests"]),
        verify(new_lane["requests"]),
    )
    verify(new_lane["stderr"])
    counts, changed, classifications = Counter(), Counter(), Counter()
    changes, nonfinite, maxima = [], [], {}
    seen = set()
    full_bin_lane = name in ["kraken_all_bins", "kraken_topology_all_bins"]
    for old, new, old_request, request in itertools.zip_longest(
        rows(old_output), rows(new_output), rows(old_requests), rows(new_requests)
    ):
        if any(value is None for value in [old, new, old_request, request]):
            raise ValueError(f"Mismatched lane lengths: {name}")
        id_ = new["id"]
        if len({old["id"], id_, old_request["id"], request["id"]}) != 1 or id_ in seen:
            raise ValueError(f"Mismatched or duplicate ID: {name} {id_}")
        seen.add(id_)
        if "topology" in name:
            projected = dict(request)
            projected["op"] = old_request.get("op", "geometry")
            if projected != old_request:
                raise ValueError("Topology lane changed more than explicit operation")
        elif request != old_request:
            raise ValueError("Original request meaning changed")
        counts["rows"] += 1
        for method in METHODS:
            if method not in old and method not in new:
                continue
            old_value, value = old.get(method), new.get(method)
            if (method in old) != (method in new):
                raise ValueError(f"Missing result section: {name} {id_} {method}")
            if error(value):
                counts[method + ".errors"] += 1
            else:
                counts[method + ".success"] += 1
            if old_value != value:
                changed[method] += 1
            if bool(error(old_value)) != bool(error(value)):
                reason = classification(old_value, value, name)
                classifications[reason] += 1
                changes.append(
                    {
                        "id": id_,
                        "method": method,
                        "before_error": error(old_value),
                        "after_error": error(value),
                        "classification": reason,
                    }
                )
            nonfinite.extend(finite_and_bits(new, method, value))
            if method != "dump":
                scalar_changes(old_value, value, maxima, id_, method)
        dump = new.get("dump")
        if isinstance(dump, dict) and not error(dump):
            frames = dump["orientations"]
            if len(frames) != 3:
                raise ValueError("Private observation lacks three frames")
            counts["private_frames"] += 3
            for frame in frames:
                if len(frame["quadrants"]) != 4 or len(frame["octants"]) != 8:
                    raise ValueError("Missing private bins")
                if full_bin_lane and any(
                    value is None or not math.isfinite(value)
                    for value in frame["quadrants"] + frame["octants"]
                ):
                    raise ValueError("Nonfinite full-corpus private bins")
                counts["private_bins"] += 12
        elif full_bin_lane and not error(new["buried_volume"]):
            raise ValueError("Successful public BV lacks private full-corpus bins")
    if counts["rows"] != new_lane["rows"] or not counts["rows"]:
        raise ValueError("Empty or false declared coverage")
    return {
        "counts": counts,
        "changed_rows": changed,
        "status_classifications": classifications,
        "status_changes": changes,
        "successful_nonfinite_outputs": nonfinite,
        "observed_changes_not_tolerances": maxima,
    }


def invariance(manifest):
    lane = manifest["lanes"]["geometry"]
    requests = list(rows(verify(lane["requests"])))
    observations = {row["id"]: row for row in rows(verify(lane["stdout"]))}
    reports = {}
    for category in ["permutation", "rigid"]:
        counts, changed, maxima, status = Counter(), Counter(), {}, []
        for request in requests:
            if request.get("category") != category:
                continue
            id_, parent = request["id"], request["parent"]
            counts[parent] += 1
            for method in METHODS[1:-1]:
                old, new = observations[parent][method], observations[id_][method]
                if old != new:
                    changed[method] += 1
                if bool(error(old)) != bool(error(new)):
                    status.append(
                        {
                            "id": id_,
                            "method": method,
                            "parent_error": error(old),
                            "case_error": error(new),
                        }
                    )
                scalar_changes(old, new, maxima, id_, method)
        reports[category] = {
            "groups": counts,
            "changed_rows": changed,
            "acceptance_changes": status,
            "observed_changes_not_tolerances": maxima,
        }
    return reports


def main(args):
    refuse_audit_output(args.output)
    args.output.mkdir(parents=True, exist_ok=False)
    baseline = args.baseline / "manifest.json"
    candidate = args.observations / "manifest.json"
    if sha(baseline) != BASELINE_SHA256:
        raise ValueError("Original baseline manifest differs from sealed v3 receipt")
    if sha(candidate) != candidate.with_suffix(".json.sha256").read_text().strip():
        raise ValueError("Candidate observation receipt hash mismatch")
    before, after = json.loads(baseline.read_text()), json.loads(candidate.read_text())
    if not after["complete_capture"] or not set(LANES) <= set(after["lanes"]):
        raise ValueError("Incomplete scientific capture")
    verify(after["binary"])
    build_path = verify(after["build_manifest"])
    build = json.loads(build_path.read_text())
    if build["binaries"]["observer"] != after["binary"]:
        raise ValueError("Observation binary is not linked to the recorded build")
    for record in build["adapter_files"]:
        verify(record)
    for source in build["sources"]:
        verify(source)
    result = {
        "kind": "corrected_geometry_review",
        "global_scientific_pass": False,
        "inputs": {
            str(p): sha(p) for p in [baseline, candidate, build_path, Path(__file__)]
        },
        "lanes": {name: compare_lane(name, before, after) for name in LANES},
        "invariance": invariance(after),
        "reference_status": (
            "Independent residual review required separately; "
            "no measured change is a tolerance"
        ),
    }
    (args.output / "review.json").write_text(
        json.dumps(result, indent=2, allow_nan=False) + "\n"
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ["baseline", "observations", "output"]:
        parser.add_argument("--" + name, type=Path, required=True)
    main(parser.parse_args())
