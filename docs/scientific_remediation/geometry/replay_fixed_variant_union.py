#!/usr/bin/env python3
"""Replay the union of sealed historical and fresh Morfeus diagnostic cases.

Only complete, byte-sealed fresh request/SUT pairs can supply reused equations.
The historical unchanged run() determines every reference variant. Its inferred
neighbor order remains explicit even in the topology SUT lane; new three-plane
means are separate derived observations, not substitutions in the raw results.
"""

import argparse
import importlib.metadata
import itertools
import json
import platform
import shutil
import sys
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

from sealed_reference import AUDIT, MANIFEST, ROOT, refuse_audit_output, sha, verify

VARIANTS = [
    "matched_sut_001",
    "matched_sut_0001",
    "primary_001",
    "primary_0001",
    "ster_primary_center_bondi",
    "ster_primary_center_paton",
    "ster_sut_center_sut_radii",
    "ster_sut_center_paton",
    "pyr_f32_topology",
    "pyr_f64_topology",
    "pyr_f64_primary_neighbors",
]


def rows(path):
    with path.open() as handle:
        for line in handle:
            yield json.loads(line)


def save(path, value):
    with path.open("x") as handle:
        json.dump(value, handle, indent=2, allow_nan=False)
        handle.write("\n")


def canonical_request(request):
    value = dict(request)
    value["op"] = "geometry"
    for name in ["dump", "include_points"]:
        if value.pop(name, False):
            raise ValueError("Unexpected private dump in historical reference input")
    return value


def sealed_campaign(path, original=False):
    paths = [
        path / name
        for name in [
            "selected_inputs.jsonl",
            "selection_frozen.json",
            "reference_raw.jsonl",
            "manifest_frozen.json",
        ]
    ]
    if original:
        for item in paths:
            verify(item)
    selection, manifest = [json.loads(paths[i].read_text()) for i in [1, 3]]
    if sha(paths[0]) != selection.get(
        "selected_inputs_sha256", selection.get("selected_input_sha256")
    ):
        raise ValueError(f"Changed selected pairs: {path}")
    if sha(paths[2]) != manifest.get(
        "reference_raw_sha256", manifest.get("raw_sha256")
    ):
        raise ValueError(f"Changed reference equations: {path}")
    if sha(paths[1]) != manifest.get(
        "selection_manifest_sha256", manifest.get("selection_sha256")
    ):
        raise ValueError(f"Changed selection manifest: {path}")
    pairs = {row["request"]["id"]: row for row in rows(paths[0])}
    results = {row["id"]: row for row in rows(paths[2])}
    if (
        len(pairs) != selection["N"]
        or len(results) != selection["N"]
        or pairs.keys() != results.keys()
    ):
        raise ValueError("Incomplete campaign identity coverage")
    return pairs, results, {str(p): sha(p) for p in paths}


def main(args):
    refuse_audit_output(args.output)
    args.output.mkdir(parents=True, exist_ok=False)
    shutil.copy2(__file__, args.output / Path(__file__).name)
    shutil.copy2(
        Path(__file__).with_name("sealed_reference.py"),
        args.output / "sealed_reference.py",
    )
    historical = AUDIT / "kraken"
    reference_script = historical / "scripts/kraken_morfeus_outliers.py"
    primary_script = historical / "raw_sources/official_PL_dft_library_201027.body"
    for path in [reference_script, primary_script]:
        verify(path)
    sys.path.insert(0, str(reference_script.parent))
    from kraken_morfeus_outliers import run

    inputs = {
        str(p): sha(p)
        for p in [
            MANIFEST,
            reference_script,
            primary_script,
            Path(__file__),
            Path(__file__).with_name("sealed_reference.py"),
            ROOT / "uv.lock",
            ROOT / "pyproject.toml",
        ]
    }
    memberships, union, fresh = {}, {}, {}
    for label, base in [
        ("historical", historical),
        ("inferred", args.inferred),
        ("topology", args.topology),
    ]:
        fresh[label] = {}
        for campaign in ["morfeus_outliers", "delta_outliers"]:
            pairs, results, identities = sealed_campaign(
                base / campaign, label == "historical"
            )
            inputs.update(identities)
            memberships[label + "/" + campaign] = sorted(pairs)
            for id_, pair in pairs.items():
                request = canonical_request(pair["request"])
                if id_ in union and union[id_] != request:
                    raise ValueError(f"Scientific request drift for {id_}")
                union[id_] = request
                if label != "historical":
                    if id_ in fresh[label] and fresh[label][id_] != (
                        pair,
                        results[id_],
                    ):
                        raise ValueError("Fresh reference sources disagree")
                    fresh[label][id_] = (pair, results[id_])
    observation_manifest = args.observations / "manifest.json"
    if (
        sha(observation_manifest)
        != observation_manifest.with_suffix(".json.sha256").read_text().strip()
    ):
        raise ValueError("Changed observation manifest")
    inputs[str(observation_manifest)] = sha(observation_manifest)
    observations = json.loads(observation_manifest.read_text())
    captured = {}
    for mode, lane_name in [("inferred", "kraken"), ("topology", "kraken_topology")]:
        lane = observations["lanes"][lane_name]
        for field in ["requests", "stdout"]:
            item = lane[field]
            path = Path(item["path"])
            if sha(path) != item["sha256"] or path.stat().st_size != item["bytes"]:
                raise ValueError("Changed observation stream")
            inputs[str(path)] = item["sha256"]
        selected = {}
        count = 0
        for req, obs in itertools.zip_longest(
            rows(Path(lane["requests"]["path"])), rows(Path(lane["stdout"]["path"]))
        ):
            if req is None or obs is None or req["id"] != obs["id"]:
                raise ValueError("Observation identity mismatch")
            count += 1
            if req["id"] in union:
                if canonical_request(req) != union[req["id"]]:
                    raise ValueError("Reference request differs from captured geometry")
                if req["id"] in selected:
                    raise ValueError("Duplicate observation ID")
                selected[req["id"]] = {"request": union[req["id"]], "sut": obs}
        if count != 31721 or selected.keys() != union.keys():
            raise ValueError("Incomplete observed union")
        captured[mode] = selected
    started = {
        "kind": "fixed_historical_and_adaptive_union_reference",
        "global_scientific_pass": False,
        "inputs": inputs,
        "membership_counts": {k: len(v) for k, v in memberships.items()},
        "memberships": memberships,
        "union_count": len(union),
        "argv": sys.argv,
        "packages": {
            k: importlib.metadata.version(k) for k in ["numpy", "scipy", "morfeus-ml"]
        },
        "python": sys.version,
        "platform": platform.platform(),
        "request_projection": (
            "Only op and false private-dump defaults normalized; all chemical "
            "inputs unchanged. Original unchanged reference function ignores "
            "op/dump. Each lane uses its exact captured SUT values."
        ),
        "reference_contract": (
            "Unchanged original run(): density0.01 and0.001; captured-center/"
            "captured-radii variants, primary raw-vector-center/Bondi variants, "
            "four Sterimol and three pyramidalization variants. Its matched_sut "
            "plane order uses inferred bonded_neighbors even for explicit-topology "
            "SUT; all plane identities retained."
        ),
    }
    save(args.output / "started.json", started)
    summaries = {}
    for mode in ["inferred", "topology"]:
        directory = args.output / mode
        directory.mkdir()
        pairs, reusable, pending = captured[mode], {}, []
        with (directory / "selected_inputs.jsonl").open("x") as handle:
            for id_ in sorted(union):
                pair = pairs[id_]
                handle.write(
                    json.dumps(pair, separators=(",", ":"), allow_nan=False) + "\n"
                )
                previous = fresh[mode].get(id_)
                if previous is not None and previous[0] == pair:
                    reusable[id_] = previous[1]
                else:
                    pending.append(pair)
        print(
            mode,
            "union",
            len(union),
            "reuse exact fresh pairs",
            len(reusable),
            "new",
            len(pending),
            flush=True,
        )
        calculations = {}
        with ProcessPoolExecutor(max_workers=args.workers) as executor:
            for index, value in enumerate(executor.map(run, pending, chunksize=1), 1):
                calculations[value["id"]] = value
                if index % 10 == 0:
                    print(mode, index, "/", len(pending), flush=True)
        results = reusable | calculations
        if results.keys() != union.keys():
            raise ValueError("Missing union reference result")
        coverage = {name: Counter() for name in VARIANTS}
        with (
            (directory / "reference_raw.jsonl").open("x") as raw,
            (directory / "mean_recombination.jsonl").open("x") as means,
        ):
            for id_ in sorted(union):
                row, pair = results[id_], pairs[id_]
                raw.write(
                    json.dumps(row, separators=(",", ":"), allow_nan=False) + "\n"
                )
                reduced = {"id": id_, "variants": {}}
                for name in VARIANTS:
                    value = row.get(name)
                    status = (
                        "missing"
                        if value is None
                        else "error"
                        if "error" in value
                        else "success"
                    )
                    coverage[name][status] += 1
                    if status != "success" or "orientations" not in value:
                        continue
                    frames = value["orientations"]
                    planes = [frame["plane_atom"] for frame in frames]
                    wanted = pair["request"]["neighbors"]
                    matches = len(frames) == 3 and set(planes) == set(wanted)
                    entry = {
                        "plane_atoms": planes,
                        "exact_three_topological_planes": matches,
                    }
                    if matches:
                        average = {
                            key: sum(frame[key] for frame in frames) / 3
                            for key in ["buried_volume", "near_vbur", "far_vbur"]
                        }
                        entry["independent_three_plane_mean"] = average
                        public = pair["sut"]["buried_volume"]
                        entry["sut_minus_variant"] = (
                            {
                                key: public[key] - reference
                                for key, reference in average.items()
                            }
                            if "error" not in public
                            else None
                        )
                    else:
                        entry["reason_not_recombined"] = (
                            "Historical reference plane set is not the three "
                            "supplied topology neighbors"
                        )
                    reduced["variants"][name] = entry
                means.write(json.dumps(reduced, allow_nan=False) + "\n")
        summaries[mode] = {
            "rows": len(results),
            "reused_exact_fresh_pairs": len(reusable),
            "new_calculations": len(calculations),
            "coverage": {k: dict(v) for k, v in coverage.items()},
            "artifacts": {p.name: sha(p) for p in directory.iterdir()},
        }
    if any(sha(Path(path)) != digest for path, digest in inputs.items()):
        raise ValueError("Input or reference bytes changed during replay")
    save(
        args.output / "complete.json", started | {"complete": True, "lanes": summaries}
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ["observations", "inferred", "topology", "output"]:
        parser.add_argument("--" + name, type=Path, required=True)
    parser.add_argument("--workers", type=int, default=3)
    main(parser.parse_args())
