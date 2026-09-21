# ruff: noqa: E501
"""Freeze raw SUT responses first, then independently calculate comparisons."""

import argparse
import datetime
import hashlib
import json
import subprocess
import traceback
import warnings
from pathlib import Path

import numpy as np
from geometry_reference import bv_refs, pyramidalization_ref, sterimol_refs

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "geometry"


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def clean(x):
    if isinstance(x, float) and not np.isfinite(x):
        return str(x)
    if isinstance(x, dict):
        return {k: clean(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [clean(v) for v in x]
    if isinstance(x, np.generic):
        return clean(x.item())
    return x


def cases():
    return (
        json.loads((OUT / "inputs.json").read_text())["cases"]
        + json.loads((OUT / "supplement_inputs.json").read_text())["cases"]
    )


def freeze(binary):
    target = OUT / "sut_outputs.jsonl"
    manifest = OUT / "sut_freeze.json"
    if target.exists():
        old = json.loads(manifest.read_text())
        assert sha(target) == old["output_sha256"]
        assert sha(Path(binary)) == old["observer_sha256"]
        for name, expected in old["input_hashes"].items():
            assert sha(OUT / name) == expected
        print("Verified existing immutable SUT capture")
        return
    requests = []
    for c in cases():
        requests.append(
            c | {"op": "geometry", "dump": not c.get("sterimol_only", False)}
        )
    for c in json.loads((OUT / "supplement_inputs.json").read_text())["raw"]:
        requests.append(c | {"op": "raw_occupancy"})
    req = OUT / "requests.jsonl"
    req.write_text("".join(json.dumps(c) + "\n" for c in requests))
    with (
        req.open("rb") as src,
        target.open("wb") as dst,
        (OUT / "sut_stderr.txt").open("wb") as err,
    ):
        p = subprocess.run([binary], stdin=src, stdout=dst, stderr=err, check=False)
    data = {
        "acquired_utc": datetime.datetime.now(datetime.UTC).isoformat(),
        "observer_binary": str(binary),
        "observer_sha256": sha(Path(binary)),
        "request_sha256": sha(req),
        "output_sha256": sha(target),
        "returncode": p.returncode,
        "line_count": sum(1 for _ in target.open()),
        "input_hashes": {
            p.name: sha(p)
            for p in [OUT / "inputs.json", OUT / "supplement_inputs.json"]
        },
    }
    manifest.write_text(json.dumps(data, indent=2) + "\n")
    print(data)
    assert p.returncode == 0 and data["line_count"] == len(requests)


def compare():
    manifest = json.loads((OUT / "sut_freeze.json").read_text())
    assert sha(OUT / "sut_outputs.jsonl") == manifest["output_sha256"]
    for name, expected in manifest["input_hashes"].items():
        assert sha(OUT / name) == expected
    sut = {
        r["id"]: r
        for r in map(json.loads, (OUT / "sut_outputs.jsonl").read_text().splitlines())
        if "id" in r
    }
    path = OUT / "reference_results.jsonl"
    with path.open("w") as fp:
        for i, c in enumerate(cases()):
            s = sut[c["id"]]
            xyz = np.array([a["position"] for a in s["atoms"]])
            row = {
                "id": c["id"],
                "category": c["category"],
                "coordinate_mode": "identical f32-rounded coordinates converted to float64",
                "reference_radii": "Morfeus Bondi or explicit input override",
            }
            with warnings.catch_warnings(record=True) as ws:
                warnings.simplefilter("always")
                for key, fun in [
                    (
                        "sterimol",
                        lambda c=c, xyz=xyz, s=s: sterimol_refs(
                            c,
                            xyz,
                            controlled_center=s.get("center")
                            if isinstance(s.get("center"), list)
                            and all(v is not None for v in s["center"])
                            else None,
                        ),
                    ),
                    (
                        "sterimol_independent_center",
                        lambda c=c, xyz=xyz, s=s: sterimol_refs(c, xyz),
                    ),
                    (
                        "pyramidalization",
                        lambda c=c, xyz=xyz: pyramidalization_ref(c, xyz),
                    ),
                ]:
                    try:
                        row[key] = fun()
                    except Exception as e:
                        row[key] = {"error": str(e), "type": type(e).__name__}
                if not c.get("sterimol_only"):
                    center = s.get("center")
                    if isinstance(center, list) and all(v is not None for v in center):
                        try:
                            row["buried_volume"] = bv_refs(c, xyz, center=center)
                        except Exception as e:
                            row["buried_volume"] = {
                                "error": str(e),
                                "type": type(e).__name__,
                                "traceback": traceback.format_exc(),
                            }
                    else:
                        row["buried_volume"] = {
                            "error": "SUT coordination center unavailable; cannot make controlled-center comparison",
                            "sut_center": center,
                        }
                row["warnings"] = [str(w.message) for w in ws]
            fp.write(json.dumps(clean(row), allow_nan=False) + "\n")
            fp.flush()
            if i % 50 == 0:
                print(i, c["id"], flush=True)
    print("Reference results", sha(path))


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("mode", choices=["freeze", "compare"])
    p.add_argument("--binary", default=str(ROOT / "frozen/bin/stericx-audit-observer"))
    args = p.parse_args()
    if args.mode == "freeze":
        freeze(args.binary)
    else:
        compare()
