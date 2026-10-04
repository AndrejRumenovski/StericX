"""Compare every scientific record/bit/private bin after optional API addition."""
import json
from itertools import zip_longest
from pathlib import Path

from freeze import FINAL, OUT, save, sha


def main():
    root = OUT / "candidates/angular_sampling"
    before_path = FINAL / "observations/manifest.json"
    after_path = root / "full_observations/manifest.json"
    before, after = [json.loads(p.read_text()) for p in (before_path, after_path)]
    assert set(before["lanes"]) == set(after["lanes"])
    rows = []
    for name in sorted(before["lanes"]):
        a, b = before["lanes"][name], after["lanes"][name]
        for r in (a["stdout"], b["stdout"], a["requests"], b["requests"]):
            assert sha(r["path"]) == r["sha256"]
        assert a["requests"]["sha256"] == b["requests"]["sha256"], name
        differences = []
        count = 0
        with open(a["stdout"]["path"]) as left, open(b["stdout"]["path"]) as right:
            for old, new in zip_longest(left, right):
                count += 1
                if old != new:
                    differences.append({"row": count, "id": json.loads(old)["id"] if old else None,
                                        "dictionary_equal": bool(old and new and json.loads(old) == json.loads(new))})
        rows.append({"lane": name, "N": count, "raw_line_differences": differences,
            "before_sha256": a["stdout"]["sha256"], "after_sha256": b["stdout"]["sha256"],
            "private_frames": b.get("private_frames", 0), "private_bins": b.get("private_bins", 0),
            "exit_code": b["returncode"]})
        print(name, count, "differences", len(differences), flush=True)
    save(root / "full_replay_comparison.json", {"script_sha256": sha(__file__),
        "before_manifest_sha256": sha(before_path), "after_manifest_sha256": sha(after_path),
        "lanes": rows, "passed": all(not r["raw_line_differences"] and r["exit_code"] == 0 for r in rows),
        "total_rows": sum(r["N"] for r in rows), "numeric_tolerance": 0,
        "scope": "All original scientific observer lanes, including complete 31721-conformer inferred/topology public/private runs; atom-order, rigid-transform, symmetric, non-P and invalid inputs; kinetics and model observations."})


if __name__ == "__main__":
    main()
