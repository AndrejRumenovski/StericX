"""Cover new final-stage top100 extremizers without changing any cohort input."""
import concurrent.futures
import gzip
import json

from analyze import read_jsonl
from freeze import OUT, save, sha
from minimal_reference import protect


def main():
    out = OUT / "minimal_reference_final"
    out.mkdir(exist_ok=False)
    needed = set(json.loads((OUT / "minimal_reference/final_top100_additional_ids.json").read_text()))
    save(out / "plan.json", {"script_sha256": sha(__file__), "ids": sorted(needed),
        "selection": "New final-stage top100 min/max/selected conformers not already covered by the sealed independent minimal calculation; same criterion as original selection"})
    requests = [r for r in read_jsonl(OUT / "campaigns/recovered_precision/requests.jsonl.gz") if r["id"] in needed]
    refs = {r["id"]: r for r in read_jsonl(OUT / "campaigns/recovered_precision/morfeus_observations.jsonl.gz") if r["id"] in needed}
    plan = json.loads((OUT / "conventions/plan.json").read_text())
    radii = next(r["value"] for r in plan["constant_audit"] if r["quantity"] == "BV_radii")
    with concurrent.futures.ProcessPoolExecutor(max_workers=4) as pool:
        rows = list(pool.map(protect, ((r, {"primary_center": refs[r["id"]]["center"], "primary_neighbors": refs[r["id"]]["neighbors"]}, radii) for r in requests)))
    with gzip.open(out / "added_observations.jsonl.gz", "wt") as dst:
        for r in rows:
            dst.write(json.dumps(r, separators=(",", ":")) + "\n")
    n = 0
    with gzip.open(out / "observations.jsonl.gz", "wt") as dst:
        for r in read_jsonl(OUT / "campaigns/recovered_precision/minimal_observations.jsonl.gz"):
            dst.write(json.dumps(r, separators=(",", ":")) + "\n")
            n += 1
        for r in rows:
            dst.write(json.dumps(r, separators=(",", ":")) + "\n")
            n += 1
    save(out / "complete.json", {"N": n, "added_N": len(rows), "errors": [r for r in rows if "error" in r],
        "observations_sha256": sha(out / "observations.jsonl.gz"), "plan_sha256": sha(out / "plan.json")})


if __name__ == "__main__":
    main()
