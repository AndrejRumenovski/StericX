"""Complete independent top-100 coverage for the separately named old target."""
import concurrent.futures
import gzip
import json

import pandas as pd

from analyze import read_jsonl
from freeze import OUT, save, sha
from minimal_reference import protect


def main():
    out = OUT / "minimal_reference_historical"
    out.mkdir(exist_ok=False)
    needed = set()
    for stage in ("baseline", "final"):
        for file in (OUT / "historical_target" / stage / "rankings").glob("*.csv"):
            for r in pd.read_csv(file).head(100).to_dict("records"):
                for key in ("min_conformer_id", "max_conformer_id", "selected_conformer_id"):
                    if pd.notna(r[key]):
                        needed.add(f"KRAKEN:{r['molecule_id']}:{int(r[key])}")
    old = list(read_jsonl(OUT / "minimal_reference_final/observations.jsonl.gz"))
    missing = needed-{r["id"] for r in old}
    save(out / "plan.json", {"script_sha256": sha(__file__), "needed": sorted(needed), "additional": sorted(missing),
        "criterion": "Every min/max/selected conformer in baseline and final top100 of every historical-reference metric; same criterion as API coverage."})
    requests = [r for r in read_jsonl(OUT / "campaigns/recovered_precision/requests.jsonl.gz") if r["id"] in missing]
    refs = {r["id"]: r for r in read_jsonl(OUT / "campaigns/recovered_precision/morfeus_observations.jsonl.gz") if r["id"] in missing}
    plan = json.loads((OUT / "conventions/plan.json").read_text())
    radii = next(r["value"] for r in plan["constant_audit"] if r["quantity"] == "BV_radii")
    with concurrent.futures.ProcessPoolExecutor(max_workers=3) as pool:
        rows = list(pool.map(protect, ((r, {"primary_center": refs[r["id"]]["center"], "primary_neighbors": refs[r["id"]]["neighbors"]}, radii) for r in requests)))
    assert all("error" not in r for r in rows)
    with gzip.open(out / "observations.jsonl.gz", "wt") as dst:
        for r in old+rows:
            dst.write(json.dumps(r, separators=(",", ":")) + "\n")
    save(out / "complete.json", {"N": len(old)+len(rows), "added_N": len(rows), "errors": [],
        "required_N": len(needed), "missing_after": sorted(needed-{r['id'] for r in old+rows}),
        "observations_sha256": sha(out / "observations.jsonl.gz")})


if __name__ == "__main__":
    main()
