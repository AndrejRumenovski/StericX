"""Recover maximum available public coordinate precision for the entire cohort."""
import concurrent.futures
import datetime
import gzip
import hashlib
import json
import threading
import time

import numpy as np
import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

from analyze import read_jsonl
from freeze import OUT, save, sha

LOCAL = threading.local()


def acquire(req):
    if not hasattr(LOCAL, "session"):
        LOCAL.session = requests.Session()
        LOCAL.session.mount("https://", HTTPAdapter(max_retries=Retry(total=3, backoff_factor=.5, status_forcelist=[429, 500, 502, 503, 504], allowed_methods=["GET"])))
    cid = req["id"].split(":")[-1]
    url = f"https://descriptor-libraries.molssi.org/api/kraken/conformers/export/xyz/{cid}"
    row = {"id": req["id"], "url": url, "started_utc": datetime.datetime.now(datetime.UTC).isoformat()}
    try:
        r = LOCAL.session.get(url, timeout=30)
        row.update(status=r.status_code, headers=dict(r.headers), body=r.text,
            sha256=hashlib.sha256(r.content).hexdigest(), bytes=len(r.content))
        r.raise_for_status()
        lines = r.text.splitlines()
        n = int(lines[0])
        data = [line.split() for line in lines[2:2+n]]
        xyz = np.array([[float(v) for v in a[1:4]] for a in data])
        expected = np.array([a["position"] for a in req["atoms"]])
        assert [a[0] for a in data] == [a["element"] for a in req["atoms"]]
        assert xyz.shape == expected.shape and np.isfinite(xyz).all()
        difference = xyz - expected
        row.update(coordinates_equal=bool(np.array_equal(xyz, expected)),
            max_coordinate_difference_A=float(np.max(np.abs(difference))),
            RMSD_A=float(np.sqrt(np.mean(np.sum(difference**2, axis=1)))),
            consistent_with_SDF_rounding=bool(np.max(abs(difference)) <= .000050000001),
            coordinates=xyz.tolist())
    except Exception as error:
        row["error"] = repr(error)
    row["finished_utc"] = datetime.datetime.now(datetime.UTC).isoformat()
    return row


def main():
    out = OUT / "recovered_xyz"
    out.mkdir(exist_ok=False)
    requests_ = list(read_jsonl(OUT / "frozen/requests.jsonl.gz"))
    save(out / "plan.json", {"created_utc": datetime.datetime.now(datetime.UTC).isoformat(),
        "script_sha256": sha(__file__), "N": len(requests_), "workers": 8,
        "reason": "Primary-source XYZ checks revealed six-decimal coordinates for ligand 821 that its SDF export rounds to four decimals.",
        "selection": "Every one of the frozen 31611 conformer IDs, not a target-residual subset.",
        "admission": "Same atom identity/order and every coordinate differs by <=0.000050000001 A: admit as recoverable export precision. Larger differences or failures are explicit and require investigation, never silently replace input or remove ligand.",
        "predicted_effect": "Smaller coordinate-rounding contribution; no expectation or selection based on R2. Exact historical precision/membership remains separately unverified.",
        "no_reference_value_access_in_acquisition": True})
    start = time.monotonic()
    errors, changed, incompatible = [], [], []
    with gzip.open(out / "responses.jsonl.gz", "wt", compresslevel=3) as dst, concurrent.futures.ThreadPoolExecutor(max_workers=8) as pool:
        for i, r in enumerate(pool.map(acquire, requests_), 1):
            dst.write(json.dumps(r, separators=(",", ":")) + "\n")
            if "error" in r:
                errors.append({"id": r["id"], "error": r["error"]})
            elif not r["consistent_with_SDF_rounding"]:
                incompatible.append(r["id"])
            elif not r["coordinates_equal"]:
                changed.append(r["id"])
            if i % 500 == 0:
                print(i, "of", len(requests_), f"{time.monotonic()-start:.1f}s", "changed", len(changed), "errors", len(errors), flush=True)
    save(out / "complete.json", {"N": len(requests_), "changed": changed, "errors": errors,
        "incompatible": incompatible, "elapsed_seconds": time.monotonic()-start,
        "responses_sha256": sha(out / "responses.jsonl.gz"), "plan_sha256": sha(out / "plan.json")})


if __name__ == "__main__":
    main()
