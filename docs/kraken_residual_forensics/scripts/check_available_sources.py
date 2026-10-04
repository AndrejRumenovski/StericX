"""Fresh bounded primary-source checks; retain responses and failures verbatim."""
import concurrent.futures
import datetime
import hashlib
import json
import time

import numpy as np
import pandas as pd
import requests

from analyze import read_jsonl
from freeze import AUDIT, OUT, save, sha


def fetch(job):
    name, url = job
    start = datetime.datetime.now(datetime.UTC).isoformat()
    path = OUT / "source_availability" / name
    try:
        r = requests.get(url, timeout=30)
        path.with_suffix(".body").write_bytes(r.content)
        row = {"name": name, "url": url, "final_url": r.url, "status": r.status_code,
            "started_utc": start, "finished_utc": datetime.datetime.now(datetime.UTC).isoformat(),
            "sha256": hashlib.sha256(r.content).hexdigest(), "bytes": len(r.content),
            "headers": dict(r.headers)}
    except requests.RequestException as error:
        row = {"name": name, "url": url, "started_utc": start, "error": repr(error)}
    save(path.with_suffix(".metadata.json"), row)
    return row


def main():
    out = OUT / "source_availability"
    out.mkdir(exist_ok=False)
    keys, mids = set(), set()
    rankings = OUT / "campaigns/primary_grid/historical_last_plane/rankings"
    for name in ("max_delta_qvbur_min", "buried_volume_delta", "sterimol_l_min", "sterimol_l_max", "sterimol_b1_min", "sterimol_b1_max", "sterimol_b5_min", "sterimol_b5_max", "pyr_p_min", "pyr_p_max", "pyr_alpha_min", "pyr_alpha_max"):
        for r in pd.read_csv(rankings / (name + ".csv")).head(5).to_dict("records"):
            mids.add(int(r["molecule_id"]))
            for k in ("min_conformer_id", "max_conformer_id", "selected_conformer_id"):
                if pd.notna(r[k]):
                    keys.add(f"KRAKEN:{r['molecule_id']}:{int(r[k])}")
    api = "https://descriptor-libraries.molssi.org/api/kraken"
    jobs = [(f"molecule_{mid}", f"{api}/molecules/{mid}") for mid in sorted(mids)]
    for key in sorted(keys):
        cid = key.split(":")[-1]
        jobs.extend([(f"xyz_{cid}", f"{api}/conformers/export/xyz/{cid}"), (f"energy_{cid}", f"{api}/conformers/data/{cid}")])
    jobs += [("original_releases", "https://api.github.com/repos/the-matter-lab/kraken/releases"),
             ("original_tree", "https://api.github.com/repos/the-matter-lab/kraken/git/trees/master?recursive=1"),
             ("updated_tree", "https://api.github.com/repos/SigmanGroup/kraken/git/trees/main?recursive=1")]
    save(out / "plan.json", {"created_utc": datetime.datetime.now(datetime.UTC).isoformat(), "script_sha256": sha(__file__),
        "selected_ids": sorted(keys), "ligands": sorted(mids), "jobs": jobs,
        "scope": "Current XYZ precision, membership, and energy availability for leading residual cases; original source repository release/tree check. Does not prove data absent from authors' private archives."})
    with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
        responses = list(pool.map(fetch, jobs))
    old = {r["id"]: r for r in read_jsonl(OUT / "frozen/requests.jsonl.gz") if r["id"] in keys}
    records = []
    for key in sorted(keys):
        cid = key.split(":")[-1]
        row = {"id": key}
        xyzpath = out / f"xyz_{cid}.body"
        try:
            lines = xyzpath.read_text().splitlines()
            n = int(lines[0])
            data = [line.split() for line in lines[2:2+n]]
            elems = [r[0] for r in data]
            xyz = np.array([[float(v) for v in r[1:4]] for r in data])
            req = old[key]
            expected = np.array([a["position"] for a in req["atoms"]])
            row["atoms_and_order_equal"] = elems == [a["element"] for a in req["atoms"]]
            row["coordinates_equal"] = bool(np.array_equal(xyz, expected))
            row["coordinate_RMSD_A"] = float(np.sqrt(np.mean(np.sum((xyz-expected)**2, axis=1))))
            row["more_than_4_decimal_precision_present"] = bool(np.any(abs(xyz-np.round(xyz, 4)) > 1e-12))
        except Exception as error:
            row["coordinate_check_error"] = repr(error)
        try:
            energy = json.loads((out / f"energy_{cid}.body").read_text())
            row["energy_data"] = energy.get("data")
            row["energy_response_keys"] = list(energy)
        except Exception as error:
            row["energy_check_error"] = repr(error)
        records.append(row)
    save(out / "coordinate_checks.json", records)
    inherited = {}
    for name in ("energy_source_search/REPORT.md", "energy_source_search/manifest.json", "raw_sources/kraken_si_zip_listing.json", "primary_si/manifest.json", "delta_interpretation/ligand_369_dossier.json"):
        p = AUDIT / name
        if p.exists():
            inherited[name] = {"path": str(p), "sha256": sha(p)}
    save(out / "complete.json", {"responses": responses, "coordinate_checks": len(records),
        "all_successfully_checked_geometries_identical": all(r.get("coordinates_equal", False) for r in records),
        "inherited_source_availability_evidence": inherited,
        "limit": "Historical Gaussian logs, retained-ID lists, per-conformer corrected free energies, and dependency pins have not been recovered. Current metadata alone cannot verify class A."})
    print(json.dumps({"requests": len(responses), "coordinate_checks": len(records), "unchanged": sum(r.get("coordinates_equal", False) for r in records), "failures": sum(r.get("status") != 200 for r in responses)}))


if __name__ == "__main__":
    main()
