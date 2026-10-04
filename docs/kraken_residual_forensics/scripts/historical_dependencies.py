"""Archive date-selected upstream code, without scoring dependency variants."""
import concurrent.futures
import datetime
import json

import requests

from freeze import OUT, save, sha


def main():
    out = OUT / "historical_dependencies"
    out.mkdir(exist_ok=False)
    save(out / "plan.json", {
        "script_sha256": sha(__file__),
        "selection": "Last upstream commit before each independently dated public milestone: preprint 2021-04-27 and journal online publication 2022-01-12. These are inspection snapshots, NOT claimed Kraken dependency pins. No variant selected using target values.",
        "files": ["sterimol.py", "buried_volume.py", "pyramidalization.py", "geometry.py", "data.py"],
    })
    records = []

    def fetch(name, url):
        r = requests.get(url, timeout=40)
        p = out / (name + ".body")
        p.write_bytes(r.content)
        record = {"name": name, "url": url, "status": r.status_code,
                  "sha256": sha(p), "fetched_utc": datetime.datetime.now(datetime.UTC).isoformat()}
        save(out / (name + ".metadata.json"), record)
        return record, r

    meta, response = fetch("updated_repo", "https://api.github.com/repos/SigmanGroup/kraken")
    records.append(meta)
    if response.status_code == 200:
        branch = response.json()["default_branch"]
        meta, _ = fetch("updated_tree", f"https://api.github.com/repos/SigmanGroup/kraken/git/trees/{branch}?recursive=1")
        records.append(meta)
    repo = "digital-chemistry-laboratory/morfeus"
    jobs = []
    selected = {}
    for label, date in [("preprint", "2021-04-27"), ("journal", "2022-01-12")]:
        meta, response = fetch(label + "_commit", f"https://api.github.com/repos/{repo}/commits?until={date}T23:59:59Z&per_page=1")
        records.append(meta)
        if response.status_code == 200 and response.json():
            commit = response.json()[0]["sha"]
            selected[label] = commit
            for file in ["sterimol.py", "buried_volume.py", "pyramidalization.py", "geometry.py", "data.py"]:
                jobs.append((label + "_" + file, f"https://raw.githubusercontent.com/{repo}/{commit}/morfeus/{file}"))
    with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
        for meta, _ in pool.map(lambda args: fetch(*args), jobs):
            records.append(meta)
    save(out / "complete.json", {"selected_commits": selected, "responses": records,
        "historical_dependency_pin": None, "no_descriptor_values_calculated": True})
    print(json.dumps({"selected_commits": selected, "statuses": [(r["name"], r["status"]) for r in records]}))


if __name__ == "__main__":
    main()
