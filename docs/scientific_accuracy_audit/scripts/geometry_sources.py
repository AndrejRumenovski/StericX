"""Acquire primary/official geometry references and preserve exact HTTP bodies."""

import datetime
import hashlib
import json
from pathlib import Path

import requests

ROOT = Path(__file__).resolve().parents[1] / "geometry" / "sources"
ROOT.mkdir(exist_ok=True)
if (ROOT / "manifest.json").exists():
    frozen = json.loads((ROOT / "manifest.json").read_text())
    for source in frozen:
        if "sha256" in source:
            body = ROOT / (source["id"] + ".body")
            assert hashlib.sha256(body.read_bytes()).hexdigest() == source["sha256"]
    print("Existing frozen geometry source captures verified; no downloads replaced.")
    raise SystemExit(0)
SOURCES = {
    "verloop1976": "https://doi.org/10.1016/B978-0-12-060307-7.50010-9",
    "radhakrishnan1991": "https://link.springer.com/article/10.1007/BF00676621",
    "bondi1964": "https://pubs.acs.org/doi/10.1021/j100785a001",
    "cordero2008": "https://pubs.rsc.org/en/content/articlehtml/2008/dt/b801115j",
    "cordero2008_si": "https://www.rsc.org/suppdata/dt/b8/b801115j/b801115j.txt",
    "sambvca_manual": "https://www.aocdweb.com/OMtools/sambvca2.1/help/help.html",
    "morfeus_sterimol": "https://digital-chemistry-laboratory.github.io/morfeus/sterimol.html",
    "morfeus_buried_volume": "https://digital-chemistry-laboratory.github.io/morfeus/buried_volume.html",
    "morfeus_pyramidalization": "https://digital-chemistry-laboratory.github.io/morfeus/pyramidalization.html",
    "aarontools_sterimol": "https://aarontools.readthedocs.io/en/latest/api/substituent.html",
}
rows = []
for name, url in SOURCES.items():
    try:
        response = requests.get(url, timeout=45)
        data = response.content
        (ROOT / (name + ".body")).write_bytes(data)
        rows.append(
            {
                "id": name,
                "requested_url": url,
                "url": response.url,
                "status": response.status_code,
                "sha256": hashlib.sha256(data).hexdigest(),
                "bytes": len(data),
                "acquired_utc": datetime.datetime.now(datetime.UTC).isoformat(),
                "content_type": response.headers.get("content-type"),
            }
        )
    except Exception as e:
        rows.append({"id": name, "requested_url": url, "error": str(e)})
(ROOT / "manifest.json").write_text(json.dumps(rows, indent=2) + "\n")
print(json.dumps(rows, indent=2))
