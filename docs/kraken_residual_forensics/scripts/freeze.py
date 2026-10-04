"""Freeze the requested historical cohort and validated current implementation.

Run once. Refuses to replace an existing freeze. No scientific calculation or
reference-value selection is performed here.
"""
from __future__ import annotations

import csv
import datetime
import gzip
import hashlib
import json
import shutil
import subprocess
import tarfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "docs/kraken_residual_forensics"
AUDIT = ROOT / "docs/scientific_accuracy_audit/kraken"
FINAL = ROOT / ".stericx/profiling/scientifically_validated_optimization/final_post_optimization_v1"


def sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def save(path, value):
    Path(path).write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n")


def main():
    out = OUT / "frozen"
    out.mkdir(exist_ok=False)
    historical = ROOT / "docs/study_004/kraken_dft_scaled_comparison.csv"
    ids = sorted(int(r["Source_ID"]) for r in csv.DictReader(historical.open()))
    assert len(ids) == len(set(ids)) == 1541
    cohort = set(ids)
    save(out / "cohort.json", ids)
    sources = {}

    def copy(src, name):
        src = Path(src)
        sources[str(src.relative_to(ROOT))] = {"sha256": sha(src), "bytes": src.stat().st_size}
        shutil.copy2(src, out / name)

    binary = ROOT / "target/release/stericx"
    assert sha(binary) == "b233501640e4d06555a36a278581fbf9cb011962311f205c81f3fc9a83ae89e6"
    copy(binary, "stericx")
    copy(historical, "historical_1541_predictions.csv")
    copy(ROOT / "data/official/ni_hda_kraken.csv", "historical_reference.csv")
    copy(ROOT / "docs/scientific_remediation/CONVENTIONS_AND_LIMITS.md", "current_conventions.md")
    copy(ROOT / "docs/performance/POST_OPTIMIZATION_SCIENTIFIC_EQUIVALENCE.md", "validated_executable_receipt.md")
    copy(FINAL / "observations/manifest.json", "current_observations_manifest.json")
    copy(AUDIT / "primary/ligands_manifest.json", "ligands_manifest.json")
    copy(AUDIT / "primary/conformers_manifest.json", "conformers_manifest.json")
    copy(AUDIT / "raw_sources/official_PL_dft_library_201027.body", "kraken_dft_source.py.txt")
    copy(AUDIT / "raw_sources/official_PL_dft_library_201027.metadata.json", "kraken_dft_source.metadata.json")
    copy(AUDIT / "full_analytic_reference/manifest_frozen.json", "prior_morfeus_manifest.json")
    manifest = json.loads((FINAL / "observations/manifest.json").read_text())
    source_observation = Path(manifest["lanes"]["kraken_topology"]["stdout"]["path"])
    assert sha(source_observation) == manifest["lanes"]["kraken_topology"]["stdout"]["sha256"]
    copy(Path(manifest["binary"]["path"]), "current_observer")
    assert sha(out / "current_observer") == manifest["binary"]["sha256"]
    stream_sources = {
        "requests.jsonl.gz": AUDIT / "prepared_all/requests.jsonl",
        "inventory.jsonl.gz": AUDIT / "prepared_all/inventory.jsonl",
        "current_predictions.jsonl.gz": source_observation,
        "prior_morfeus_analytic.jsonl.gz": AUDIT / "full_analytic_reference/reference_raw.jsonl",
    }
    counts = {}
    for name, src in stream_sources.items():
        sources[str(src.relative_to(ROOT))] = {"sha256": sha(src), "bytes": src.stat().st_size}
        n = 0
        with gzip.open(out / name, "wt") as dst, src.open() as stream:
            for line in stream:
                row = json.loads(line)
                if "id" in row and int(row["id"].split(":")[1]) in cohort:
                    dst.write(line)
                    n += 1
        counts[name] = n
    assert len(set(counts.values())) == 1, counts
    geometries = []
    with tarfile.open(out / "inputs_and_reference.tar.gz", "w:gz") as archive:
        for row in map(json.loads, gzip.open(out / "inventory.jsonl.gz", "rt")):
            src = AUDIT / row["sdf_path"]
            assert sha(src) == row["sdf_sha256"]
            archive.add(src, arcname=row["sdf_path"], recursive=False)
            geometries.append({"id": row["id"], "path": row["sdf_path"], "sha256": row["sdf_sha256"]})
        for mid in ids:
            for directory in ("primary/published_dft", "primary/molecule"):
                src = AUDIT / directory / f"{mid}.body"
                if src.exists():
                    archive.add(src, arcname=f"{directory}/{mid}.body", recursive=False)
                    sources[str(src.relative_to(ROOT))] = {"sha256": sha(src), "bytes": src.stat().st_size}
    save(out / "geometry_hashes.json", geometries)
    tracked = subprocess.check_output(["git", "ls-files", "-z"], cwd=ROOT).decode().split("\0")
    with tarfile.open(out / "source_and_tests.tar.gz", "w:gz") as archive:
        for name in tracked:
            if name and (name.startswith(("src/", "tests/", "scripts/", "studies/", ".github/")) or name in ("Cargo.toml", "Cargo.lock", "pyproject.toml", "uv.lock")):
                src = ROOT / name
                archive.add(src, arcname=name, recursive=False)
                sources[name] = {"sha256": sha(src), "bytes": src.stat().st_size}
    save(out / "manifest.json", {
        "created_utc": datetime.datetime.now(datetime.UTC).isoformat(),
        "source_commit": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
        "git_status": subprocess.check_output(["git", "status", "--short"], cwd=ROOT, text=True),
        "cohort_source": str(historical.relative_to(ROOT)), "ligands": 1541,
        "rows": counts, "source_files": sources,
        "baseline": "Current validated accepted-C2 executable; previously measured full precision explicit-topology outputs authenticated against final observation receipt. Historical pre-remediation predictions preserved separately.",
        "scientific_scope": "14 native descriptor fields x min/max/delta/vburminconf; percent Vbur is explicitly derived from absolute reference. No weights inferred from published values.",
        "exact_geometry_warning": "API SDF exports have 4-decimal coordinates; their identity to historical full-precision calculation geometries is not established by their DFT label.",
        "files": {p.name: {"sha256": sha(p), "bytes": p.stat().st_size} for p in sorted(out.iterdir()) if p.is_file()},
    })
    print(json.dumps({"ligands": len(ids), "rows": counts, "manifest_sha256": sha(out / "manifest.json")}))


if __name__ == "__main__":
    main()
