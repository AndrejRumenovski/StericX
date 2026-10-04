"""Final integrity and coverage audit; no target-based acceptance decision."""
import datetime
import hashlib
import json
import re
import subprocess

import numpy as np
import pandas as pd

from analyze import FIELDS, metrics, read_jsonl
from freeze import OUT, ROOT, save, sha


def main():
    assert not (OUT / "verification.json").exists()
    checks = []
    frozen = json.loads((OUT / "frozen/manifest.json").read_text())
    assert sha(OUT / "frozen/manifest.json") == "73ee1518d8ea548662307f3b78fa96435509ef69343d0909e0eae87a6463bc69"
    for name, info in frozen["files"].items():
        assert sha(OUT / "frozen" / name) == info["sha256"], name
    changed = []
    for name, info in frozen["source_files"].items():
        if sha(ROOT / name) != info["sha256"]:
            changed.append(name)
    assert changed == ["src/geometry/sterimol.rs"], changed
    checks.append({"check": "All frozen files and source receipts intact", "source_files": len(frozen["source_files"]), "authorized_changed_source": changed})
    cohort = set(json.loads((OUT / "frozen/cohort.json").read_text()))
    assert len(cohort) == 1541
    requests = {r["id"] for r in read_jsonl(OUT / "frozen/requests.jsonl.gz")}
    recovered = {r["id"] for r in read_jsonl(OUT / "campaigns/recovered_precision/requests.jsonl.gz")}
    refs = list(read_jsonl(OUT / "campaigns/recovered_precision/morfeus_observations.jsonl.gz"))
    assert len(refs) == len(requests) == 31611 and recovered == requests == {r["id"] for r in refs}
    assert all("error" not in r for r in refs)
    native = pd.read_csv(OUT / "campaigns/recovered_precision/analysis/conformers.csv.gz", float_precision="round_trip")
    assert len(native) == 31611 and set(native.id) == requests and set(native.molecule_id) == cohort
    assert np.isfinite(native[list(FIELDS)]).all().all()
    checks.append({"check": "Full native and independent-reference input coverage", "ligands": 1541, "conformers": 31611, "failures": 0})
    needed = set()
    for name in ["baseline", "results/final_native", "historical_target/baseline", "historical_target/final"]:
        p = OUT / name
        d = pd.read_csv(p / "all_residuals.csv.gz", float_precision="round_trip")
        met = pd.read_csv(p / "metrics.csv", float_precision="round_trip").set_index("metric")
        assert len(d) == 1541*56 and len(met) == 56
        for metric, g in d.groupby("metric"):
            assert len(g) == 1541 and set(g.molecule_id) == cohort
            calculated = metrics(g.kraken, g.stericx)
            for key, value in calculated.items():
                assert np.isclose(met.loc[metric, key], value, rtol=1e-13, atol=1e-14), (name, metric, key)
            r = pd.read_csv(p / "rankings" / (metric+".csv"), float_precision="round_trip")
            assert len(r) == 1541 and np.array_equal(r["rank"], np.arange(1, 1542))
            assert (np.diff(r.absolute_residual) <= 0).all()
            for row in r.head(100).to_dict("records"):
                for k in ["min_conformer_id", "max_conformer_id", "selected_conformer_id"]:
                    if pd.notna(row[k]):
                        needed.add(f"KRAKEN:{row['molecule_id']}:{int(row[k])}")
        assert len(list((p / "plots").glob("*.png"))) == 56
        top = json.loads((p / "top_10_25_50_100.json").read_text())
        assert len(top) == 224 and {r["top_n"] for r in top} == {10, 25, 50, 100}
    independent = list(read_jsonl(OUT / "minimal_reference_historical/observations.jsonl.gz"))
    ids = {r["id"] for r in independent}
    assert len(independent) == len(ids) == 5491 and not (needed-ids)
    assert all("error" not in r for r in independent)
    checks.append({"check": "Both reference versions, all rankings/metrics/plots, all baseline/final top100 extremizers independently covered", "independent_minimal_N": len(ids), "required_extremizers": len(needed)})
    for name in ["results", "historical_target"]:
        d = pd.read_csv(OUT / name / "all_residual_dossiers.csv.gz", float_precision="round_trip")
        assert len(d) == 86296 and len(list((OUT / name / "dossiers").glob("*_top100.json"))) == 56
        v = pd.read_csv(OUT / name / "remaining_residual_variance_components.csv")
        assert v.identity_error.abs().max() < 1e-10
    for name in ["results/stage_metrics.csv", "historical_target/stage_metrics.csv"]:
        s = pd.read_csv(OUT / name)
        assert len(s) == 7*56 and (s.N == 1541).all()
    p = pd.read_csv(OUT / "provenance_analysis/ligand_provenance.csv")
    assert len(p) == 1541 and set(p.provenance_class) == {"F"} and not p.historical_exact_geometry_verified.any()
    checks.append({"check": "Every residual has a dossier; variance identity closes; no fabricated historical exact subset", "dossiers_per_reference": 86296, "historical_exact_N": 0})
    convergence = json.loads((OUT / "convergence/analysis/complete.json").read_text())
    assert convergence["conformer_grid_records"] == 3180 and convergence["ligands"] == 55
    sampling = json.loads((OUT / "campaigns/primary_sampling_3600/acceptance_checks.json").read_text())
    assert sampling["N"] == 31611 and sampling["default_bit_mismatches"] == sampling["L_B5_changes"] == 0 and not sampling["bound_failures"]
    c = OUT / "candidates/angular_sampling"
    replay = json.loads((c / "full_replay_comparison.json").read_text())
    assert replay["passed"] and replay["total_rows"] == 128021 and len(replay["lanes"]) == 13
    for lane in replay["lanes"]:
        assert lane["exit_code"] == 0 and not lane["raw_line_differences"] and lane["before_sha256"] == lane["after_sha256"]
    rust = (c / "final_rust_tests.log").read_text()
    counts = re.findall(r"test result: ok\. (\d+) passed; (\d+) failed; (\d+) ignored", rust)
    assert sum(int(r[0]) for r in counts) == 322 and all(r[1:] == ("0", "0") for r in counts)
    python = (c / "full_python_tests.log").read_text()
    assert "Ran 124 tests" in python and "\nOK\n" in python
    assert "Finished" in (c / "final_clippy.log").read_text() and "Generated" in (c / "rustdoc.log").read_text()
    subprocess.run(["cargo", "fmt", "--check"], cwd=ROOT, check=True, capture_output=True)
    checks.append({"check": "Scientific and generality acceptance", "Rust_tests": 322, "Python_tests": 124, "exact_scientific_records": 128021, "optional_sampling_bound_failures": 0})
    # Save immutable copies of live scripts before sealing the final inventory.
    history = OUT / "scripts_history"
    for f in (OUT / "scripts").glob("*.py"):
        path = history / (sha(f)+".py.txt")
        if not path.exists():
            path.write_bytes(f.read_bytes())
    for plan in OUT.rglob("plan*.json"):
        if "target" in plan.parts:
            continue
        d = json.loads(plan.read_text())
        if d.get("script_sha256"):
            assert (history / (d["script_sha256"]+".py.txt")).exists(), plan
    for doc in ["KRAKEN_RESIDUAL_FORENSICS.md", "REPRODUCE.md", "REQUIREMENTS_AUDIT.md"]:
        text = (OUT / doc).read_text()
        for link in re.findall(r"\]\(([^)]+)\)", text):
            if not link.startswith(("https://", "http://", "#")):
                assert (OUT / link.split("#")[0]).exists(), (doc, link)
    assert len(list((OUT / "family_plots_v2").glob("*.png"))) == 56
    assert json.loads((OUT / "historical_target/reference_precision_v2/complete.json").read_text())["beyond_rounding_ligands"].keys() == {"821", "1036"}
    save(OUT / "verification.json", {"passed": True, "created_utc": datetime.datetime.now(datetime.UTC).isoformat(),
        "script_sha256": sha(__file__), "checks": checks,
        "limits": "Integrity/coverage acceptance, not an assertion that unavailable historical causes were resolved. No R2 threshold used."})
    acceptance = {"accepted": True, "reason": "Independent envelope bound and source-declared sampling convention; default exact replay and scientific tests pass. R2 is not an acceptance criterion.",
        "production_source_sha256": sha(ROOT / "src/geometry/sterimol.rs"), "focused_tests_sha256": sha(ROOT / "tests/sterimol_sampling.rs"),
        "candidate_plan_sha256": sha(c / "plan_before_edit.json"), "scientific_replay_sha256": sha(c / "full_replay_comparison.json"),
        "Rust_tests": 322, "Python_tests": 124, "default_changes": 0}
    save(c / "acceptance.json", acceptance)
    files = {}
    for path in sorted(OUT.rglob("*")):
        rel = path.relative_to(OUT)
        if not path.is_file() or "target" in rel.parts or "__pycache__" in rel.parts or path.name == "BUNDLE_MANIFEST.json":
            continue
        files[str(rel)] = {"bytes": path.stat().st_size, "sha256": sha(path)}
    save(OUT / "BUNDLE_MANIFEST.json", {"created_utc": datetime.datetime.now(datetime.UTC).isoformat(),
        "source_commit": frozen["source_commit"], "files": files,
        "exclusions": ["rebuildable **/target/", "**/__pycache__/", "this self-referential manifest"],
        "source_changes": {"src/geometry/sterimol.rs": sha(ROOT / "src/geometry/sterimol.rs"), "tests/sterimol_sampling.rs": sha(ROOT / "tests/sterimol_sampling.rs")}})
    print(json.dumps({"passed": True, "files": len(files), "bytes": sum(v["bytes"] for v in files.values()), "bundle_manifest_sha256": sha(OUT / "BUNDLE_MANIFEST.json")}))


if __name__ == "__main__":
    main()
