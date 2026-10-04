"""Fresh replay of the validated scientific/invariance lanes, without mutation."""
import datetime
import gzip
import json
import subprocess

from freeze import FINAL, OUT, ROOT, save, sha


def main():
    out = OUT / "validation/generality"
    out.mkdir(exist_ok=False)
    manifest_path = FINAL / "observations/manifest.json"
    manifest = json.loads(manifest_path.read_text())
    names = [n for n in manifest["lanes"] if not n.startswith("kraken")]
    save(out / "plan.json", {"created_utc": datetime.datetime.now(datetime.UTC).isoformat(),
        "script_sha256": sha(__file__), "source_manifest_sha256": sha(manifest_path), "lanes": names,
        "scope": "all non-Kraken validated native observer lanes; default Kraken cohort replay has its separate full 31611-row receipt; no tolerance or target rescaling"})
    results = []
    for name in names:
        lane = manifest["lanes"][name]
        requests = lane["requests"]
        expected = lane["stdout"]
        assert sha(requests["path"]) == requests["sha256"]
        assert sha(expected["path"]) == expected["sha256"]
        with open(requests["path"]) as inp, open(expected["path"]) as oracle, gzip.open(out / f"{name}.jsonl.gz", "wt") as dst, open(out / f"{name}.stderr", "w") as err:
            p = subprocess.Popen([str(OUT / "frozen/current_observer")], stdin=inp, stdout=subprocess.PIPE, stderr=err, text=True)
            differences = []
            count = 0
            for actual, old in zip(p.stdout, oracle, strict=True):
                row = json.loads(actual)
                reference = json.loads(old)
                count += 1
                dst.write(actual)
                if row != reference:
                    differences.append({"row": count, "id": row.get("id")})
            code = p.wait()
        results.append({"lane": name, "rows": count, "exit_code": code, "exact_differences": differences,
            "request_sha256": requests["sha256"], "output_sha256": sha(out / f"{name}.jsonl.gz")})
        print(name, count, code, len(differences), flush=True)
    frozen = json.loads((OUT / "frozen/manifest.json").read_text())
    source = {name: sha(ROOT / name) == row["sha256"] for name, row in frozen["source_files"].items() if name.startswith("src/") or name in ("Cargo.toml", "Cargo.lock")}
    save(out / "complete.json", {"lanes": results, "production_source_unchanged": source,
        "passed": all(v["exit_code"] == 0 and not v["exact_differences"] for v in results) and all(source.values())})


if __name__ == "__main__":
    main()
