"""Verify published evidence, sample completeness, hashes and exact count gate.

Standard library only. This is integrity verification, not a scientific rerun.
"""
import argparse
import gzip
import hashlib
import json
from pathlib import Path
import struct

HERE=Path(__file__).resolve().parent


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):return json.loads(Path(p).read_text())


def inventory():
    records=[]
    for p in sorted(HERE.rglob("*")):
        if not p.is_file():continue
        rel=p.relative_to(HERE)
        if any(part in ["target","venv","__pycache__","repeated_10000"] for part in rel.parts):continue
        if rel.as_posix() in ["evidence_manifest.json","evidence_manifest.json.sha256","verification.json"]:continue
        records.append(dict(path=rel.as_posix(),bytes=p.stat().st_size,sha256=sha(p)))
    return records


def verify():
    manifest=read(HERE/"evidence_manifest.json")
    assert sha(HERE/"evidence_manifest.json")== (HERE/"evidence_manifest.json.sha256").read_text().strip()
    for r in manifest["files"]:
        p=HERE/r["path"]
        assert p.stat().st_size==r["bytes"] and sha(p)==r["sha256"],r["path"]
    assert manifest["files"]==inventory(),"Uninventoried or missing published file"
    assert sha(HERE/"build/stericx")==read(HERE/"freeze.json")["executable_sha256"]
    for receipt in [read(HERE/"preflight.json"),read(HERE/"pyramidalization/scientific_gate.json")]:
        for path,digest in receipt["identities"].items():
            archived=HERE/path
            if path.startswith("adapter/target/release/"):archived=HERE/"build"/Path(path).name
            assert sha(archived)==digest,path
    corpus=read(HERE/"input_manifest.json")
    assert len(corpus["corpus"])==56 and len(corpus["repeated_10000"])==10000
    for row in corpus["corpus"]:assert sha(HERE/row["filename"])==row["sha256"]
    for i,row in enumerate(corpus["repeated_10000"]):assert row["corpus_index"]==i%56 and row["sha256"]==corpus["corpus"][i%56]["sha256"]
    gate=read(HERE/"scientific_gate.json")
    assert gate["passed"] and gate["numerical_tolerance"] is None and gate["failures"]==[]
    evidence=read(HERE/"raw/exact_volume_evidence.json")
    assert len(evidence)==56
    for row in evidence:
        assert not row["issues"]
        for a,b in zip(row["native"]["frames"],row["morfeus"]["frames"],strict=True):
            assert a["occupied_counts"]==b["occupied_counts"] and a["grid_populations"]==b["grid_populations"]
    counts={}
    for experiment,rawpath in [("volume",HERE/"raw_timings.json"),("pyramidalization",HERE/"pyramidalization/raw_timings.json")]:
        raw=read(rawpath)
        assert raw["status"]=="complete" and raw["exclusions"]==[] and len(raw["samples"])==120
        expected={(n,w,t,i) for n in [56,1000,10000] for w in [1,2,4,6] for t in ["stericx","morfeus"] for i in [-1,0,1,2,3]}
        seen=set()
        for s in raw["samples"]:
            key=(s["structures"],s["workers"],s["tool"],s["pair"])
            assert key not in seen;seen.add(key)
            assert s["status"]=="passed" and s["metrics"]["returncode"]==0
            assert s["metrics"]["wall_ns"]>0
            prefix=HERE/s["raw_prefix"]
            contents=gzip.decompress(prefix.with_suffix(".stdout.gz").read_bytes())
            assert hashlib.sha256(contents).hexdigest()==s["stdout_sha256"]
            assert sha(prefix.with_suffix(".stderr"))==s["stderr_sha256"]
            assert read(prefix.with_suffix(".metrics.json"))==s["metrics"]
            assert len(json.loads(contents))==s["structures"]
        assert seen==expected
        counts[experiment]=dict(launches=len(seen),measured=96,warmups=24)
    py=read(HERE/"pyramidalization/scientific_gate.json")
    assert py["passed"] and py["numerical_tolerance"] is None
    native=read(HERE/"pyramidalization/preflight-stericx-w1.stdout")
    ref=read(HERE/"pyramidalization/preflight-morfeus-w1.stdout")
    for a,b in zip(native,ref,strict=True):
        for k in ["pyr_p","pyr_alpha"]:assert struct.pack("f",a[k])==struct.pack("f",b[k])
    memory=read(HERE/"memory_results.json");assert memory["complete"] and len(memory["results"])==8
    result=dict(passed=True,published_files=len(manifest["files"]),samples=counts,memory_runs=8,scope="Historical evidence integrity, original regional-count/rounding checks and launch completeness only; scientific admission withdrawn; not a prospective scientific gate or benchmark rerun")
    (HERE/"verification.json").write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps(result,indent=2))


if __name__=="__main__":
    parser=argparse.ArgumentParser();parser.add_argument("--seal",action="store_true");args=parser.parse_args()
    if args.seal:
        p=HERE/"evidence_manifest.json"
        assert not p.exists(),"Refuse to replace an existing evidence seal"
        p.write_text(json.dumps(dict(files=inventory(),excluded="Disposable target/venv/cache and regenerated repeated_10000; original 56 input files and manifests retained"),indent=2)+"\n")
        p.with_suffix(".json.sha256").write_text(sha(p)+"\n")
    verify()
