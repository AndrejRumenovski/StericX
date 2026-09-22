"""Postflight source, dependency, topology, radius and measurement review."""
import importlib.metadata
import json
from pathlib import Path
import shutil
import tomllib
from benchmark import HERE,ROOT,sha,read,write,verify_freeze
import morfeus
import morfeus_batch


def main():
    verify_freeze()
    assert read(HERE/"raw_timings.json")["status"]=="complete"
    reference_sources=read(HERE/"environment/reference-sources.stdout")
    package=Path(morfeus.__file__).parent
    for name,digest in reference_sources.items():assert sha(package/name)==digest,name
    packages={}
    for line in (HERE/"environment/requirements.txt").read_text().splitlines():
        name,version=line.split("==")
        assert importlib.metadata.version(name)==version,name
        packages[name]=version
    root_lock=tomllib.loads((ROOT/"Cargo.lock").read_text())
    adapter_lock=tomllib.loads((HERE/"adapter/Cargo.lock").read_text())
    base={(p["name"],p["version"]):p.get("checksum") for p in root_lock["package"]}
    adapted={(p["name"],p["version"]):p.get("checksum") for p in adapter_lock["package"] if p["name"]!="morfeus-current-observer"}
    assert base==adapted,"Adapter changed dependency versions"
    observations=[json.loads(x) for x in (HERE/"raw/stericx-observer.stdout").read_text().splitlines()]
    refs=read(HERE/"raw/morfeus-reference.json")["source_float64"]
    topology=[]
    for o,r in zip(observations,refs,strict=True):
        elements,xyz,d,ns,center,radii=morfeus_batch.geometry(o["file"])
        assert d==o["donor"] and ns==o["neighbors"]
        assert center.tolist()==r["center"]
        assert list(elements)==o["elements"]
        topology.append(dict(filename=str(Path(o["file"]).relative_to(HERE)),donor=d,neighbors=ns,center=center.tolist(),scaled_volume_radii=radii.tolist(),passed=True))
    preflight=read(HERE/"preflight.json")
    for path,digest in preflight["identities"].items():assert sha(HERE/path)==digest,path
    for name in ["volume-batch","morfeus-current-observer","grid-observer","pyramid-batch"]:
        src=HERE/"adapter/target/release"/name
        if src.exists():shutil.copy2(src,HERE/"build"/name)
    s=read(HERE/"scientific_summary.json")
    s["scientific_gate"]="Residual report only; scientific_gate.json admits volume by exact integer equality and bitwise rounding reconstruction. No historical error maximum was used as a tolerance."
    write(HERE/"scientific_summary.json",s)
    write(HERE/"protocol_review.json",dict(passed=True,production_sources_unchanged=True,reference_python_sources_unchanged=True,packages=packages,adapter_dependencies_identical=True,topology=topology,
        timing_scope="Fresh XYZ parsing, sole-P donor and donor-neighbor inference, independent unit-vector virtual center, same 3 oriented buried-volume calls and all 9 outputs; native current public API and reference public BuriedVolume API. No timing of visibility observer.",
        wrapper_notes="Current shipping descriptors CLI includes Sterimol and pyramidalization. Timing volume-only API drivers is explicitly scoped; no full shipping-CLI speedup claim."))


if __name__=="__main__":main()
