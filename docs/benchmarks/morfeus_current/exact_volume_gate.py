if __name__ == "__main__":
    raise SystemExit("Historical exploratory runner retired. Use reproduce.py and the frozen EQUIVALENCE_POLICY.md.")
"""Require exact discrete populations, then account for every f32 rounding step.

No fitted tolerance or historical maximum is used. The volume estimator and
octant conventions are those already specified in the independent audit METHODS.
An equality failure stops this gate; it is never converted to an allowed error.
"""
import hashlib
import json
import os
from pathlib import Path
import numpy as np
from compare import HERE, capture, reference, write

F = np.float32


def bits(x):
    return F(x).tobytes()


def total(values):
    result = F(0)
    for v in values:
        result = F(result+F(v))
    return result


def mean(values):
    return total(F(v)/F(3) for v in sorted(values))


def reconstruct(frames, volume):
    results = []
    for frame in frames:
        occupied = frame["occupied_counts"]
        population = frame["grid_populations"]
        octants = [F(F(F(n)/F(d))*volume)/F(8) for n,d in zip(occupied,population,strict=True)]
        quadrants = [F(F(F(occupied[i]+occupied[i+4])/F(population[i]+population[i+4]))*volume)/F(4) for i in range(4)]
        results.append(dict(buried_volume=F(F(sum(occupied))/F(sum(population)))*volume,
                            near_vbur=total(octants[4:]),far_vbur=total(octants[:4]),octants=octants,quadrants=quadrants))
    vals = {k:mean([r[k] for r in results]) for k in ["buried_volume","near_vbur","far_vbur"]}
    vals.update(percent_buried_volume=F(100)*F(vals["buried_volume"]/volume),
                qvbur_min=min(x for r in results for x in r["quadrants"]),
                qvbur_max=max(x for r in results for x in r["quadrants"]),
                ovbur_min=min(x for r in results for x in r["octants"]),
                ovbur_max=max(x for r in results for x in r["octants"]),
                max_delta_qvbur=max(abs(F(r["quadrants"][i]-r["quadrants"][(i+3)%4])) for r in results for i in range(4)))
    return results,vals


def main():
    os.sched_setaffinity(0,{2})
    root=HERE.parents[2]
    source=(root/"src/geometry/buried_volume.rs").read_bytes()
    assert (HERE/"adapter/src/copied_buried_volume.rs").read_bytes()==source+b"\n"+(HERE/"adapter/grid_append.rs").read_bytes()
    obs=[json.loads(x) for x in (HERE/"raw/stericx-observer.stdout").read_text().splitlines()]
    paths=[r["file"] for r in obs]
    native=[json.loads(x) for x in capture("native-grid-observer",[str(HERE/"adapter/target/release/grid-observer"),*paths]).splitlines()]
    assert len(native)==len(obs)==56
    records=[]
    failures=[]
    for path,o,s in zip(paths,obs,native,strict=True):
        m=reference(path,o)
        frames,values=reconstruct(m["frames"],F(s["sphere_volume_f32"]))
        issues=[]
        for a,b,c in zip(s["frames"],m["frames"],frames,strict=True):
            for key in ["grid_populations","occupied_counts"]:
                if a[key]!=b[key]: issues.append(dict(plane=a["plane"],field=key,native=a[key],morfeus=b[key]))
            for key in ["buried_volume","near_vbur","far_vbur"]:
                if bits(a[key])!=bits(c[key]): issues.append(dict(plane=a["plane"],field=key))
            for key in ["quadrants","octants"]:
                if any(bits(x)!=bits(y) for x,y in zip(a[key],c[key],strict=True)): issues.append(dict(plane=a["plane"],field=key))
        for key,value in values.items():
            if bits(value)!=bits(o[key]): issues.append(dict(field=key,native=o[key],reconstructed=float(value)))
        records.append(dict(filename=str(Path(path).relative_to(HERE)),morfeus=m,native=s,issues=issues))
        failures.extend(dict(filename=path,issue=i) for i in issues)
    write(HERE/"raw/exact_volume_evidence.json",records)
    gate=dict(passed=not failures,scope="56 unique conformers; three planes each; all 8 octant populations and occupied counts exactly equal; each native f32 plane and aggregate output exactly reconstructed from independent morfeus counts",
              numerical_tolerance=None,reason="Exact integer equality and bitwise reconstruction; stricter than a residual tolerance, no historical maximum promoted to a bound",
              scientific_definition="Independent audit geometry/METHODS.md buried-volume equations and finite region normalization; corrected three-plane means in remediation geometry/RESULTS.md",
              descriptors=list(values),excluded={"Sterimol":"Angular grids and transverse phase differ; no existing authorized error threshold for the observed scan residual. Excluded from headline and combined timing.","pyramidalization":"Compared and retained, not part of volume-only timing.","combined":"Not timed because Sterimol equivalence at matching convergence is not established."},
              structures=56,planes=168,failures=failures,source_sha256=hashlib.sha256(source).hexdigest())
    write(HERE/"scientific_gate.json",gate)
    if failures: raise RuntimeError(f"Exact gate failed: {len(failures)} issues")


if __name__=="__main__": main()
