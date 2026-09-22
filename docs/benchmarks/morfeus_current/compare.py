"""Untimed complete-corpus comparison. Records residuals; never fits tolerances."""
import csv
import hashlib
import json
import os
from pathlib import Path
import subprocess
from datetime import datetime, timezone

for key in ["OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "BLIS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"]:
    os.environ[key] = "1"

import numpy as np
from morfeus import BuriedVolume, Pyramidalization, Sterimol
from morfeus.io import read_xyz
from morfeus.utils import get_radii

HERE = Path(__file__).resolve().parent
FIELDS = ["sterimol_l", "sterimol_b1", "sterimol_b5", "pyr_p", "pyr_alpha", "buried_volume", "percent_buried_volume", "qvbur_min", "qvbur_max", "max_delta_qvbur", "ovbur_min", "ovbur_max", "near_vbur", "far_vbur"]


def write(p, obj):
    Path(p).write_text(json.dumps(obj, indent=2, allow_nan=False) + "\n")


def capture(name, argv):
    p = subprocess.run(argv, capture_output=True, env=os.environ | {"RAYON_NUM_THREADS":"1"})
    (HERE / "raw" / (name + ".stdout")).write_bytes(p.stdout)
    (HERE / "raw" / (name + ".stderr")).write_bytes(p.stderr)
    write(HERE / "raw" / (name + ".command.json"), dict(argv=argv,returncode=p.returncode, timing=False))
    p.check_returncode()
    return p.stdout.decode()


def reference(path, observed, controlled=False, rotations=3600):
    elements, xyz = read_xyz(path)
    donor = observed["donor"]
    neighbors = observed["neighbors"]
    radii = np.array(get_radii(elements, radii_type="bondi"))
    # A shared explicit radius convention, not reliance on matching table names.
    native_radii = np.array(observed["radii"])
    radius_table_differences = (native_radii - radii).tolist()
    if controlled:
        xyz = np.array(observed["coordinates"])
        center = np.array(observed["center"])
        radii = native_radii
        volume_radii = (native_radii.astype(np.float32) * np.float32(1.17)).astype(float)
    else:
        # Frozen corpus uses matching elements; use explicitly declared native
        # radii to handle any table extension without excluding a molecule.
        radii = np.round(native_radii, 2)
        unit = xyz[neighbors] - xyz[donor]
        unit /= np.linalg.norm(unit, axis=1)[:,None]
        direction = -unit.sum(axis=0)
        assert np.linalg.norm(direction) > 0.01, "Planar center requires separate review"
        center = xyz[donor] + 2.28 * direction / np.linalg.norm(direction)
        volume_radii = radii * 1.17
    sterimol = Sterimol(["H",*elements], np.vstack([center,xyz]),1,donor+2,radii=np.r_[0.,radii],n_rot_vectors=rotations)
    pyr = Pyramidalization(xyz,donor+1,neighbor_indices=[i+1 for i in neighbors])
    frames = []
    # Same center, z toward the virtual metal, x toward each donor substituent.
    # This explicit frame avoids the reference constructor's translated-center
    # aliasing issue; see the pre-existing independent audit adapter.
    for plane in neighbors:
        z = center - xyz[donor]
        z /= np.linalg.norm(z)
        v = xyz[plane] - center
        x = v - np.dot(v,z)*z
        x /= np.linalg.norm(x)
        y = np.cross(z,x)
        aligned = (xyz-center) @ np.array([x,y,z]).T
        b = BuriedVolume(["H",*elements],np.vstack([np.zeros(3),aligned]),1,
                         radii=np.r_[0.,volume_radii],radius=3.5,density=0.01,include_hs=False)
        b.octant_analysis()
        octants = [float(b.octants["buried_volume"][i]) for i in [0,1,2,3,7,6,5,4]]
        quadrants = [octants[i]+octants[i+4] for i in range(4)]
        def populations(points):
            return [int(np.sum((points[:,0]*sx>0)&(points[:,1]*sy>0)&(points[:,2]*sz>0)))
                    for sz in [1,-1] for sx,sy in [(1,1),(-1,1),(-1,-1),(1,-1)]]
        frames.append(dict(plane=plane, buried_volume=float(b.buried_volume), percent_buried_volume=float(b.fraction_buried_volume*100),octants=octants,quadrants=quadrants,near_vbur=sum(octants[4:]),far_vbur=sum(octants[:4]),
                           grid_populations=populations(b._sphere.points),occupied_counts=populations(b._buried_points)))
    q = np.array([f["quadrants"] for f in frames])
    octs = np.array([f["octants"] for f in frames])
    values = dict(sterimol_l=sterimol.L_value, sterimol_b1=sterimol.B_1_value,sterimol_b5=sterimol.B_5_value,
                  pyr_p=pyr.P,pyr_alpha=pyr.alpha, qvbur_min=q.min(),qvbur_max=q.max(),
                  max_delta_qvbur=np.abs(q-np.roll(q,1,axis=1)).max(),ovbur_min=octs.min(),ovbur_max=octs.max())
    values.update({key:np.mean([f[key] for f in frames]) for key in ["buried_volume","percent_buried_volume","near_vbur","far_vbur"]})
    return dict(values={k:float(v) for k,v in values.items()}, frames=frames, center=center.tolist(), radius_table_differences=radius_table_differences,
                coordinate_mode="same f32 represented coordinates/center/radii" if controlled else "same XYZ bytes, float64 reference parsing and center",
                sterimol_rot_vectors=rotations)


def metrics(reference_values, observed_values):
    ref, sut = np.array(reference_values), np.array(observed_values)
    delta = sut-ref
    denom = float(np.sum((ref-ref.mean())**2))
    slope, intercept = np.polyfit(ref,sut,1) if denom else (None,None)
    return dict(N=len(ref), MAE=float(np.mean(abs(delta))),RMSE=float(np.sqrt(np.mean(delta**2))),maximum_absolute_error=float(max(abs(delta))),median_absolute_error=float(np.median(abs(delta))),R_squared=1-float(np.sum(delta**2))/denom if denom else None,slope=slope,intercept=intercept)


def main():
    os.sched_setaffinity(0,{2})
    manifest = json.loads((HERE / "input_manifest.json").read_text())
    paths = [str(HERE / r["filename"]) for r in manifest["corpus"]]
    for row,path in zip(manifest["corpus"],paths,strict=True):
        assert hashlib.sha256(Path(path).read_bytes()).hexdigest() == row["sha256"]
    executable = HERE / "build/stericx"
    freeze = json.loads((HERE / "freeze.json").read_text())
    assert hashlib.sha256(executable.read_bytes()).hexdigest() == freeze["executable_sha256"]
    cli = json.loads(capture("stericx-cli",[str(executable),"descriptors",*paths,"--sterimol-axis","coordination","--format","json"]))
    observations = [json.loads(line) for line in capture("stericx-observer",[str(HERE / "adapter/target/release/morfeus-current-observer"),*paths]).splitlines()]
    assert len(cli) == len(observations) == 56
    for a,b in zip(cli,observations,strict=True):
        assert a["file"] == b["file"]
        for field in set(FIELDS) & set(a):
            assert np.float32(a[field]).tobytes() == np.float32(b[field]).tobytes(), (field,a[field],b[field])
    references = {"source_float64":[], "controlled_f32_inputs":[], "matched_angular_count_361":[]}
    for path,obs in zip(paths,observations,strict=True):
        references["source_float64"].append(reference(path,obs))
        references["controlled_f32_inputs"].append(reference(path,obs,controlled=True))
        references["matched_angular_count_361"].append(reference(path,obs,controlled=True,rotations=361))
    write(HERE / "raw/morfeus-reference.json",references)
    write(HERE / "geometry_configuration.json",[dict(filename=r["filename"], donor=o["donor"],neighbors=o["neighbors"],center=o["center"],radii=o["radii"],elements=o["elements"],precision="f32 native inputs; explicit native center/radii for controlled comparison") for r,o in zip(manifest["corpus"],observations,strict=True)])
    all_metrics = {}
    with (HERE / "scientific_comparisons.csv").open("w",newline="") as f:
        out = csv.DictWriter(f,fieldnames=["lane","filename","descriptor","stericx","morfeus","absolute_difference","relative_difference"])
        out.writeheader()
        for lane,rows in references.items():
            all_metrics[lane] = {}
            for field in FIELDS:
                a,b = [r["values"][field] for r in rows],[o[field] for o in observations]
                all_metrics[lane][field] = metrics(a,b)
                for src,x,y in zip(manifest["corpus"],a,b,strict=True):
                    out.writerow(dict(lane=lane,filename=src["filename"],descriptor=field,stericx=y,morfeus=x,absolute_difference=abs(y-x),relative_difference=abs((y-x)/x) if x else ""))
    write(HERE / "scientific_summary.json",dict(created_utc=datetime.now(timezone.utc).isoformat(),metrics=all_metrics,cli_observer_bitwise_agreement=True,scientific_gate="Residual report only; separate scientific_gate.json determines admission by exact integer counts and bitwise rounding reconstruction. Existing audit maxima are not acceptance tolerances.",tolerances_invented=False,all_56_retained=True))


if __name__ == "__main__":
    main()
