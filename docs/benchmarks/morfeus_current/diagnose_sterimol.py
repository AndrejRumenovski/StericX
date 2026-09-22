"""Explain excluded Sterimol residuals with the unchanged independent audit math."""
import importlib.util
import json
import os
from pathlib import Path
for key in ["OPENBLAS_NUM_THREADS","OMP_NUM_THREADS","MKL_NUM_THREADS"]:os.environ[key]="1"
import numpy as np
from compare import HERE,write,metrics
from freeze import sha


def main():
    assert json.loads((HERE/"raw_timings.json").read_text())["status"]=="complete","Do not compete with native timings"
    source=HERE.parents[2]/"docs/scientific_accuracy_audit/scripts/geometry_reference.py"
    spec=importlib.util.spec_from_file_location("independent_audit_reference",source)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    obs=[json.loads(x) for x in (HERE/"raw/stericx-observer.stdout").read_text().splitlines()]
    refs=json.loads((HERE/"raw/morfeus-reference.json").read_text())["controlled_f32_inputs"]
    rows=[]
    for o,m in zip(obs,refs,strict=True):
        exact=module.sterimol_exact(np.array(o["coordinates"]),np.array(o["radii"]),np.array(o["center"]),o["donor"])
        exact["l"]+=0.4
        rows.append(dict(filename=str(Path(o["file"]).relative_to(HERE)),analytic=exact,stericx={k:o["sterimol_"+k] for k in exact},morfeus={k:m["values"]["sterimol_"+k] for k in exact},
                         signed_stericx_error={k:o["sterimol_"+k]-v for k,v in exact.items()},signed_morfeus_error={k:m["values"]["sterimol_"+k]-v for k,v in exact.items()}))
    stats={tool:{k:metrics([r["analytic"][k] for r in rows],[r[tool][k] for r in rows]) for k in ["l","b1","b5"]} for tool in ["stericx","morfeus"]}
    write(HERE/"sterimol_diagnosis.json",dict(reference_source=str(source.relative_to(HERE.parents[2])),reference_sha256=sha(source),rows=rows,metrics=stats,
            conclusion="No timed Sterimol claim. The independent continuous support-envelope minimum distinguishes angular sampling residuals from input/center/radius conventions. Native B1 scans 360 directions in a shortest-arc Z frame; morfeus scans 3600 endpoint-inclusive directions in a Kabsch X frame. The separate 361-direction diagnostic matches one-degree spacing but not transverse phase. Native B5 is a radial maximum; morfeus B5 is sampled angular support. No phase tuning, program changes, or fitted tolerances."))


if __name__=="__main__":main()
