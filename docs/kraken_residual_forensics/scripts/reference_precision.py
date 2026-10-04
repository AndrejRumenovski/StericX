"""Separate decimal serialization differences from reference-version changes."""
import collections
import csv
import decimal
import json
import tarfile
from pathlib import Path

from analyze import FIELDS
from freeze import OUT, save, sha


def main():
    out = OUT / "historical_target/reference_precision_v2"
    out.mkdir(exist_ok=False)
    save(out / "plan.json", {"script_sha256": sha(__file__),
        "model": "Half a unit in each serialized last decimal place; decimal quantization, not scientific uncertainty. No numerical threshold optimized against outputs."})
    D = decimal.Decimal
    hist = {int(r[""]): r for r in csv.DictReader((OUT / "frozen/historical_reference.csv").open())}
    keys = {prop+"_"+red for _, _, prop, _ in FIELDS.values() for red in ["min", "max", "delta", "vburminconf"]}
    records = []
    with tarfile.open(OUT / "frozen/inputs_and_reference.tar.gz") as archive:
        for member in archive:
            if not member.name.startswith("primary/published_dft/"):
                continue
            mid = int(Path(member.name).stem)
            for r in json.load(archive.extractfile(member), parse_float=D, parse_int=D):
                for red in ["min", "max", "delta", "vburminconf"]:
                    key = r["property"]+"_"+red
                    if key not in keys or key not in hist[mid] or r[red] is None:
                        continue
                    try:
                        h, a = D(hist[mid][key]), D(r[red])
                        if not h.is_finite() or not a.is_finite():
                            continue
                    except decimal.InvalidOperation:
                        continue
                    delta = a-h
                    if not delta:
                        continue
                    bound = (D(10)**h.as_tuple().exponent+D(10)**a.as_tuple().exponent)/2
                    records.append({"molecule_id": mid, "property": key, "historical_literal": str(h),
                        "API_literal": str(a), "difference": str(delta), "sum_half_last_decimal_units": str(bound),
                        "beyond_serialization_rounding": abs(delta) > bound})
    with (out / "all_differing_values.csv").open("w") as f:
        writer = csv.DictWriter(f, fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)
    save(out / "complete.json", {"beyond_rounding_ligands": dict(collections.Counter(
        r["molecule_id"] for r in records if r["beyond_serialization_rounding"])),
        "all_differing_rows": len(records), "scope": "All implemented native geometry reference fields; percent reference derives from absolute volume and is not counted twice."})
    print((out / "complete.json").read_text())


if __name__ == "__main__":
    main()
