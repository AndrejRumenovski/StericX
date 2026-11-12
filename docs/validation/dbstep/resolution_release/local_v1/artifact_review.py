"""Check RC archive records and source identity without executing artifacts."""

import base64
import csv
import hashlib
import io
import json
import pathlib
import tarfile
import tomllib
import zipfile

root = pathlib.Path.cwd()
local = root / ".stericx/dbstep_resolution/release_candidate_v1/local"
receipt = json.loads((local / "receipt.json").read_text())
wheel = pathlib.Path(receipt["artifacts"]["wheel"]["path"])
sdist = pathlib.Path(receipt["artifacts"]["sdist"]["path"])
checks = []


def digest(data):
    return hashlib.sha256(data).hexdigest()


with zipfile.ZipFile(wheel) as z:
    names = z.namelist()
    assert len(names) == len(set(names))
    record_name = "stericx-0.4.0rc1.dist-info/RECORD"
    rows = list(csv.reader(io.StringIO(z.read(record_name).decode())))
    assert set(r[0] for r in rows) == set(names)
    for name, encoded, length in rows:
        data = z.read(name)
        if name == record_name:
            assert not encoded and not length
        else:
            expected = "sha256=" + base64.urlsafe_b64encode(
                hashlib.sha256(data).digest()
            ).decode().rstrip("=")
            assert encoded == expected and int(length) == len(data)
        checks.append(
            {
                "archive": "wheel",
                "member": name,
                "sha256": digest(data),
                "record_verified": True,
            }
        )
    for name in ["__init__.py", "_types.py", "py.typed"]:
        assert (
            z.read("stericx/" + name) == (root / "python/stericx" / name).read_bytes()
        )
    metadata = z.read("stericx-0.4.0rc1.dist-info/METADATA").decode()
    assert "\nVersion: 0.4.0rc1\n" in metadata

with tarfile.open(sdist) as t:
    prefix = "stericx-0.4.0rc1/"
    files = {
        m.name.removeprefix(prefix): t.extractfile(m).read()
        for m in t.getmembers()
        if m.isfile()
    }
    required = ["Cargo.toml", "Cargo.lock", "python/Cargo.lock", "python/rust/lib.rs"]
    required += [p.relative_to(root).as_posix() for p in (root / "src").rglob("*.rs")]
    required += [
        p.relative_to(root).as_posix()
        for p in (root / "python/tests").rglob("*")
        if p.is_file()
        and "__pycache__" not in p.parts
        and ".pytest_cache" not in p.parts
    ]
    for name in required:
        assert files[name] == (root / name).read_bytes(), name
        checks.append(
            {
                "archive": "sdist",
                "member": name,
                "sha256": digest(files[name]),
                "source_exact": True,
            }
        )
    for name in ["__init__.py", "_types.py", "py.typed"]:
        assert files["stericx/" + name] == (root / "python/stericx" / name).read_bytes()
        checks.append(
            {
                "archive": "sdist",
                "member": "stericx/" + name,
                "sha256": digest(files["stericx/" + name]),
                "source_exact": True,
            }
        )
    manifest = tomllib.loads(files["python/Cargo.toml"].decode())
    expected = tomllib.loads((root / "python/Cargo.toml").read_text())
    expected["package"]["readme"] = "README.md"
    assert manifest == expected
    assert tomllib.loads(files["Cargo.toml"].decode())["package"]["version"] == "0.3.0"
    project = tomllib.loads(files["pyproject.toml"].decode())
    expected = tomllib.loads((root / "python/pyproject.toml").read_text())
    expected["tool"]["maturin"]["manifest-path"] = "python/Cargo.toml"
    assert project == expected
    assert project["project"]["version"] == "0.4.0rc1"
result = {
    "passed": True,
    "wheel_record_entries": len(rows),
    "sdist_exact_source_members": len(required) + 3,
    "checks": checks,
    "documented_packaging_transformations": [
        "python/Cargo.toml adds package.readme = README.md",
        "root pyproject.toml is python/pyproject.toml with "
        "tool.maturin.manifest-path = python/Cargo.toml",
    ],
    "artifacts": receipt["artifacts"],
    "controller_sha256": digest(pathlib.Path(__file__).read_bytes()),
}
(local / "artifact_review.json").write_text(json.dumps(result, indent=2) + "\n")
print(
    json.dumps(
        {
            k: result[k]
            for k in ["passed", "wheel_record_entries", "sdist_exact_source_members"]
        }
    )
)
