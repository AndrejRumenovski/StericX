"""The current observer may add a direct edge, never change pinned crates."""
from pathlib import Path
import sys
import tomllib
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from replay_scientific_remediation import AUDIT, REPO, bind_observer_dependencies


class ObserverDependencyBinding(unittest.TestCase):
    def test_historical_source_keeps_the_original_template_exactly(self):
        cargo = (AUDIT / "observer/Cargo.toml").read_text()
        lock = (AUDIT / "observer/Cargo.lock").read_text()
        result = bind_observer_dependencies(cargo, lock, (REPO / "Cargo.lock").read_text(), b"no direct derives")
        self.assertEqual(result, (cargo, lock, []))

    def test_current_derives_add_only_the_direct_serde_edge(self):
        cargo = (AUDIT / "observer/Cargo.toml").read_text()
        lock = (AUDIT / "observer/Cargo.lock").read_text()
        changed, changed_lock, adjustments = bind_observer_dependencies(cargo, lock, (REPO / "Cargo.lock").read_text(), b"serde::Serialize")
        old = tomllib.loads(lock); new = tomllib.loads(changed_lock)
        for before, after in zip(old["package"], new["package"], strict=True):
            if before["name"] == "stericx-audit-observer":
                self.assertEqual(set(after["dependencies"]) - set(before["dependencies"]), {"serde"})
                after["dependencies"] = before["dependencies"]
            self.assertEqual(before, after)
        self.assertEqual(old, new)
        source_serde = next(p for p in tomllib.loads((REPO / "Cargo.lock").read_text())["package"] if p["name"] == "serde")
        self.assertEqual(tomllib.loads(changed)["dependencies"]["serde"]["version"], "=" + source_serde["version"])
        self.assertEqual(adjustments[0]["version"], source_serde["version"])


if __name__ == "__main__":
    unittest.main()
