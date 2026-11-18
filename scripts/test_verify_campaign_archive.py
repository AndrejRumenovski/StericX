"""Tamper checks beside the CI maintenance tool, outside frozen science sources."""

import hashlib
import tempfile
import unittest
from functools import partial
from pathlib import Path, PurePosixPath

from scripts.verify_campaign_archive import (
    published_records,
    verify_archive,
    verify_index,
    verify_records,
)


class CampaignArchiveTests(unittest.TestCase):
    campaign = PurePosixPath("docs/ultimate_engine")

    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.root = Path(self.directory.name)
        self.path = self.root / self.campaign / "recorded.py"
        self.path.parent.mkdir(parents=True)
        self.content = b"value = 1\n"
        self.path.write_bytes(self.content)
        self.records = [
            {
                "path": str(self.campaign / "recorded.py"),
                "bytes": len(self.content),
                "sha256": hashlib.sha256(self.content).hexdigest(),
            }
        ]
        self.verify = partial(verify_records, campaign=self.campaign)

    def test_published_verification_does_not_require_private_measurements(self):
        private = {"path": ".stericx/ultimate_engine/private.json"}
        records = published_records(
            {"records": [*self.records, private]}, self.campaign
        )
        self.assertEqual(records, self.records)
        self.assertEqual(self.verify(self.root, records, ["recorded.py"]), 1)

    def test_same_size_tampering_is_rejected(self):
        self.path.write_bytes(b"value = 2\n")
        with self.assertRaisesRegex(ValueError, "Hash mismatch"):
            self.verify(self.root, self.records, ["recorded.py"])

    def test_missing_published_file_is_rejected(self):
        self.path.unlink()
        with self.assertRaises(FileNotFoundError):
            self.verify(self.root, self.records, ["recorded.py"])

    def test_exclusion_cannot_hide_new_scripts_or_use_a_wildcard(self):
        for exclusions in [
            ["*.py"],
            ["recorded.py", "new.py"],
            ["recorded.py", "recorded.py"],
            [],
        ]:
            with self.subTest(exclusions=exclusions):
                with self.assertRaisesRegex(ValueError, "exactly the sealed"):
                    self.verify(self.root, self.records, exclusions)

    def test_manifest_rejects_duplicate_and_traversing_paths(self):
        for rows in [
            self.records * 2,
            [{"path": str(self.campaign / "../other")}],
            [{"path": "/absolute/path"}],
        ]:
            with self.subTest(rows=rows), self.assertRaises(ValueError):
                published_records({"records": rows}, self.campaign)

    def test_archived_python_must_still_parse(self):
        content = b"def broken(\n"
        self.path.write_bytes(content)
        record = dict(
            self.records[0],
            bytes=len(content),
            sha256=hashlib.sha256(content).hexdigest(),
        )
        with self.assertRaises(SyntaxError):
            self.verify(self.root, [record], ["recorded.py"])

    def test_symlink_cannot_replace_an_archived_file(self):
        target = self.root / "replacement.py"
        target.write_bytes(self.content)
        self.path.unlink()
        self.path.symlink_to(target)
        with self.assertRaisesRegex(ValueError, "replaced with symlink"):
            self.verify(self.root, self.records, ["recorded.py"])

    def test_deleted_index_entry_is_rejected(self):
        with self.assertRaisesRegex(ValueError, "missing from the Git index"):
            verify_index(self.records, set(), ["recorded.py"])

    def test_exclusion_cannot_hide_new_nested_helper(self):
        tracked = {self.records[0]["path"], str(self.campaign / "new/recorded.py")}
        with self.assertRaisesRegex(ValueError, "hide a new helper"):
            verify_index(self.records, tracked, ["recorded.py"])

    def test_new_helper_keeps_repository_lint_rules(self):
        tracked = {self.records[0]["path"], str(self.campaign / "new.py")}
        verify_index(self.records, tracked, ["recorded.py"])

    def test_changed_seal_is_rejected_before_reading_records(self):
        (self.root / self.campaign / "seal.json").write_text("{}\n")
        with self.assertRaisesRegex(ValueError, "Original campaign seal changed"):
            verify_archive(self.root, self.campaign, "seal.json", "0" * 64)


class OccupancyArchiveTests(CampaignArchiveTests):
    campaign = PurePosixPath("docs/occupancy_precision")


if __name__ == "__main__":
    unittest.main()
