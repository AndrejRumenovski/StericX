"""Tamper checks for the CI gate protecting archived campaign evidence."""

import hashlib
import tempfile
import unittest
from pathlib import Path

from scripts.verify_campaign_archive import (
    published_records,
    verify_records,
)


class CampaignArchiveTests(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.root = Path(self.directory.name)
        self.path = self.root / "docs/ultimate_engine/recorded.py"
        self.path.parent.mkdir(parents=True)
        self.content = b"value = 1\n"
        self.path.write_bytes(self.content)
        self.records = [
            {
                "path": "docs/ultimate_engine/recorded.py",
                "bytes": len(self.content),
                "sha256": hashlib.sha256(self.content).hexdigest(),
            }
        ]

    def test_published_verification_does_not_require_private_measurements(self):
        private = {"path": ".stericx/ultimate_engine/private.json"}
        records = published_records({"records": [*self.records, private]})
        self.assertEqual(records, self.records)
        self.assertEqual(verify_records(self.root, records, ["recorded.py"]), 1)

    def test_same_size_tampering_is_rejected(self):
        self.path.write_bytes(b"value = 2\n")
        with self.assertRaisesRegex(ValueError, "Hash mismatch"):
            verify_records(self.root, self.records, ["recorded.py"])

    def test_missing_published_file_is_rejected(self):
        self.path.unlink()
        with self.assertRaises(FileNotFoundError):
            verify_records(self.root, self.records, ["recorded.py"])

    def test_exclusion_cannot_hide_new_scripts_or_use_a_wildcard(self):
        for exclusions in [["*.py"], ["recorded.py", "new.py"], []]:
            with self.subTest(exclusions=exclusions):
                with self.assertRaisesRegex(ValueError, "exactly the sealed"):
                    verify_records(self.root, self.records, exclusions)

    def test_manifest_rejects_duplicate_and_traversing_paths(self):
        for rows in [self.records * 2, [{"path": "docs/ultimate_engine/../other"}]]:
            with self.subTest(rows=rows), self.assertRaises(ValueError):
                published_records({"records": rows})

    def test_archived_python_must_still_parse(self):
        content = b"def broken(\n"
        self.path.write_bytes(content)
        record = dict(
            self.records[0],
            bytes=len(content),
            sha256=hashlib.sha256(content).hexdigest(),
        )
        with self.assertRaises(SyntaxError):
            verify_records(self.root, [record], ["recorded.py"])


if __name__ == "__main__":
    unittest.main()
