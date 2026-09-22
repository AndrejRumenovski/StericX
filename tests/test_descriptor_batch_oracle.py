"""Descriptor batches must preserve values, duplicates, order, and error policy."""

from __future__ import annotations

import copy
import csv
import io
import json
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

import check_scientific_cli_equivalence as oracle


class DescriptorBatchOracleTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        source = oracle.AUDIT / "geometry/focused/cli_outputs.json"
        cls.audited = [
            {
                "id": f"geometry_{row['id']}",
                "argv": row["argv"][1:],
                "expected_status": row["returncode"],
                "historical_stdout": row["stdout"],
                "historical_stderr": row["stderr"],
            }
            for row in json.loads(source.read_text())
        ]
        cls.cases = {
            case["id"]: case for case in oracle.descriptor_batch_cases(cls.audited)
        }

    def test_fixed_coverage_has_each_failure_policy_in_each_format(self):
        self.assertEqual(
            set(self.cases),
            {
                f"batch_{kind}_{output_format}"
                for kind in ("mixed", "all_failures", "donor_index_precheck")
                for output_format in ("json", "text", "csv")
            },
        )
        for case in self.cases.values():
            self.assertEqual(case["argv"][0], "descriptors")
            self.assertEqual(case["argv"].count("--format"), 1)

    def test_success_order_duplicates_and_audited_numeric_spellings(self):
        case = self.cases["batch_mixed_json"]
        records = json.loads(case["historical_stdout"])
        self.assertEqual(
            [Path(record["file"]).stem for record in records],
            [
                "triphenylphosphine",
                "methylphosphine",
                "triphenylphosphine",
                "asymmetric_phosphine",
            ],
        )
        self.assertEqual(records[0], records[2])
        originals = {row["id"]: row for row in self.audited}
        for record in records:
            name = f"geometry_{Path(record['file']).stem}"
            self.assertEqual(
                record, json.loads(originals[name]["historical_stdout"])[0]
            )
            inner = originals[name]["historical_stdout"][2:-3]
            self.assertIn(inner, case["historical_stdout"])
        self.assertEqual(case["expected_status"], 0)

    def test_text_and_csv_have_four_ordered_records_without_changed_values(self):
        text = self.cases["batch_mixed_text"]["historical_stdout"]
        blocks = text.rstrip("\n").split("\n\n")
        self.assertEqual(len(blocks), 4)
        self.assertEqual(blocks[0], blocks[2])
        self.assertIn("  Sterimol      L 6.93   B1 2.49   B5 6.76   Å\n", blocks[0])
        self.assertIn("  pyramidalization pyr_P 0.891   pyr_alpha 23.29°", blocks[0])
        rows = list(
            csv.DictReader(
                io.StringIO(self.cases["batch_mixed_csv"]["historical_stdout"])
            )
        )
        self.assertEqual(len(rows), 4)
        self.assertEqual(len(rows[0]), 16)
        self.assertEqual(rows[0], rows[2])
        self.assertEqual(rows[0]["sterimol_l"], "6.9254")
        self.assertEqual(rows[0]["pyr_p"], "0.8909")
        self.assertEqual(
            [Path(row["file"]).stem for row in rows],
            [Path(block.splitlines()[0]).stem for block in blocks],
        )

    def test_text_csv_round_original_f32_not_short_decimal_approximation(self):
        row = copy.deepcopy(
            next(
                row
                for row in self.audited
                if row["id"] == "geometry_triphenylphosphine"
            )
        )
        values = json.loads(row["historical_stdout"])
        values[0]["sterimol_l"] = 1.23445
        row["historical_stdout"] = json.dumps(values)
        output = oracle.descriptor_batch_stdout([row], "csv")
        parsed = next(csv.DictReader(io.StringIO(output)))
        self.assertEqual(parsed["sterimol_l"], "1.2344")
        self.assertNotEqual(parsed["sterimol_l"], f"{1.23445:.4f}")

    def test_failure_order_duplicates_terminal_error_and_precheck(self):
        for output_format in ("json", "text", "csv"):
            case = self.cases[f"batch_mixed_{output_format}"]
            lines = case["historical_stderr"].splitlines()
            self.assertIn("nearly_collinear.xyz:", lines[0])
            self.assertIn("multiple_phosphorus.xyz:", lines[1])
            self.assertEqual(lines[1], lines[2])
            self.assertEqual(lines[3], "featurized 4 of 7 files (3 skipped)")
            case = self.cases[f"batch_all_failures_{output_format}"]
            self.assertEqual(case["expected_status"], 2)
            self.assertEqual(case["historical_stdout"], "")
            lines = case["historical_stderr"].splitlines()
            self.assertEqual(len(lines), 4)
            self.assertEqual(lines[0], lines[1])
            self.assertIn("multiple_phosphorus.xyz:", lines[2])
            self.assertEqual(lines[3], "error: no ligand files could be featurized")
            case = self.cases[f"batch_donor_index_precheck_{output_format}"]
            self.assertEqual(case["expected_status"], 2)
            self.assertEqual(case["historical_stdout"], "")
            self.assertEqual(
                case["historical_stderr"],
                "error: --donor-index applies to a single file; "
                "omit it for batch runs\n",
            )
            self.assertTrue(case["batch_contract"]["precheck_before_file_processing"])

    def test_anchored_contract_rejects_equal_but_wrong_batch_observations(self):
        case = {
            **self.cases["batch_mixed_json"],
            "outputs": [],
            "expected_artifacts": {},
        }
        records = json.loads(case["historical_stdout"])
        stderr = case["historical_stderr"].splitlines(keepends=True)
        changes = (
            ("stdout", json.dumps(list(reversed(records)))),
            ("stdout", json.dumps(records[:2] + records[3:])),
            ("stdout", case["historical_stdout"].replace("6.9253626", "6.925363")),
            ("stderr", "".join([stderr[1], stderr[0], *stderr[2:]])),
            ("stderr", ""),
        )
        for stream, wrong in changes:
            with self.subTest(stream=stream), tempfile.TemporaryDirectory() as name:
                directory = Path(name)
                (directory / "result.json").write_text('{"returncode": 0}')
                for key in ("stdout", "stderr"):
                    (directory / key).write_text(case[f"historical_{key}"])
                oracle.validate_outputs(case, directory)
                (directory / stream).write_text(wrong)
                with self.assertRaisesRegex(
                    ValueError, f"batch {stream} contract changed"
                ):
                    oracle.validate_outputs(case, directory)

    def test_derivation_rejects_changed_audited_options(self):
        audited = copy.deepcopy(self.audited)
        row = next(row for row in audited if row["id"] == "geometry_triphenylphosphine")
        row["argv"].extend(["--sterimol-axis", "coordination"])
        with self.assertRaisesRegex(ValueError, "audited descriptor options changed"):
            oracle.descriptor_batch_cases(audited)


if __name__ == "__main__":
    unittest.main()
