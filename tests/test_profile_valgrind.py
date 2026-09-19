"""Validate profiler units, denominators, and the native scientific-output gate."""

import json
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

import profile_valgrind as profile


class CachegrindParserTests(unittest.TestCase):
    def parse(self, text):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "events.out"
            path.write_text(text)
            return profile.parse_cachegrind(path)

    def test_event_names_control_denominators_even_when_columns_are_reordered(self):
        events = {
            "Ir": 1000,
            "I1mr": 20,
            "ILmr": 4,
            "Dr": 300,
            "D1mr": 30,
            "DLmr": 5,
            "Dw": 100,
            "D1mw": 10,
            "DLmw": 3,
            "Bc": 200,
            "Bcm": 20,
            "Bi": 50,
            "Bim": 10,
        }
        order = list(reversed(events))
        result = self.parse(
            "desc: I1 cache: 32768 B, 64 B, 8-way associative\n"
            "events: "
            + " ".join(order)
            + "\nsummary: "
            + " ".join(str(events[name]) for name in order)
            + "\n"
        )
        self.assertEqual(result["counts"], events)
        self.assertEqual(result["derived_counts"]["LL_requests_after_L1_misses"], 60)
        rates = result["simulated_rates"]
        expected = {
            "I1_misses_per_instruction": (20, 1000),
            "D1_misses_per_data_reference": (40, 400),
            "LL_misses_per_LL_request_local": (12, 60),
            "LL_misses_per_all_reference_global": (12, 1400),
            "LL_instruction_misses_per_instruction": (4, 1000),
            "LL_data_misses_per_data_reference": (8, 400),
            "branch_mispredictions_per_branch": (30, 250),
            "conditional_mispredictions_per_conditional_branch": (20, 200),
            "indirect_mispredictions_per_indirect_branch": (10, 50),
        }
        for name, (numerator, denominator) in expected.items():
            with self.subTest(rate=name):
                self.assertEqual(rates[name]["numerator"], numerator)
                self.assertEqual(rates[name]["denominator"], denominator)
                self.assertAlmostEqual(rates[name]["fraction"], numerator / denominator)

    def test_zero_activity_has_no_defined_miss_rate(self):
        result = self.parse(
            "events: "
            + " ".join(profile.CACHE_EVENTS)
            + "\nsummary: "
            + " ".join("0" for _ in profile.CACHE_EVENTS)
            + "\n"
        )
        for rate in result["simulated_rates"].values():
            self.assertEqual(rate["denominator"], 0)
            self.assertIsNone(rate["fraction"])

    def test_incomplete_or_ambiguous_profile_is_rejected(self):
        header = "events: " + " ".join(profile.CACHE_EVENTS) + "\n"
        totals = "summary: " + " ".join("1" for _ in profile.CACHE_EVENTS) + "\n"
        for text in (
            header,
            totals,
            header + "summary: 1 2\n",
            header + header + totals,
            header + totals + totals,
            header + totals.replace("summary: 1", "summary: -1", 1),
        ):
            with self.subTest(text=text), self.assertRaises(RuntimeError):
                self.parse(text)


class DhatParserTests(unittest.TestCase):
    def test_global_peak_uses_simultaneous_live_bytes_not_site_local_maxima(self):
        fixture = {
            "dhatFileVersion": 2,
            "mode": "heap",
            "tu": "instrs",
            "te": 10000,
            "tg": 2000,
            "ftbl": ["[root]", "allocate_first", "allocate_second"],
            "pps": [
                {
                    "tb": 120,
                    "tbk": 2,
                    "tl": 500,
                    "mb": 100,
                    "mbk": 2,
                    "gb": 30,
                    "gbk": 1,
                    "eb": 0,
                    "ebk": 0,
                    "rb": 50,
                    "wb": 100,
                    "fs": [1],
                },
                {
                    "tb": 220,
                    "tbk": 3,
                    "tl": 1000,
                    "mb": 200,
                    "mbk": 3,
                    "gb": 40,
                    "gbk": 2,
                    "eb": 10,
                    "ebk": 1,
                    "rb": 20,
                    "wb": 200,
                    "fs": [2],
                },
            ],
        }
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "dhat.json"
            path.write_text(json.dumps(fixture))
            result = profile.parse_dhat(path)
        self.assertEqual(result["total_requested_bytes"], 340)
        self.assertEqual(result["total_allocation_blocks_including_realloc"], 5)
        self.assertEqual(result["peak_live_heap_bytes"], 70)
        self.assertEqual(result["blocks_at_global_heap_peak"], 3)
        self.assertEqual(result["end_live_heap_bytes"], 10)
        self.assertEqual(result["average_lifetime_instructions"], 300)
        self.assertEqual(result["global_heap_peak_at_instruction"], 2000)
        self.assertEqual(
            result["top_allocation_sites_by_requested_bytes"][0]["allocation_stack"],
            ["allocate_second"],
        )

    def test_other_modes_and_non_instruction_units_are_rejected(self):
        for mode, units in (("copy", "instrs"), ("heap", "seconds")):
            with self.subTest(mode=mode, units=units):
                with tempfile.TemporaryDirectory() as temporary:
                    path = Path(temporary) / "dhat.json"
                    path.write_text(
                        json.dumps(
                            {
                                "dhatFileVersion": 2,
                                "mode": mode,
                                "tu": units,
                                "pps": [],
                            }
                        )
                    )
                    with self.assertRaises(RuntimeError):
                        profile.parse_dhat(path)


class ProfilerOutputOracleTests(unittest.TestCase):
    def test_small_scientific_change_rejects_profile_before_acceptance(self):
        workload = {"kind": "descriptors", "expected_records": 1}
        with tempfile.TemporaryDirectory() as temporary:
            folder = Path(temporary)
            stdout = folder / "stdout.txt"
            stdout.write_text('[{"conformers":1,"sterimol_l":1.0000001}]\n')
            (folder / "stderr.txt").write_text("")
            native = profile.check_result(workload, folder)
            self.assertEqual(
                profile.check_against_native(workload, folder, native), native
            )
            stdout.write_text('[{"conformers":1,"sterimol_l":1.0000002}]\n')
            with self.assertRaisesRegex(RuntimeError, "differs from native oracle"):
                profile.check_against_native(workload, folder, native)

    def test_matching_stdout_cannot_hide_a_skipped_input(self):
        workload = {"kind": "descriptors", "expected_records": 1}
        with tempfile.TemporaryDirectory() as temporary:
            folder = Path(temporary)
            (folder / "stdout.txt").write_text('[{"conformers":1}]\n')
            stderr = folder / "stderr.txt"
            stderr.write_text("")
            native = profile.check_result(workload, folder)
            stderr.write_text("skipped ligand.xyz: no donor\n")
            with self.assertRaisesRegex(RuntimeError, "skipped inputs"):
                profile.check_against_native(workload, folder, native)


if __name__ == "__main__":
    unittest.main()
