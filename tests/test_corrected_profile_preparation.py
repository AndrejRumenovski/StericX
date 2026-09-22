"""Admission and preservation contracts for corrected profiling inputs."""

from __future__ import annotations

import argparse
import copy
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import prepare_corrected_profile as prep  # noqa: E402


def conformer(molecule, conformer_id, *, degree=3, donors=1, status="submitted"):
    return {
        "molecule_id": molecule,
        "conformer_id": conformer_id,
        "status": status,
        "phosphorus_atoms": [{"degree": degree}] * donors,
    }


class CorrectedProfilePreparationTests(unittest.TestCase):
    def test_metadata_policy_allows_relocation_but_not_scientific_token_changes(self):
        original = {
            "Source_ID": "2064",
            "Temp_K": "298.15",
            "Ligand_XYZ_Path": "xyz/ligand.xyz",
            "Conformer_XYZ_Paths": "conformers/first.xyz;conformers/second.xyz",
            "Conformer_Boltzmann_Weights": "0.1234500000;0.8765500000",
            "Exp_ddG_kcal_mol": "0.056166692",
            "Dataset_Split": "train",
        }
        corrected = {
            **original,
            "Temp_K": "353.15",
            "Ligand_XYZ_Path": "/relocated/ligand.xyz",
            "Conformer_XYZ_Paths": "/relocated/first.xyz;/relocated/second.xyz",
            "Historical_Recorded_Temp_K": "298.15",
            "Conformer_Population_Source": (
                "Supplied historical MMFF94 populations at 298.15 K; "
                "retained without reweighting"
            ),
            "Target_Source_Note": (
                "Published ddG_abs retained; source ee/DDG inconsistency "
                "remains unresolved"
            ),
        }
        fields, new_fields = list(original), list(corrected)
        prep.metadata_policy(fields, [original], new_fields, [corrected])
        for key, replacement in (
            ("Conformer_Boltzmann_Weights", "0.1;0.9"),
            ("Exp_ddG_kcal_mol", "0.05616669"),
            ("Dataset_Split", "blind"),
            ("Historical_Recorded_Temp_K", "353.15"),
            ("Target_Source_Note", "Published ddG_abs retained"),
        ):
            changed = copy.deepcopy(corrected)
            changed[key] = replacement
            with self.subTest(key=key), self.assertRaises(ValueError):
                prep.metadata_policy(fields, [original], new_fields, [changed])

    def test_bound_copy_refuses_stale_source_and_corrupted_copy(self):
        with tempfile.TemporaryDirectory() as directory:
            source, target = Path(directory) / "source", Path(directory) / "copy"
            source.write_bytes(b"original")
            binding = prep.record(source)
            copied = prep.copy_bound(binding, target)
            self.assertTrue(prep.same_content(binding, copied))
            with self.assertRaises(FileExistsError):
                prep.copy_bound(binding, target)
            source.write_bytes(b"modified")
            with self.assertRaises(ValueError):
                prep.copy_bound(binding, Path(directory) / "stale")
            binding = prep.record(source)

            def corrupt(_source, destination):
                destination.write_bytes(b"tampered")

            with (
                patch.object(prep.shutil, "copymode", side_effect=corrupt),
                self.assertRaisesRegex(ValueError, "copied input"),
            ):
                prep.copy_bound(binding, Path(directory) / "corrupted")

    def test_failed_full_build_proof_aborts_before_any_native_execution(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            with (
                patch.object(
                    prep.replay, "verify_build", side_effect=ValueError("changed proof")
                ) as verifier,
                patch.object(prep.subprocess, "run") as run,
                self.assertRaisesRegex(ValueError, "changed proof"),
            ):
                prep.prepare(argparse.Namespace(root=root, build=root / "build"))
            verifier.assert_called_once_with(root / "build/manifest.json")
            run.assert_not_called()
            self.assertFalse((root / "workloads").exists())

    def test_admission_is_whole_ensemble_and_independent_of_input_order(self):
        rows = [
            conformer(1, 2),
            conformer(2, 0),
            conformer(3, 0, degree=4),
            conformer(1, 1),
            conformer(2, 1, donors=2),
            conformer(4, 0, status="no_DFT_geometry_available"),
        ]
        selected, excluded = prep.eligible_ensembles(rows)
        self.assertEqual(list(selected), [1])
        self.assertEqual([r["conformer_id"] for r in selected[1]], [1, 2])
        self.assertEqual([r["molecule_id"] for r in excluded], [2, 3, 4])
        self.assertEqual(len(excluded[0]["inventory_rows"]), 2)
        self.assertEqual(prep.eligible_ensembles(rows[::-1])[0], selected)

    def test_csv_preserves_decimal_tokens_weights_and_declared_split(self):
        fields = ["coordinates", "weights", "split", "temperature"]
        rows = [
            {
                "coordinates": "1.0000000001;-0.0000000000001",
                "weights": "0.1234500000;0.8765500000",
                "split": "blind",
                "temperature": "353.15",
            }
        ]
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "rows.csv"
            prep.write_rows(path, fields, rows)
            self.assertEqual(prep.read_rows(path), (fields, rows))
            with self.assertRaises(FileExistsError):
                prep.write_rows(path, fields, rows)

    def test_bound_symlink_payload_mutation_is_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            source = Path(directory) / "source.body"
            source.write_bytes(b"original")
            access = Path(directory) / "conformer.sdf"
            access.symlink_to(source)
            item = prep.record(access)
            prep.verify(item)
            manifest = {"files": [{**item, "path": str(access)}]}
            prep.profile.verify_manifest(manifest)
            source.write_bytes(b"modified")
            with self.assertRaises(ValueError):
                prep.verify(item)
            with self.assertRaises(RuntimeError):
                prep.profile.verify_manifest(manifest)

    def test_existing_preparation_is_never_overwritten_or_executed(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            work = root / "workloads"
            work.mkdir()
            sentinel = work / "failed-evidence"
            sentinel.write_text("preserved")
            with (
                patch.object(prep.subprocess, "run") as run,
                self.assertRaises(FileExistsError),
            ):
                prep.prepare(argparse.Namespace(root=root))
            run.assert_not_called()
            self.assertEqual(sentinel.read_text(), "preserved")

    def test_original_audit_destination_rejected_before_execution(self):
        with (
            patch.object(prep.subprocess, "run") as run,
            self.assertRaisesRegex(ValueError, "immutable"),
        ):
            prep.prepare(argparse.Namespace(root=prep.AUDIT / "new_work"))
        run.assert_not_called()


if __name__ == "__main__":
    unittest.main()
