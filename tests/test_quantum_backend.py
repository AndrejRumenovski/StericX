"""Unit tests for deterministic quantum-backend parsing and center selection."""

from __future__ import annotations

import tempfile
import time
import unittest
from pathlib import Path
from unittest.mock import patch

import numpy as np

from stericx_quantum import (
    ANGSTROM_TO_BOHR,
    CacheLock,
    LmoCenterResult,
    QuantumBackend,
    QuantumBackendError,
    QuantumConfig,
    donor_neighbor_indices,
    parse_crest_summary,
    parse_lmocent_coord,
    parse_xyz_frames,
    select_kraken_lmo_center,
)


class QuantumBackendParsingTests(unittest.TestCase):
    """Exercise file formats without requiring CREST or xTB executables."""

    def test_parses_multiframe_crest_xyz_energies(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "crest_conformers.xyz"
            path.write_text(
                "2\n-10.000000\nP 0 0 0\nC 1 0 0\n2\n-9.998000\nP 0 0 0\nC 0 1 0\n",
                encoding="utf-8",
            )
            frames = parse_xyz_frames(path)
        self.assertEqual(len(frames), 2)
        self.assertEqual(frames[0].elements, ("P", "C"))
        self.assertEqual(frames[0].energy_hartree, -10.0)
        self.assertEqual(frames[1].energy_hartree, -9.998)

    def test_parses_xtb_lmo_centers_from_bohr(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "lmocent.coord"
            path.write_text(
                "$coord\n"
                f" 0.0 0.0 {-ANGSTROM_TO_BOHR:.12f} he\n"
                f" {ANGSTROM_TO_BOHR:.12f} 0.0 0.0 he\n"
                " 0.0 0.0 0.0 p\n"
                "$end\n",
                encoding="utf-8",
            )
            centers = parse_lmocent_coord(path)
        np.testing.assert_allclose(
            centers,
            [[0.0, 0.0, -1.0], [1.0, 0.0, 0.0]],
            atol=1.0e-10,
        )

    def test_selects_lone_pair_lmo_farthest_from_substituents(self) -> None:
        elements = ("P", "C", "C", "C")
        coordinates = np.asarray(
            [
                [0.0, 0.0, 0.0],
                [1.4, 0.0, 0.45],
                [-0.7, 1.212, 0.45],
                [-0.7, -1.212, 0.45],
            ]
        )
        lmo_centers = np.asarray(
            [
                [0.0, 0.0, -1.0],
                [0.7, 0.0, 0.2],
                [-0.35, 0.6, 0.2],
                [-0.35, -0.6, 0.2],
                [5.0, 5.0, 5.0],
            ]
        )
        center, selected, neighbors = select_kraken_lmo_center(
            elements,
            coordinates,
            0,
            lmo_centers,
        )
        self.assertEqual(selected, 0)
        self.assertEqual(set(neighbors), {1, 2, 3})
        np.testing.assert_allclose(center, [0.0, 0.0, -2.1], atol=1.0e-12)

    def test_neighbor_selection_is_distance_then_index_stable(self) -> None:
        elements = ("P", "C", "C", "H", "C", "C")
        coordinates = np.asarray(
            [
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [-1.0, 0.0, 0.0],
                [0.0, 0.1, 0.0],
                [0.0, 1.0, 0.0],
                [0.0, 2.0, 0.0],
            ]
        )
        self.assertEqual(
            donor_neighbor_indices(elements, coordinates, 0),
            (1, 2, 4),
        )

    def test_parses_crest_212_population_table(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "crest.log"
            path.write_text(
                "CREST Version 2.12\n"
                "Erel/kcal        Etot weight/tot  conformer     set   degen"
                "     origin\n"
                " 1 0.000 -10.000 1 0.750 1 3 MTD\n"
                " 2 1.000 -9.998 2 0.250 2 1 GC\n"
                "T /K 298.15\n"
                "CREST terminated normally.\n",
                encoding="utf-8",
            )
            rows = parse_crest_summary(path)
        self.assertEqual(len(rows), 2)
        self.assertEqual(rows[0]["relative_energy_kcal_mol"], 0.0)
        self.assertEqual(rows[0]["boltzmann_weight"], 0.75)
        self.assertEqual(rows[0]["degeneracy"], 3)


class CacheRecoveryTests(unittest.TestCase):
    """Validate lock ownership, stale recovery, and duplicate-safe batching."""

    def test_live_lock_cannot_be_stolen(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "calculation.lock"
            lock = CacheLock.acquire(path, stale_after_seconds=60.0)
            try:
                with self.assertRaisesRegex(
                    QuantumBackendError,
                    "already locked",
                ):
                    CacheLock.acquire(path, stale_after_seconds=60.0)
            finally:
                lock.release()
            self.assertFalse(path.exists())

    def test_dead_local_owner_is_reclaimed(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "calculation.lock"
            path.write_text(
                '{"hostname": "'
                + __import__("socket").gethostname()
                + '", "pid": 999999999, "process_start_token": "missing"}\n',
                encoding="utf-8",
            )
            lock = CacheLock.acquire(path, stale_after_seconds=60.0)
            lock.release()
            self.assertFalse(path.exists())

    def test_owner_does_not_remove_replacement_lock(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "calculation.lock"
            lock = CacheLock.acquire(path, stale_after_seconds=60.0)
            path.unlink()
            path.write_text(
                '{"token": "replacement", "hostname": "remote"}\n',
                encoding="utf-8",
            )
            lock.release()
            self.assertTrue(path.is_file())
            self.assertIn("replacement", path.read_text(encoding="utf-8"))

    def test_parallel_batch_deduplicates_identical_geometries(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            first = root / "first.xyz"
            duplicate = root / "duplicate.xyz"
            distinct = root / "distinct.xyz"
            first.write_text("1\nfirst\nP 0 0 0\n", encoding="utf-8")
            duplicate.write_bytes(first.read_bytes())
            distinct.write_text("1\ndistinct\nP 1 0 0\n", encoding="utf-8")
            backend = QuantumBackend.__new__(QuantumBackend)
            backend.config = QuantumConfig(
                cache_dir=root / "cache",
                threads=1,
                lmo_workers=8,
            )
            calls: list[Path] = []

            def fake_lmo(path: Path, _: int) -> LmoCenterResult:
                calls.append(path)
                time.sleep(0.01)
                return LmoCenterResult(
                    center_angstrom=(0.0, 0.0, -2.1),
                    selected_lmo_angstrom=(0.0, 0.0, -1.0),
                    selected_lmo_index=0,
                    donor_neighbors=(1, 2, 3),
                    cache_key=path.read_text(encoding="utf-8"),
                    cache_hit=False,
                    manifest_path=str(path) + ".json",
                )

            with (
                patch.object(backend, "lmo_center", side_effect=fake_lmo),
                patch("stericx_quantum.os.cpu_count", return_value=2),
            ):
                batch = backend.lmo_centers(
                    [first, duplicate, distinct],
                    donor_idx=0,
                )
            self.assertEqual(len(calls), 2)
            self.assertEqual(batch.effective_workers, 2)
            self.assertEqual(len(batch.results), 3)
            self.assertEqual(
                batch.results[0].cache_key,
                batch.results[1].cache_key,
            )


class SplitQuantumCacheIntegrationTests(unittest.TestCase):
    """Exercise split-stage caching and interrupted-job recovery end to end."""

    @staticmethod
    def _write_executables(root: Path) -> tuple[Path, Path]:
        xtb = root / "fake_xtb.py"
        crest = root / "fake_crest.py"
        xtb.write_text(
            """#!/usr/bin/env python3
import sys
from pathlib import Path
if "--version" in sys.argv:
    print("xtb version 6.4.0")
    raise SystemExit(0)
Path("lmocent.coord").write_text(
    "$coord\\n"
    " 0.000000000 0.000000000 -1.889725989 he\\n"
    " 1.322808192 0.000000000 0.377945198 he\\n"
    "-0.661404096 1.133835593 0.377945198 he\\n"
    "-0.661404096 -1.133835593 0.377945198 he\\n"
    "$end\\n"
)
print("normal xTB property completion")
""",
            encoding="utf-8",
        )
        crest.write_text(
            """#!/usr/bin/env python3
import sys
from pathlib import Path
if "--version" in sys.argv:
    print("Version 2.12")
    raise SystemExit(0)
Path("crest_conformers.xyz").write_text(
    "4\\n-10.000000\\n"
    "P 0 0 0\\nC 1.4 0 0.45\\nC -0.7 1.212 0.45\\nC -0.7 -1.212 0.45\\n"
    "4\\n-9.999000\\n"
    "P 0 0 0\\nC 1.4 0 0.55\\nC -0.7 1.212 0.55\\nC -0.7 -1.212 0.55\\n"
)
print("Erel/kcal Etot weight/tot conformer set degen origin")
print("1 0.000 -10.000 1 0.750 1 1 MTD")
print("2 0.628 -9.999 2 0.250 2 1 MTD")
print("T /K : 298.15")
print("CREST terminated normally.")
""",
            encoding="utf-8",
        )
        xtb.chmod(0o755)
        crest.chmod(0o755)
        return xtb, crest

    @staticmethod
    def _write_input(root: Path) -> Path:
        path = root / "input.xyz"
        path.write_text(
            "4\ninput\nP 0 0 0\nC 1.4 0 0.45\nC -0.7 1.212 0.45\nC -0.7 -1.212 0.45\n",
            encoding="utf-8",
        )
        return path

    def _backend(self, root: Path) -> QuantumBackend:
        xtb, crest = self._write_executables(root)
        return QuantumBackend(
            QuantumConfig(
                cache_dir=root / "cache",
                xtb_executable=xtb,
                crest_executable=crest,
                threads=1,
                lmo_workers=2,
            )
        )

    def test_split_stages_cache_and_replay(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            backend = self._backend(root)
            input_xyz = self._write_input(root)
            first = backend.conformer_ensemble(input_xyz, donor_idx=0)
            second = backend.conformer_ensemble(input_xyz, donor_idx=0)
            self.assertFalse(first.cache_hit)
            self.assertFalse(first.crest_cache_hit)
            self.assertEqual(first.lmo_cache_misses, 2)
            self.assertEqual(first.effective_lmo_workers, 2)
            self.assertTrue(second.cache_hit)
            self.assertTrue(second.crest_cache_hit)
            self.assertEqual(second.lmo_cache_hits, 2)
            self.assertEqual(second.lmo_cache_misses, 0)
            self.assertEqual(len(second.conformers), 2)
            self.assertTrue((root / "cache" / "crest" / first.crest_cache_key).is_dir())
            self.assertTrue((root / "cache" / "jobs" / first.cache_key).is_dir())

    def test_failed_lmo_job_resumes_from_independent_caches(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            backend = self._backend(root)
            input_xyz = self._write_input(root)
            original = backend.lmo_center
            calls = 0

            def fail_second(path: Path, donor_idx: int) -> LmoCenterResult:
                nonlocal calls
                calls += 1
                if calls == 2:
                    raise QuantumBackendError("injected LMO interruption")
                return original(path, donor_idx)

            backend.config = QuantumConfig(
                **{
                    **backend.config.__dict__,
                    "lmo_workers": 1,
                }
            )
            with patch.object(backend, "lmo_center", side_effect=fail_second):
                with self.assertRaisesRegex(
                    QuantumBackendError,
                    "injected LMO interruption",
                ):
                    backend.conformer_ensemble(input_xyz, donor_idx=0)

            states = list((root / "cache" / "jobs").glob("*/state.json"))
            self.assertEqual(len(states), 1)
            self.assertIn('"status": "failed"', states[0].read_text())
            resumed = backend.conformer_ensemble(input_xyz, donor_idx=0)
            self.assertTrue(resumed.job_resumed)
            self.assertEqual(resumed.lmo_cache_hits, 1)
            self.assertEqual(resumed.lmo_cache_misses, 1)


if __name__ == "__main__":
    unittest.main()
