"""Production CREST/xTB execution and Kraken lone-pair-center extraction.

This module deliberately contains no chemistry-specific model fitting. It owns
the expensive, reproducible boundary between Cartesian input structures and
quantum-derived conformer ensembles:

* explicit executable discovery and version/hash provenance;
* shell-free subprocess execution with timeouts and bounded thread settings;
* content-addressed, immutable result caches;
* CREST multi-XYZ and Boltzmann population parsing;
* xTB ``lmocent.coord`` parsing and the exact public Kraken center-selection
  rule.
"""

from __future__ import annotations

import hashlib
import json
import math
import os
import re
import shutil
import socket
import subprocess
import threading
import time
import uuid
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import asdict, dataclass
from datetime import UTC, datetime
from pathlib import Path
from typing import Final

import numpy as np

ANGSTROM_TO_BOHR: Final[float] = 1.889725989
HARTREE_TO_KCAL_MOL: Final[float] = 627.5094740631
GAS_CONSTANT_KCAL_MOL_K: Final[float] = 0.00198720425864083
CACHE_SCHEMA_VERSION: Final[int] = 1
NORMAL_TERMINATION: Final[str] = "CREST terminated normally."


class QuantumBackendError(RuntimeError):
    """Raised when a tool, calculation, or cached artifact is invalid."""


@dataclass(frozen=True)
class QuantumConfig:
    """Frozen settings for the Kraken-compatible quantum backend."""

    cache_dir: Path
    xtb_executable: Path | None = None
    crest_executable: Path | None = None
    threads: int = 1
    charge: int = 0
    uhf: int = 0
    temperature_k: float = 298.15
    energy_window_kcal_mol: float = 6.0
    center_distance_angstrom: float = 2.1
    crest_timeout_seconds: float = 21_600.0
    xtb_timeout_seconds: float = 1_800.0
    solvent: str = "toluene"
    quick: bool = False
    lmo_workers: int = 4
    stale_lock_seconds: float = 86_400.0

    def validate(self) -> None:
        """Reject settings that cannot define a physical deterministic run."""
        if self.threads <= 0:
            raise ValueError("threads must be positive")
        if self.uhf < 0:
            raise ValueError("uhf must be non-negative")
        if not math.isfinite(self.temperature_k) or self.temperature_k <= 0.0:
            raise ValueError("temperature must be positive and finite")
        if (
            not math.isfinite(self.energy_window_kcal_mol)
            or self.energy_window_kcal_mol <= 0.0
        ):
            raise ValueError("energy window must be positive and finite")
        if (
            not math.isfinite(self.center_distance_angstrom)
            or self.center_distance_angstrom <= 0.0
        ):
            raise ValueError("center distance must be positive and finite")
        if self.crest_timeout_seconds <= 0.0 or self.xtb_timeout_seconds <= 0.0:
            raise ValueError("tool timeouts must be positive")
        if self.lmo_workers <= 0:
            raise ValueError("LMO workers must be positive")
        if self.stale_lock_seconds <= 0.0:
            raise ValueError("stale-lock timeout must be positive")


@dataclass(frozen=True)
class ToolInfo:
    """Executable identity included in every content-addressed cache key."""

    name: str
    path: str
    version: str
    sha256: str


@dataclass(frozen=True)
class XyzFrame:
    """One XYZ structure and optional comment-line electronic energy."""

    elements: tuple[str, ...]
    coordinates: np.ndarray
    comment: str
    energy_hartree: float | None


@dataclass(frozen=True)
class LmoCenterResult:
    """Selected xTB LMO direction and normalized virtual metal center."""

    center_angstrom: tuple[float, float, float]
    selected_lmo_angstrom: tuple[float, float, float]
    selected_lmo_index: int
    donor_neighbors: tuple[int, int, int]
    cache_key: str
    cache_hit: bool
    manifest_path: str


@dataclass(frozen=True)
class LmoBatchResult:
    """Ordered results and execution statistics for a parallel LMO batch."""

    results: tuple[LmoCenterResult, ...]
    cache_hits: int
    cache_misses: int
    effective_workers: int


@dataclass(frozen=True)
class QuantumConformer:
    """One retained CREST conformer with an xTB-derived coordination center."""

    index: int
    xyz_path: str
    energy_hartree: float
    relative_energy_kcal_mol: float
    boltzmann_weight: float
    degeneracy: int
    coordination_center_angstrom: tuple[float, float, float]
    lmo_cache_key: str


@dataclass(frozen=True)
class CrestConformer:
    """One durable conformer produced by the CREST-only cache stage."""

    index: int
    xyz_path: str
    xyz_sha256: str
    energy_hartree: float
    relative_energy_kcal_mol: float
    boltzmann_weight: float
    degeneracy: int


@dataclass(frozen=True)
class CrestEnsemble:
    """A CREST ensemble committed before downstream xTB property jobs."""

    cache_key: str
    cache_hit: bool
    manifest_path: str
    conformers: tuple[CrestConformer, ...]


@dataclass(frozen=True)
class QuantumEnsemble:
    """Immutable content-addressed CREST/xTB result."""

    cache_key: str
    cache_hit: bool
    manifest_path: str
    conformers: tuple[QuantumConformer, ...]
    xtb: ToolInfo
    crest: ToolInfo
    crest_cache_key: str
    crest_cache_hit: bool
    lmo_cache_hits: int
    lmo_cache_misses: int
    job_resumed: bool
    effective_lmo_workers: int


def sha256_file(path: Path) -> str:
    """Hash a file without loading an executable or trajectory into memory."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def atomic_write_json(path: Path, value: object) -> None:
    """Atomically replace one JSON manifest."""
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(
        json.dumps(value, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    os.replace(temporary, path)


def parse_xyz_frames(path: Path) -> list[XyzFrame]:
    """Parse one or more concatenated XYZ frames with stable atom ordering."""
    lines = path.read_text(encoding="utf-8").splitlines()
    frames: list[XyzFrame] = []
    offset = 0
    while offset < len(lines):
        while offset < len(lines) and not lines[offset].strip():
            offset += 1
        if offset >= len(lines):
            break
        try:
            atom_count = int(lines[offset].strip())
        except ValueError as exc:
            raise QuantumBackendError(
                f"{path} has invalid XYZ atom count on line {offset + 1}"
            ) from exc
        if atom_count <= 0 or offset + atom_count + 2 > len(lines):
            raise QuantumBackendError(f"{path} contains a truncated XYZ frame")
        comment = lines[offset + 1].strip()
        energy = _comment_energy(comment)
        elements: list[str] = []
        coordinates: list[list[float]] = []
        for line_number, line in enumerate(
            lines[offset + 2 : offset + atom_count + 2],
            start=offset + 3,
        ):
            fields = line.split()
            if len(fields) < 4:
                raise QuantumBackendError(
                    f"{path} atom line {line_number} has fewer than four fields"
                )
            try:
                coordinate = [float(value) for value in fields[1:4]]
            except ValueError as exc:
                raise QuantumBackendError(
                    f"{path} atom line {line_number} has invalid coordinates"
                ) from exc
            if not np.isfinite(coordinate).all():
                raise QuantumBackendError(
                    f"{path} atom line {line_number} is non-finite"
                )
            elements.append(fields[0])
            coordinates.append(coordinate)
        frames.append(
            XyzFrame(
                elements=tuple(elements),
                coordinates=np.asarray(coordinates, dtype=float),
                comment=comment,
                energy_hartree=energy,
            )
        )
        offset += atom_count + 2
    if not frames:
        raise QuantumBackendError(f"{path} contains no XYZ frames")
    reference_elements = frames[0].elements
    if any(frame.elements != reference_elements for frame in frames[1:]):
        raise QuantumBackendError(
            f"{path} changes atom identity or order between frames"
        )
    return frames


def write_xyz_frame(path: Path, frame: XyzFrame, comment: str | None = None) -> None:
    """Write one XYZ frame while preserving atom order."""
    lines = [
        str(len(frame.elements)),
        comment if comment is not None else frame.comment,
    ]
    lines.extend(
        f"{element:<3s} {coordinate[0]: .12f} "
        f"{coordinate[1]: .12f} {coordinate[2]: .12f}"
        for element, coordinate in zip(
            frame.elements,
            frame.coordinates,
            strict=True,
        )
    )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def parse_lmocent_coord(path: Path) -> np.ndarray:
    """Extract helium LMO centers from xTB's Turbomole-format output.

    xTB writes coordinates in bohr. Kraken converts these to ångströms by
    dividing by 1.889725989 before choosing the phosphorus lone-pair center.
    """
    if not path.is_file():
        raise QuantumBackendError(f"xTB LMO center file is missing: {path}")
    centers: list[list[float]] = []
    for line_number, line in enumerate(
        path.read_text(encoding="utf-8").splitlines(),
        start=1,
    ):
        fields = line.split()
        if len(fields) != 4 or fields[0].startswith("$"):
            continue
        try:
            coordinate = [
                float(fields[0]) / ANGSTROM_TO_BOHR,
                float(fields[1]) / ANGSTROM_TO_BOHR,
                float(fields[2]) / ANGSTROM_TO_BOHR,
            ]
        except ValueError as exc:
            raise QuantumBackendError(
                f"{path} line {line_number} has invalid LMO coordinates"
            ) from exc
        if fields[3].casefold() == "he":
            centers.append(coordinate)
    if not centers:
        raise QuantumBackendError(f"{path} contains no helium LMO centers")
    result = np.asarray(centers, dtype=float)
    if not np.isfinite(result).all():
        raise QuantumBackendError(f"{path} contains non-finite LMO centers")
    return result


def donor_neighbor_indices(
    elements: tuple[str, ...] | list[str],
    coordinates: np.ndarray,
    donor_idx: int,
) -> tuple[int, int, int]:
    """Select the three nearest heavy atoms, matching Kraken's XYZ topology."""
    if not 0 <= donor_idx < len(elements):
        raise QuantumBackendError(f"donor index {donor_idx} is out of bounds")
    candidates = sorted(
        (
            (
                float(np.sum((coordinates[index] - coordinates[donor_idx]) ** 2)),
                index,
            )
            for index, element in enumerate(elements)
            if index != donor_idx and element.casefold() != "h"
        ),
        key=lambda item: (item[0], item[1]),
    )
    if len(candidates) < 3:
        raise QuantumBackendError("three heavy donor neighbors are required")
    return tuple(index for _, index in candidates[:3])  # type: ignore[return-value]


def select_kraken_lmo_center(
    elements: tuple[str, ...] | list[str],
    coordinates: np.ndarray,
    donor_idx: int,
    lmo_centers_angstrom: np.ndarray,
    center_distance_angstrom: float = 2.1,
) -> tuple[np.ndarray, int, tuple[int, int, int]]:
    """Apply Kraken's exact free-phosphine LMO-center selection rule.

    The four LMO centers nearest phosphorus are considered. For each candidate,
    its nearest distance to any of the three P substituents is calculated. The
    candidate maximizing that minimum distance defines the lone-pair direction.
    The returned virtual center is normalized to the configured P-metal distance.
    """
    neighbors = donor_neighbor_indices(elements, coordinates, donor_idx)
    donor = coordinates[donor_idx]
    if lmo_centers_angstrom.ndim != 2 or lmo_centers_angstrom.shape[1] != 3:
        raise QuantumBackendError("LMO centers must have shape (n, 3)")
    if len(lmo_centers_angstrom) < 4:
        raise QuantumBackendError("at least four LMO centers are required")
    distances = np.linalg.norm(lmo_centers_angstrom - donor, axis=1)
    nearest_indices = np.argsort(distances, kind="stable")[:4]
    neighbor_coordinates = coordinates[np.asarray(neighbors)]
    minimum_distances = np.asarray(
        [
            np.linalg.norm(
                neighbor_coordinates - lmo_centers_angstrom[index], axis=1
            ).min()
            for index in nearest_indices
        ]
    )
    selected_index = int(nearest_indices[int(np.argmax(minimum_distances))])
    direction = lmo_centers_angstrom[selected_index] - donor
    norm = float(np.linalg.norm(direction))
    if not math.isfinite(norm) or norm <= np.finfo(float).eps:
        raise QuantumBackendError("selected LMO center coincides with donor atom")
    center = donor + center_distance_angstrom * direction / norm
    return center, selected_index, neighbors


def parse_crest_summary(path: Path) -> list[dict[str, float | int | str | None]]:
    """Parse CREST 2.12 conformer energies, populations, and degeneracies."""
    text = path.read_text(encoding="utf-8", errors="replace")
    rows: list[dict[str, float | int | str | None]] = []
    reading = False
    for line in text.splitlines():
        if "T /K" in line and reading:
            break
        if reading:
            fields = re.sub(r"\s+", " ", line).strip().split(" ")
            if len(fields) not in {5, 8}:
                continue
            if len(fields) == 8:
                try:
                    rows.append(
                        {
                            "relative_energy_kcal_mol": float(fields[1]),
                            "energy_hartree": float(fields[2]),
                            "boltzmann_weight": float(fields[4]),
                            "degeneracy": int(fields[6]),
                            "origin": fields[7],
                        }
                    )
                except ValueError:
                    continue
        if "Erel/kcal" in line and "weight/tot" in line:
            reading = True
    return rows


class QuantumBackend:
    """Content-addressed CREST/xTB runner."""

    def __init__(self, config: QuantumConfig):
        config.validate()
        self.config = config
        self.config.cache_dir.mkdir(parents=True, exist_ok=True)
        self.xtb = self._tool_info(
            "xtb",
            config.xtb_executable,
            "STERICX_XTB",
        )
        self.crest = self._tool_info(
            "crest",
            config.crest_executable,
            "STERICX_CREST",
        )

    def lmo_center(self, xyz_path: Path, donor_idx: int) -> LmoCenterResult:
        """Calculate or load one conformer's xTB LMO-derived virtual center."""
        input_sha = sha256_file(xyz_path)
        payload = {
            "schema_version": CACHE_SCHEMA_VERSION,
            "operation": "xtb_lmo_center",
            "input_sha256": input_sha,
            "donor_idx": donor_idx,
            "charge": self.config.charge,
            "uhf": self.config.uhf,
            "threads": self.config.threads,
            "solvent": self.config.solvent,
            "center_distance_angstrom": self.config.center_distance_angstrom,
            "xtb": asdict(self.xtb),
        }
        cache_key = _cache_key(payload)
        result_dir = self.config.cache_dir / "lmo" / cache_key
        manifest_path = result_dir / "manifest.json"
        cached = _load_valid_manifest(manifest_path, cache_key)
        if cached is not None:
            return _lmo_result_from_manifest(cached, manifest_path, cache_hit=True)

        lock_path = self.config.cache_dir / "lmo" / f"{cache_key}.lock"
        lock_path.parent.mkdir(parents=True, exist_ok=True)
        cache_lock = CacheLock.acquire(
            lock_path,
            self.config.stale_lock_seconds,
        )
        work_dir = (
            self.config.cache_dir
            / "lmo"
            / f".{cache_key}.{os.getpid()}.{cache_lock.token}.work"
        )
        started = time.monotonic()
        try:
            cached = _load_valid_manifest(manifest_path, cache_key)
            if cached is not None:
                return _lmo_result_from_manifest(
                    cached,
                    manifest_path,
                    cache_hit=True,
                )
            work_dir.mkdir(parents=False, exist_ok=False)
            input_copy = work_dir / "input.xyz"
            shutil.copy2(xyz_path, input_copy)
            command = [
                self.xtb.path,
                "--gbsa",
                self.config.solvent,
                "--lmo",
                "--vfukui",
                "--esp",
                "-P",
                str(self.config.threads),
                "--chrg",
                str(self.config.charge),
            ]
            if self.config.uhf:
                command.extend(["--uhf", str(self.config.uhf)])
            command.append(input_copy.name)
            log_path = work_dir / "xtb.log"
            self._run(command, work_dir, log_path, self.config.xtb_timeout_seconds)
            lmo_path = work_dir / "lmocent.coord"
            centers = parse_lmocent_coord(lmo_path)
            frame = parse_xyz_frames(input_copy)[0]
            center, selected_index, neighbors = select_kraken_lmo_center(
                frame.elements,
                frame.coordinates,
                donor_idx,
                centers,
                self.config.center_distance_angstrom,
            )
            manifest = {
                **payload,
                "cache_key": cache_key,
                "status": "complete",
                "created_at_utc": datetime.now(UTC).isoformat(),
                "elapsed_seconds": time.monotonic() - started,
                "command": command,
                "selected_lmo_index": selected_index,
                "selected_lmo_angstrom": centers[selected_index].tolist(),
                "coordination_center_angstrom": center.tolist(),
                "donor_neighbors": list(neighbors),
                "artifacts": {
                    "input_xyz_sha256": sha256_file(input_copy),
                    "lmocent_coord_sha256": sha256_file(lmo_path),
                    "xtb_log_sha256": sha256_file(log_path),
                },
            }
            atomic_write_json(work_dir / "manifest.json", manifest)
            result_dir.parent.mkdir(parents=True, exist_ok=True)
            os.rename(work_dir, result_dir)
            return _lmo_result_from_manifest(
                manifest,
                manifest_path,
                cache_hit=False,
            )
        except BaseException:
            if work_dir.is_dir():
                shutil.rmtree(work_dir)
            raise
        finally:
            cache_lock.release()

    def lmo_centers(
        self,
        xyz_paths: list[Path] | tuple[Path, ...],
        donor_idx: int,
    ) -> LmoBatchResult:
        """Calculate an ordered LMO batch with bounded duplicate-safe parallelism."""
        if not xyz_paths:
            raise QuantumBackendError("LMO batch cannot be empty")
        hashes = [sha256_file(path) for path in xyz_paths]
        unique_paths: dict[str, Path] = {}
        for xyz_sha, path in zip(hashes, xyz_paths, strict=True):
            unique_paths.setdefault(xyz_sha, path)
        effective_workers = self._effective_lmo_workers(len(unique_paths))
        results_by_hash: dict[str, LmoCenterResult] = {}
        with ThreadPoolExecutor(
            max_workers=effective_workers,
            thread_name_prefix="stericx-lmo",
        ) as executor:
            futures = {
                executor.submit(self.lmo_center, path, donor_idx): xyz_sha
                for xyz_sha, path in unique_paths.items()
            }
            for future in as_completed(futures):
                results_by_hash[futures[future]] = future.result()
        unique_results = tuple(results_by_hash.values())
        return LmoBatchResult(
            results=tuple(results_by_hash[xyz_sha] for xyz_sha in hashes),
            cache_hits=sum(result.cache_hit for result in unique_results),
            cache_misses=sum(not result.cache_hit for result in unique_results),
            effective_workers=effective_workers,
        )

    def conformer_ensemble(self, input_xyz: Path, donor_idx: int) -> QuantumEnsemble:
        """Join a durable CREST stage with resumable parallel xTB LMO jobs."""
        crest_ensemble = self.crest_ensemble(
            input_xyz,
            legacy_donor_idx=donor_idx,
        )
        payload = {
            "schema_version": CACHE_SCHEMA_VERSION,
            "operation": "quantum_ensemble_job",
            "crest_cache_key": crest_ensemble.cache_key,
            "donor_idx": donor_idx,
            "charge": self.config.charge,
            "uhf": self.config.uhf,
            "threads_per_lmo": self.config.threads,
            "solvent": self.config.solvent,
            "center_distance_angstrom": self.config.center_distance_angstrom,
            "xtb": asdict(self.xtb),
        }
        cache_key = _cache_key(payload)
        job_dir = self.config.cache_dir / "jobs" / cache_key
        manifest_path = job_dir / "manifest.json"
        cached = _load_valid_manifest(manifest_path, cache_key)
        if cached is not None:
            return _ensemble_from_job_manifest(
                cached,
                manifest_path,
                self.xtb,
                self.crest,
                cache_hit=True,
            )

        lock_path = self.config.cache_dir / "jobs" / f"{cache_key}.lock"
        lock_path.parent.mkdir(parents=True, exist_ok=True)
        cache_lock = CacheLock.acquire(
            lock_path,
            self.config.stale_lock_seconds,
        )
        job_dir.mkdir(parents=True, exist_ok=True)
        state_path = job_dir / "state.json"
        prior_state = _load_json_object(state_path)
        job_resumed = bool(
            prior_state
            and prior_state.get("status") in {"failed", "interrupted", "lmo_running"}
        )
        started = time.monotonic()
        state: dict[str, object] = {
            **payload,
            "cache_key": cache_key,
            "status": "lmo_running",
            "attempt_id": cache_lock.token,
            "started_at_utc": datetime.now(UTC).isoformat(),
            "crest_manifest_path": crest_ensemble.manifest_path,
            "crest_cache_hit": crest_ensemble.cache_hit,
            "completed_lmo_indices": [],
        }
        try:
            cached = _load_valid_manifest(manifest_path, cache_key)
            if cached is not None:
                return _ensemble_from_job_manifest(
                    cached,
                    manifest_path,
                    self.xtb,
                    self.crest,
                    cache_hit=True,
                )

            unique_paths: dict[str, Path] = {}
            conformer_hashes: list[str] = []
            for conformer in crest_ensemble.conformers:
                conformer_hashes.append(conformer.xyz_sha256)
                unique_paths.setdefault(
                    conformer.xyz_sha256,
                    Path(conformer.xyz_path),
                )
            effective_workers = self._effective_lmo_workers(len(unique_paths))
            state["requested_lmo_workers"] = self.config.lmo_workers
            state["effective_lmo_workers"] = effective_workers
            state["unique_lmo_jobs"] = len(unique_paths)
            atomic_write_json(state_path, state)

            results_by_hash: dict[str, LmoCenterResult] = {}
            completed_indices: set[int] = set()
            with ThreadPoolExecutor(
                max_workers=effective_workers,
                thread_name_prefix="stericx-lmo",
            ) as executor:
                futures = {
                    executor.submit(self.lmo_center, path, donor_idx): xyz_sha
                    for xyz_sha, path in unique_paths.items()
                }
                for future in as_completed(futures):
                    xyz_sha = futures[future]
                    result = future.result()
                    results_by_hash[xyz_sha] = result
                    completed_indices.update(
                        index
                        for index, conformer_sha in enumerate(conformer_hashes)
                        if conformer_sha == xyz_sha
                    )
                    state["completed_lmo_indices"] = sorted(completed_indices)
                    state["last_checkpoint_at_utc"] = datetime.now(UTC).isoformat()
                    atomic_write_json(state_path, state)

            quantum_conformers = tuple(
                QuantumConformer(
                    index=conformer.index,
                    xyz_path=conformer.xyz_path,
                    energy_hartree=conformer.energy_hartree,
                    relative_energy_kcal_mol=conformer.relative_energy_kcal_mol,
                    boltzmann_weight=conformer.boltzmann_weight,
                    degeneracy=conformer.degeneracy,
                    coordination_center_angstrom=(
                        results_by_hash[conformer.xyz_sha256].center_angstrom
                    ),
                    lmo_cache_key=results_by_hash[conformer.xyz_sha256].cache_key,
                )
                for conformer in crest_ensemble.conformers
            )
            lmo_hits = sum(result.cache_hit for result in results_by_hash.values())
            lmo_misses = len(results_by_hash) - lmo_hits
            manifest = {
                **payload,
                "cache_key": cache_key,
                "status": "complete",
                "created_at_utc": datetime.now(UTC).isoformat(),
                "elapsed_seconds": time.monotonic() - started,
                "crest_manifest_path": crest_ensemble.manifest_path,
                "crest_cache_hit": crest_ensemble.cache_hit,
                "lmo_cache_hits": lmo_hits,
                "lmo_cache_misses": lmo_misses,
                "job_resumed": job_resumed,
                "requested_lmo_workers": self.config.lmo_workers,
                "effective_lmo_workers": effective_workers,
                "conformers": [
                    {
                        **asdict(conformer),
                        "coordination_center_angstrom": list(
                            conformer.coordination_center_angstrom
                        ),
                    }
                    for conformer in quantum_conformers
                ],
                "lmo_manifest_paths": sorted(
                    {result.manifest_path for result in results_by_hash.values()}
                ),
            }
            atomic_write_json(manifest_path, manifest)
            state.update(
                {
                    "status": "complete",
                    "completed_lmo_indices": list(
                        range(len(crest_ensemble.conformers))
                    ),
                    "finished_at_utc": datetime.now(UTC).isoformat(),
                    "manifest_path": str(manifest_path),
                }
            )
            atomic_write_json(state_path, state)
            return _ensemble_from_job_manifest(
                manifest,
                manifest_path,
                self.xtb,
                self.crest,
                cache_hit=False,
            )
        except BaseException as exc:
            state.update(
                {
                    "status": (
                        "interrupted"
                        if isinstance(exc, (KeyboardInterrupt, SystemExit))
                        else "failed"
                    ),
                    "error_type": type(exc).__name__,
                    "error": str(exc),
                    "finished_at_utc": datetime.now(UTC).isoformat(),
                }
            )
            atomic_write_json(state_path, state)
            raise
        finally:
            cache_lock.release()

    def _effective_lmo_workers(self, job_count: int) -> int:
        """Bound process fan-out by requested jobs and available CPU threads."""
        if job_count <= 0:
            raise QuantumBackendError("LMO worker pool requires at least one job")
        cpu_bound = max(
            1,
            (os.cpu_count() or 1) // self.config.threads,
        )
        return min(self.config.lmo_workers, job_count, cpu_bound)

    def crest_ensemble(
        self,
        input_xyz: Path,
        legacy_donor_idx: int | None = None,
    ) -> CrestEnsemble:
        """Run or load CREST independently of downstream property calculations."""
        input_sha = sha256_file(input_xyz)
        payload = {
            "schema_version": CACHE_SCHEMA_VERSION,
            "operation": "crest_ensemble",
            "input_sha256": input_sha,
            "charge": self.config.charge,
            "threads": self.config.threads,
            "temperature_k": self.config.temperature_k,
            "energy_window_kcal_mol": self.config.energy_window_kcal_mol,
            "solvent": self.config.solvent,
            "quick": self.config.quick,
            "xtb": asdict(self.xtb),
            "crest": asdict(self.crest),
        }
        cache_key = _cache_key(payload)
        result_dir = self.config.cache_dir / "crest" / cache_key
        manifest_path = result_dir / "manifest.json"
        cached = _load_valid_manifest(manifest_path, cache_key)
        if cached is not None:
            return _crest_ensemble_from_manifest(
                cached,
                result_dir,
                cache_hit=True,
            )

        lock_path = self.config.cache_dir / "crest" / f"{cache_key}.lock"
        lock_path.parent.mkdir(parents=True, exist_ok=True)
        cache_lock = CacheLock.acquire(
            lock_path,
            self.config.stale_lock_seconds,
        )
        work_dir = (
            self.config.cache_dir
            / "crest"
            / f".{cache_key}.{os.getpid()}.{cache_lock.token}.work"
        )
        started = time.monotonic()
        try:
            cached = _load_valid_manifest(manifest_path, cache_key)
            if cached is not None:
                return _crest_ensemble_from_manifest(
                    cached,
                    result_dir,
                    cache_hit=True,
                )
            if legacy_donor_idx is not None:
                migrated = self._migrate_legacy_ensemble(
                    input_xyz,
                    legacy_donor_idx,
                    payload,
                    cache_key,
                    result_dir,
                    work_dir,
                )
                if migrated is not None:
                    return migrated

            work_dir.mkdir(parents=False, exist_ok=False)
            input_copy = work_dir / "input.xyz"
            shutil.copy2(input_xyz, input_copy)
            command = [
                self.crest.path,
                str(input_copy.resolve()),
                "--gbsa",
                self.config.solvent,
                "-metac",
                "-nozs",
                "-T",
                str(self.config.threads),
                "--chrg",
                str(self.config.charge),
                "-xnam",
                self.xtb.path,
            ]
            if self.config.quick:
                command.append("-quick")
            crest_log = work_dir / "crest.log"
            self._run(
                command,
                work_dir,
                crest_log,
                self.config.crest_timeout_seconds,
            )
            log_text = crest_log.read_text(encoding="utf-8", errors="replace")
            if NORMAL_TERMINATION not in log_text:
                raise QuantumBackendError(
                    f"CREST did not report normal termination: {crest_log}"
                )
            trajectory = work_dir / "crest_conformers.xyz"
            frames = parse_xyz_frames(trajectory)
            conformer_data = self._conformer_thermodynamics(frames, crest_log)
            conformer_dir = work_dir / "conformers"
            conformer_dir.mkdir()
            conformers: list[dict[str, object]] = []
            for output_index, (frame, thermo) in enumerate(
                zip(frames, conformer_data, strict=True)
            ):
                if (
                    float(thermo["relative_energy_kcal_mol"])
                    > self.config.energy_window_kcal_mol
                ):
                    continue
                destination = conformer_dir / f"conf_{output_index:03d}.xyz"
                write_xyz_frame(destination, frame)
                energy_hartree = (
                    frame.energy_hartree
                    if frame.energy_hartree is not None
                    else thermo["energy_hartree"]
                )
                conformers.append(
                    {
                        "index": output_index,
                        "xyz_path": str(
                            Path("conformers") / f"conf_{output_index:03d}.xyz"
                        ),
                        "xyz_sha256": sha256_file(destination),
                        "energy_hartree": float(energy_hartree),
                        "relative_energy_kcal_mol": float(
                            thermo["relative_energy_kcal_mol"]
                        ),
                        "boltzmann_weight": float(thermo["boltzmann_weight"]),
                        "degeneracy": int(thermo["degeneracy"]),
                    }
                )
            if not conformers:
                raise QuantumBackendError(
                    "CREST produced no conformers in energy window"
                )
            weight_sum = sum(float(row["boltzmann_weight"]) for row in conformers)
            for row in conformers:
                row["boltzmann_weight"] = float(row["boltzmann_weight"]) / weight_sum
            manifest = {
                **payload,
                "cache_key": cache_key,
                "status": "complete",
                "created_at_utc": datetime.now(UTC).isoformat(),
                "elapsed_seconds": time.monotonic() - started,
                "command": command,
                "conformers": conformers,
                "artifacts": {
                    "input_xyz_sha256": sha256_file(input_copy),
                    "crest_log_sha256": sha256_file(crest_log),
                    "crest_conformers_sha256": sha256_file(trajectory),
                },
            }
            atomic_write_json(work_dir / "manifest.json", manifest)
            os.rename(work_dir, result_dir)
            return _crest_ensemble_from_manifest(
                manifest,
                result_dir,
                cache_hit=False,
            )
        except BaseException:
            if work_dir.is_dir():
                shutil.rmtree(work_dir)
            raise
        finally:
            cache_lock.release()

    def _migrate_legacy_ensemble(
        self,
        input_xyz: Path,
        donor_idx: int,
        crest_payload: dict[str, object],
        crest_cache_key: str,
        result_dir: Path,
        work_dir: Path,
    ) -> CrestEnsemble | None:
        """Promote a valid schema-v1 combined cache into the split CREST stage."""
        legacy_payload = {
            "schema_version": CACHE_SCHEMA_VERSION,
            "operation": "crest_xtb_ensemble",
            "input_sha256": sha256_file(input_xyz),
            "donor_idx": donor_idx,
            "charge": self.config.charge,
            "uhf": self.config.uhf,
            "threads": self.config.threads,
            "temperature_k": self.config.temperature_k,
            "energy_window_kcal_mol": self.config.energy_window_kcal_mol,
            "center_distance_angstrom": self.config.center_distance_angstrom,
            "solvent": self.config.solvent,
            "quick": self.config.quick,
            "xtb": asdict(self.xtb),
            "crest": asdict(self.crest),
        }
        legacy_key = _cache_key(legacy_payload)
        legacy_dir = self.config.cache_dir / "ensembles" / legacy_key
        legacy_manifest_path = legacy_dir / "manifest.json"
        legacy = _load_valid_manifest(legacy_manifest_path, legacy_key)
        if legacy is None or not isinstance(legacy.get("conformers"), list):
            return None

        work_dir.mkdir(parents=False, exist_ok=False)
        conformer_dir = work_dir / "conformers"
        conformer_dir.mkdir()
        conformers: list[dict[str, object]] = []
        for row in legacy["conformers"]:
            if not isinstance(row, dict):
                raise QuantumBackendError(
                    f"invalid legacy ensemble manifest: {legacy_manifest_path}"
                )
            source = legacy_dir / str(row["xyz_path"])
            destination = conformer_dir / Path(str(row["xyz_path"])).name
            shutil.copy2(source, destination)
            conformers.append(
                {
                    "index": int(row["index"]),
                    "xyz_path": str(Path("conformers") / destination.name),
                    "xyz_sha256": sha256_file(destination),
                    "energy_hartree": float(row["energy_hartree"]),
                    "relative_energy_kcal_mol": float(row["relative_energy_kcal_mol"]),
                    "boltzmann_weight": float(row["boltzmann_weight"]),
                    "degeneracy": int(row["degeneracy"]),
                }
            )
        manifest = {
            **crest_payload,
            "cache_key": crest_cache_key,
            "status": "complete",
            "created_at_utc": datetime.now(UTC).isoformat(),
            "elapsed_seconds": 0.0,
            "command": legacy.get("command", []),
            "conformers": conformers,
            "migration": {
                "source": "schema_v1_combined_ensemble",
                "legacy_cache_key": legacy_key,
                "legacy_manifest_path": str(legacy_manifest_path),
                "legacy_manifest_sha256": sha256_file(legacy_manifest_path),
                "legacy_artifacts": dict(legacy.get("artifacts", {})),
            },
            "artifacts": {
                "migrated_conformer_count": len(conformers),
            },
        }
        atomic_write_json(work_dir / "manifest.json", manifest)
        os.rename(work_dir, result_dir)
        return _crest_ensemble_from_manifest(
            manifest,
            result_dir,
            cache_hit=False,
        )

    def _conformer_thermodynamics(
        self,
        frames: list[XyzFrame],
        crest_log: Path,
    ) -> list[dict[str, float | int]]:
        summary = parse_crest_summary(crest_log)
        if len(summary) == len(frames):
            weights = np.asarray(
                [float(row["boltzmann_weight"]) for row in summary],
                dtype=float,
            )
            if np.isfinite(weights).all() and weights.sum() > 0.0:
                weights /= weights.sum()
                return [
                    {
                        "energy_hartree": float(row["energy_hartree"]),
                        "relative_energy_kcal_mol": float(
                            row["relative_energy_kcal_mol"]
                        ),
                        "boltzmann_weight": float(weight),
                        "degeneracy": int(row["degeneracy"]),
                    }
                    for row, weight in zip(summary, weights, strict=True)
                ]
        if any(frame.energy_hartree is None for frame in frames):
            raise QuantumBackendError(
                "CREST population table is unavailable and XYZ comments lack energies"
            )
        energies = np.asarray(
            [float(frame.energy_hartree) for frame in frames],
            dtype=float,
        )
        relative = (energies - energies.min()) * HARTREE_TO_KCAL_MOL
        raw = np.exp(-relative / (GAS_CONSTANT_KCAL_MOL_K * self.config.temperature_k))
        weights = raw / raw.sum()
        return [
            {
                "energy_hartree": float(frame.energy_hartree),
                "relative_energy_kcal_mol": float(energy),
                "boltzmann_weight": float(weight),
                "degeneracy": 1,
            }
            for frame, energy, weight in zip(
                frames,
                relative,
                weights,
                strict=True,
            )
        ]

    def _run(
        self,
        command: list[str],
        cwd: Path,
        log_path: Path,
        timeout_seconds: float,
    ) -> None:
        environment = self._environment()
        with log_path.open("w", encoding="utf-8") as log:
            try:
                completed = subprocess.run(
                    command,
                    cwd=cwd,
                    env=environment,
                    stdin=subprocess.DEVNULL,
                    stdout=log,
                    stderr=subprocess.STDOUT,
                    timeout=timeout_seconds,
                    check=False,
                )
            except subprocess.TimeoutExpired as exc:
                raise QuantumBackendError(
                    f"calculation exceeded {timeout_seconds:.0f} seconds: "
                    f"{' '.join(command)}"
                ) from exc
        if completed.returncode != 0:
            tail = "\n".join(
                log_path.read_text(
                    encoding="utf-8",
                    errors="replace",
                ).splitlines()[-30:]
            )
            raise QuantumBackendError(
                f"tool exited with status {completed.returncode}: "
                f"{' '.join(command)}\n{tail}"
            )

    def _environment(self) -> dict[str, str]:
        environment = os.environ.copy()
        environment.update(
            {
                "OMP_NUM_THREADS": str(self.config.threads),
                "MKL_NUM_THREADS": str(self.config.threads),
                "OMP_STACKSIZE": "4G",
            }
        )
        xtb_path = Path(self.xtb.path).resolve()
        xtb_root = xtb_path.parent.parent
        share = xtb_root / "share" / "xtb"
        library = xtb_root / "lib"
        if share.is_dir():
            environment["XTBHOME"] = str(share)
        if library.is_dir():
            existing = environment.get("LD_LIBRARY_PATH")
            environment["LD_LIBRARY_PATH"] = (
                f"{library}:{existing}" if existing else str(library)
            )
        environment["PATH"] = f"{xtb_path.parent}:{environment.get('PATH', '')}"
        return environment

    def _tool_info(
        self,
        name: str,
        configured: Path | None,
        environment_name: str,
    ) -> ToolInfo:
        root = Path(__file__).resolve().parent
        candidates = [
            configured,
            Path(os.environ[environment_name])
            if os.environ.get(environment_name)
            else None,
            root / ".stericx" / "tools" / "bin" / name,
            Path(found) if (found := shutil.which(name)) else None,
        ]
        executable = next(
            (
                candidate.resolve()
                for candidate in candidates
                if candidate is not None
                and candidate.is_file()
                and os.access(candidate, os.X_OK)
            ),
            None,
        )
        if executable is None:
            raise QuantumBackendError(
                f"{name} executable not found; run ./install_quantum_tools.sh "
                f"or set {environment_name}"
            )
        probe = subprocess.run(
            [str(executable), "--version"],
            env=self._probe_environment(executable),
            text=True,
            capture_output=True,
            timeout=30.0,
            check=False,
        )
        version_text = "\n".join(
            part.strip() for part in (probe.stdout, probe.stderr) if part.strip()
        )
        if probe.returncode != 0 or not version_text:
            raise QuantumBackendError(
                f"could not determine {name} version from {executable}"
            )
        version_match = re.search(
            r"(?:version|Version)\s*[=:]?\s*v?(\d+(?:\.\d+){1,2})",
            version_text,
        )
        version = (
            version_match.group(1) if version_match else version_text.splitlines()[0]
        )
        return ToolInfo(
            name=name,
            path=str(executable),
            version=version,
            sha256=sha256_file(executable),
        )

    @staticmethod
    def _probe_environment(executable: Path) -> dict[str, str]:
        environment = os.environ.copy()
        root = executable.resolve().parent.parent
        library = root / "lib"
        share = root / "share" / "xtb"
        if library.is_dir():
            environment["LD_LIBRARY_PATH"] = str(library)
        if share.is_dir():
            environment["XTBHOME"] = str(share)
        return environment


def _comment_energy(comment: str) -> float | None:
    if not comment:
        return None
    try:
        value = float(comment.split()[0])
    except (IndexError, ValueError):
        return None
    return value if math.isfinite(value) else None


def _cache_key(payload: dict[str, object]) -> str:
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(encoded).hexdigest()


class CacheLock:
    """Owned cache lock with PID/start-time validation and stale recovery."""

    def __init__(self, path: Path, descriptor: int, token: str, inode: int):
        self.path = path
        self.descriptor = descriptor
        self.token = token
        self.inode = inode
        self._released = False

    @classmethod
    def acquire(cls, path: Path, stale_after_seconds: float) -> CacheLock:
        """Atomically acquire a lock, reclaiming only a demonstrably stale owner."""
        path.parent.mkdir(parents=True, exist_ok=True)
        for _ in range(4):
            token = uuid.uuid4().hex
            try:
                descriptor = os.open(
                    path,
                    os.O_CREAT | os.O_EXCL | os.O_WRONLY,
                    0o600,
                )
            except FileExistsError as exc:
                if _reclaim_stale_lock(path, stale_after_seconds):
                    continue
                owner = _load_json_object(path)
                owner_text = (
                    f" owner_pid={owner.get('pid')} owner_host={owner.get('hostname')}"
                    if owner
                    else ""
                )
                raise QuantumBackendError(
                    f"cache calculation is already locked: {path}{owner_text}"
                ) from exc
            metadata = {
                "schema_version": 1,
                "token": token,
                "pid": os.getpid(),
                "process_start_token": _process_start_token(os.getpid()),
                "hostname": socket.gethostname(),
                "thread_id": threading.get_ident(),
                "created_at_utc": datetime.now(UTC).isoformat(),
                "created_unix_seconds": time.time(),
            }
            try:
                os.write(
                    descriptor,
                    (json.dumps(metadata, sort_keys=True) + "\n").encode(),
                )
                os.fsync(descriptor)
            except BaseException:
                os.close(descriptor)
                path.unlink(missing_ok=True)
                raise
            return cls(
                path,
                descriptor,
                token,
                os.fstat(descriptor).st_ino,
            )
        raise QuantumBackendError(f"could not acquire cache lock: {path}")

    def release(self) -> None:
        """Release only this owner's inode/token, never a replacement lock."""
        if self._released:
            return
        self._released = True
        try:
            try:
                stat = self.path.stat()
                metadata = _load_json_object(self.path)
                if (
                    stat.st_ino == self.inode
                    and metadata is not None
                    and metadata.get("token") == self.token
                ):
                    self.path.unlink(missing_ok=True)
            except FileNotFoundError:
                pass
        finally:
            os.close(self.descriptor)

    def __enter__(self) -> CacheLock:
        return self

    def __exit__(self, *_: object) -> None:
        self.release()


def _process_start_token(pid: int) -> str | None:
    """Return the Linux process start tick used to detect PID reuse."""
    try:
        fields = Path(f"/proc/{pid}/stat").read_text(encoding="utf-8").split()
    except (FileNotFoundError, OSError):
        return None
    return fields[21] if len(fields) > 21 else None


def _pid_is_live(pid: int, expected_start_token: object) -> bool:
    if pid <= 0:
        return False
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        return True
    actual_start = _process_start_token(pid)
    return (
        actual_start is None
        or expected_start_token is None
        or actual_start == str(expected_start_token)
    )


def _reclaim_stale_lock(path: Path, stale_after_seconds: float) -> bool:
    """Remove a stale lock only after verifying its inode did not change."""
    try:
        initial_stat = path.stat()
    except FileNotFoundError:
        return True
    metadata = _load_json_object(path)
    age_seconds = max(0.0, time.time() - initial_stat.st_mtime)
    if metadata is None:
        stale = age_seconds >= min(5.0, stale_after_seconds)
    elif metadata.get("hostname") == socket.gethostname():
        try:
            owner_pid = int(metadata["pid"])
        except (KeyError, TypeError, ValueError):
            stale = age_seconds >= min(5.0, stale_after_seconds)
        else:
            stale = not _pid_is_live(
                owner_pid,
                metadata.get("process_start_token"),
            )
    else:
        stale = age_seconds >= stale_after_seconds
    if not stale:
        return False
    try:
        if path.stat().st_ino != initial_stat.st_ino:
            return False
        path.unlink()
        return True
    except FileNotFoundError:
        return True


def _load_json_object(path: Path) -> dict[str, object] | None:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (FileNotFoundError, json.JSONDecodeError, OSError):
        return None
    return value if isinstance(value, dict) else None


def _load_valid_manifest(path: Path, cache_key: str) -> dict[str, object] | None:
    value = _load_json_object(path)
    if value is None:
        return None
    if (
        value.get("schema_version") != CACHE_SCHEMA_VERSION
        or value.get("cache_key") != cache_key
        or value.get("status") != "complete"
    ):
        return None
    return value


def _lmo_result_from_manifest(
    manifest: dict[str, object],
    manifest_path: Path,
    cache_hit: bool,
) -> LmoCenterResult:
    center = manifest["coordination_center_angstrom"]
    selected = manifest["selected_lmo_angstrom"]
    neighbors = manifest["donor_neighbors"]
    if not (
        isinstance(center, list)
        and isinstance(selected, list)
        and isinstance(neighbors, list)
    ):
        raise QuantumBackendError(f"invalid cached LMO manifest: {manifest_path}")
    return LmoCenterResult(
        center_angstrom=tuple(float(value) for value in center),  # type: ignore[arg-type]
        selected_lmo_angstrom=tuple(float(value) for value in selected),  # type: ignore[arg-type]
        selected_lmo_index=int(manifest["selected_lmo_index"]),
        donor_neighbors=tuple(int(value) for value in neighbors),  # type: ignore[arg-type]
        cache_key=str(manifest["cache_key"]),
        cache_hit=cache_hit,
        manifest_path=str(manifest_path),
    )


def _crest_ensemble_from_manifest(
    manifest: dict[str, object],
    result_dir: Path,
    cache_hit: bool,
) -> CrestEnsemble:
    rows = manifest.get("conformers")
    if not isinstance(rows, list):
        raise QuantumBackendError(f"invalid CREST manifest: {result_dir}")
    conformers: list[CrestConformer] = []
    for row in rows:
        if not isinstance(row, dict):
            raise QuantumBackendError(f"invalid CREST conformer: {result_dir}")
        relative_path = Path(str(row["xyz_path"]))
        xyz_path = result_dir / relative_path
        if not xyz_path.is_file():
            raise QuantumBackendError(f"cached conformer is missing: {xyz_path}")
        expected_sha = str(row["xyz_sha256"])
        if sha256_file(xyz_path) != expected_sha:
            raise QuantumBackendError(f"cached conformer hash mismatch: {xyz_path}")
        conformers.append(
            CrestConformer(
                index=int(row["index"]),
                xyz_path=str(xyz_path),
                xyz_sha256=expected_sha,
                energy_hartree=float(row["energy_hartree"]),
                relative_energy_kcal_mol=float(row["relative_energy_kcal_mol"]),
                boltzmann_weight=float(row["boltzmann_weight"]),
                degeneracy=int(row["degeneracy"]),
            )
        )
    if not conformers:
        raise QuantumBackendError(f"cached CREST ensemble is empty: {result_dir}")
    return CrestEnsemble(
        cache_key=str(manifest["cache_key"]),
        cache_hit=cache_hit,
        manifest_path=str(result_dir / "manifest.json"),
        conformers=tuple(conformers),
    )


def _ensemble_from_job_manifest(
    manifest: dict[str, object],
    manifest_path: Path,
    xtb: ToolInfo,
    crest: ToolInfo,
    cache_hit: bool,
) -> QuantumEnsemble:
    rows = manifest.get("conformers")
    if not isinstance(rows, list):
        raise QuantumBackendError(f"invalid quantum job manifest: {manifest_path}")
    conformers: list[QuantumConformer] = []
    for row in rows:
        if not isinstance(row, dict):
            raise QuantumBackendError(
                f"invalid quantum conformer manifest: {manifest_path}"
            )
        xyz_path = Path(str(row["xyz_path"]))
        center = row.get("coordination_center_angstrom")
        if not xyz_path.is_file():
            raise QuantumBackendError(f"cached conformer is missing: {xyz_path}")
        if not isinstance(center, list) or len(center) != 3:
            raise QuantumBackendError(f"invalid coordination center: {manifest_path}")
        center_values = tuple(float(value) for value in center)
        if not all(math.isfinite(value) for value in center_values):
            raise QuantumBackendError(
                f"non-finite coordination center: {manifest_path}"
            )
        conformers.append(
            QuantumConformer(
                index=int(row["index"]),
                xyz_path=str(xyz_path),
                energy_hartree=float(row["energy_hartree"]),
                relative_energy_kcal_mol=float(row["relative_energy_kcal_mol"]),
                boltzmann_weight=float(row["boltzmann_weight"]),
                degeneracy=int(row["degeneracy"]),
                coordination_center_angstrom=center_values,  # type: ignore[arg-type]
                lmo_cache_key=str(row["lmo_cache_key"]),
            )
        )
    if not conformers:
        raise QuantumBackendError(f"cached quantum ensemble is empty: {manifest_path}")
    return QuantumEnsemble(
        cache_key=str(manifest["cache_key"]),
        cache_hit=cache_hit,
        manifest_path=str(manifest_path),
        conformers=tuple(conformers),
        xtb=xtb,
        crest=crest,
        crest_cache_key=str(manifest["crest_cache_key"]),
        crest_cache_hit=True if cache_hit else bool(manifest["crest_cache_hit"]),
        lmo_cache_hits=(
            len({conformer.lmo_cache_key for conformer in conformers})
            if cache_hit
            else int(manifest["lmo_cache_hits"])
        ),
        lmo_cache_misses=0 if cache_hit else int(manifest["lmo_cache_misses"]),
        job_resumed=False if cache_hit else bool(manifest["job_resumed"]),
        effective_lmo_workers=int(manifest["effective_lmo_workers"]),
    )
