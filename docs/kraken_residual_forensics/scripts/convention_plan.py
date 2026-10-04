"""Seal independently sourced conventions BEFORE new effect calculations."""
import datetime
import importlib.metadata
import inspect
from pathlib import Path

from morfeus import BuriedVolume, Sterimol
from morfeus.utils import get_radii
import morfeus.buried_volume
import morfeus.sterimol
import morfeus.data

from freeze import OUT, save, sha


def main():
    out = OUT / "conventions"
    out.mkdir(exist_ok=False)
    primary = "https://github.com/the-matter-lab/kraken/blob/4eaad505c1343e6083032b4a3fda47e004e19734/conf_selection_and_DFT/PL_dft_library_201027.py"
    sources = [{"url": primary, "local": "frozen/kraken_dft_source.py.txt", "sha256": sha(OUT / "frozen/kraken_dft_source.py.txt")}]
    for module in (morfeus.buried_volume, morfeus.sterimol, morfeus.data):
        path = Path(module.__file__)
        target = out / (path.stem + ".py.txt")
        target.write_bytes(path.read_bytes())
        sources.append({"installed_path": str(path), "local": str(target.relative_to(OUT)), "sha256": sha(target)})
    elements = sorted({a["element"] for r in __import__("analyze").read_jsonl(OUT / "frozen/requests.jsonl.gz") for a in r["atoms"]})
    radii = dict(zip(elements, map(float, get_radii(elements, radii_type="bondi"))))
    save(out / "plan.json", {
        "created_utc": datetime.datetime.now(datetime.UTC).isoformat(),
        "sources_sealed_before_effect_calculations": True,
        "no_target_fitting": True, "sources": sources,
        "morfeus_version": importlib.metadata.version("morfeus-ml"),
        "morfeus_BV_signature": str(inspect.signature(BuriedVolume)),
        "morfeus_Sterimol_signature": str(inspect.signature(Sterimol)),
        "constant_audit": [
            {"quantity": "virtual_metal_distance_A", "value": 2.28, "source": "original add_valence distpx[Pd], line 66"},
            {"quantity": "center", "value": "P + 2.28 normalize(sum(P-neighbor))", "source": "original add_valence lines 73-78; raw vectors, no individual normalization"},
            {"quantity": "donor_connectivity", "value": "original Pyykko-Atsumi rcov and damp>0.85", "source": "original get_conmat lines 24-60; not fitted"},
            {"quantity": "sphere_radius_A", "value": 3.5, "source": "original implicit Morfeus default; SI audit and captured Morfeus signature", "historical_dependency_pin": "unrecovered"},
            {"quantity": "BV_radii_scale", "value": 1.17, "source": "original implicit Morfeus default; captured signature/SI audit"},
            {"quantity": "BV_include_hydrogens", "value": False, "source": "original implicit Morfeus default; captured signature/SI audit; donor H still define frame"},
            {"quantity": "BV_density_A3", "value": 0.001, "source": "original implicit Morfeus default; captured signature/SI audit", "historical_dependency_pin": "unrecovered"},
            {"quantity": "BV_radii", "value": radii, "source": "captured Morfeus Bondi table and explicit 2.0 fallback where absent; compatibility table not physical ground truth"},
            {"quantity": "Sterimol_radii", "value": "Bondi, every radius equal to 1.20 replaced by 1.09", "source": "original morfeus_properties lines 245-252"},
            {"quantity": "Sterimol_L_correction_A", "value": 0.4, "source": "Morfeus Sterimol.calculate; documented Verloop convention"},
            {"quantity": "Sterimol_directions", "value": 3600, "source": "original morfeus_properties line 252; contemporary np.linspace endpoint included; historical phase/version unknown"},
            {"quantity": "Sterimol_axis", "value": "virtual Pd to P", "source": "original add_valence reordering and Sterimol(dummy=2,attached=1)"},
            {"quantity": "BV_frame", "value": "center at Pd; P toward -z; each of three P-neighbors defines +xz half-plane", "source": "original morfeus_properties lines 179-193; Morfeus coordinate frame"},
            {"quantity": "quadrant_octant_reduction", "value": "all 3 frames extrema; maximum adjacent quadrant difference", "source": "original morfeus_properties"},
            {"quantity": "total_near_far_reduction", "value": "last neighbor in original atom ordering", "source": "original loop overwrites scalars; this finite-grid order dependence is not adopted as a general StericX default"},
            {"quantity": "ensemble_reduction", "value": "complete ensemble min/max/max-minus-min; vburminconf selected by minimum total Vbur", "source": "original read_ligand lines 553-567"},
            {"quantity": "Boltzmann", "value": "historical g_tz_gas energies, 298.15 K, R=0.0019872036", "source": "original lines 30-33,533-537", "status": "unavailable original per-conformer energies; not inferred"},
        ],
        "fixed_stages": [
            {"name": "baseline_replay", "change": "none; verify frozen current outputs"},
            {"name": "raw_center", "change": "original primary donor connectivity, raw-vector center only; current radii/grid"},
            {"name": "primary_radii", "change": "also original Morfeus Bondi table for BV, Paton substitution for Sterimol; default grid retained"},
            {"name": "primary_grid", "change": "also documented 0.001 BV density; preserve public mean and separately extract historical last-frame scalar"},
        ],
        "preserved_limits": "Native Sterimol remains a 360-direction f32 scan. A later optional sampling change requires a separately recorded mathematical witness and validation. No production default changes are authorized by an improved R² alone.",
        "convergence_densities_A3": [0.1, 0.01, 0.001, 0.000125, 0.000015625],
        "convergence_selection": "all conformers of 3 largest baseline headline ligands plus minimum-ID ligand per donor environment; further final large residuals may be added with selection recorded before calculating finer results",
        "convergence_interpretation": "Successive differences provide an empirical envelope, not a rigorous bound or theoretical ceiling; no grid selected against Kraken values.",
    })
    print(sha(out / "plan.json"))


if __name__ == "__main__":
    main()
