//! Measurement adapter: calls public Rust APIs; contains no descriptor oracle.
use glam::Vec3;
use rayon::prelude::*;
use serde::Serialize;
use serde_json::json;
use std::{error::Error, fs, hint::black_box, time::Instant};
use steric_x::density::{DensityDeclaration, DensityField, DensityGrid, DensityUnit};
use steric_x::ensemble::{DescriptorOptions, Ensemble, GeometrySnapshot, StatisticsOptions};
use steric_x::geometry::{
    BuriedVolumeCalculator, PyramidalizationCalculator, SterimolCalculator,
    coordination_center_with_neighbors, detect_trivalent_donor,
};

#[cfg(feature = "profiling")]
#[global_allocator]
static ALLOCATOR: steric_x::profiling::TrackingAllocator = steric_x::profiling::TrackingAllocator;

fn measure<T: Serialize>(kind: &str, iterations: usize, run: impl Fn() -> T) {
    // The timed value is destroyed each iteration, and serialization is outside
    // the measurement. The separate serialization lane explicitly includes it.
    black_box(run());
    #[cfg(feature = "profiling")]
    let before = steric_x::profiling::allocation_snapshot();
    let start = Instant::now();
    for _ in 0..iterations {
        black_box(run());
    }
    let elapsed_ns = start.elapsed().as_nanos();
    #[cfg(feature = "profiling")]
    let allocations = {
        let after = steric_x::profiling::allocation_snapshot();
        json!({
            "allocation_calls": after.allocation_calls - before.allocation_calls,
            "reallocation_calls": after.reallocation_calls - before.reallocation_calls,
            "deallocation_calls": after.deallocation_calls - before.deallocation_calls,
            "requested_bytes": after.requested_bytes - before.requested_bytes,
            "live_bytes_delta": after.live_requested_bytes as i128 - before.live_requested_bytes as i128,
            "failed_calls": after.failed_allocation_calls - before.failed_allocation_calls,
            "process_peak_live_bytes": after.peak_live_requested_bytes,
        })
    };
    #[cfg(not(feature = "profiling"))]
    let allocations = serde_json::Value::Null;
    println!("{}", json!({"kind":kind, "iterations":iterations,
        "elapsed_ns":elapsed_ns,"allocations":allocations,"result":run()}));
}

fn declaration() -> DensityDeclaration {
    DensityDeclaration {
        field: DensityField::TotalElectronDensity,
        unit: DensityUnit::ElectronsPerBohrCubed,
        provenance: "Frozen campaign input; explicit total-density declaration".into(),
    }
}

fn main() -> Result<(), Box<dyn Error>> {
    let args: Vec<_> = std::env::args().collect();
    if args.len() < 5 {
        return Err("usage: observer KIND INPUT OPTIONS_JSON ITERATIONS [THREADS]".into());
    }
    let kind = args[1].as_str();
    let iterations: usize = args[4].parse()?;
    if iterations == 0 {
        return Err("iterations must be positive".into());
    }
    let text = fs::read_to_string(&args[2])?;
    if kind.starts_with("density") {
        let grid = DensityGrid::from_cube_str("frozen.cube", &text, declaration())?;
        let center = grid.metadata().atoms[0].position_angstrom;
        match kind {
            "density" => measure(kind, iterations, || grid.buried_volume(center, 1.5, 0.0016).unwrap()),
            "density_profile" => measure(kind, iterations, || grid.buried_volume_scan(center, &[1.0, 1.5, 2.0], 0.0016).unwrap()),
            "density_parse" => measure(kind, iterations, || {
                let g = DensityGrid::from_cube_str("frozen.cube", &text, declaration()).unwrap();
                (g.metadata().dimensions, black_box(g.values()).len())
            }),
            _ => return Err("unknown density lane".into()),
        }
        return Ok(());
    }
    let ensemble = Ensemble::from_file(&args[2])?;
    let options: DescriptorOptions = serde_json::from_str(&fs::read_to_string(&args[3])?)?;
    if kind == "parse" {
        measure(kind, iterations, || {
            let e = if args[2].ends_with(".xyz") {
                Ensemble::from_xyz_str(&text, "frozen.xyz")
            } else {
                Ensemble::from_sdf_str(&text, "frozen.sdf")
            }.unwrap();
            // Touch every parsed geometry without timing a JSON conversion.
            e.conformers.iter().map(|c| black_box(c.frame.as_ref().unwrap().molecule.atoms.len())).sum::<usize>()
        });
        return Ok(());
    }
    let c = &ensemble.conformers[0];
    let frame = c.frame.as_ref().map_err(Clone::clone)?;
    let mol = &frame.molecule;
    let topology = detect_trivalent_donor(mol, frame.bonds.as_deref(), &options.donor_element,
        options.donor_index, None, true)?;
    let donor = topology.donor_idx;
    let neighbors = topology.substituents;
    let center = match options.center {
        Some(c) => Vec3::from_array(c),
        None => coordination_center_with_neighbors(mol, donor, neighbors, options.volume)?,
    };
    match kind {
        "sterimol" => measure(kind, iterations, || SterimolCalculator::compute_with_dummy(mol, donor, center).unwrap()),
        "volume" => measure(kind, iterations, || BuriedVolumeCalculator::compute_with_center_and_neighbors(mol, donor, neighbors, center, options.volume).unwrap()),
        "regions" => measure(kind, iterations, || BuriedVolumeCalculator::compute_detailed_with_center_and_neighbors(mol, donor, neighbors, center, options.volume).unwrap()),
        "pyramidalization" => measure(kind, iterations, || PyramidalizationCalculator::compute(mol, donor, neighbors).unwrap()),
        "center" => measure(kind, iterations, || coordination_center_with_neighbors(mol, donor, neighbors, options.volume).unwrap().to_array()),
        "snapshot" => measure(kind, iterations, || GeometrySnapshot::from(frame)),
        "radius_profile" => measure(kind, iterations, || c.buried_volume_scan(&[2.5, 3.5, 4.5], &options).unwrap()),
        "ensemble" => measure(kind, iterations, || ensemble.descriptors(&options)),
        "statistics" => {
            let descriptors = ensemble.descriptors(&options);
            measure(kind, iterations, || descriptors.statistics(&StatisticsOptions::default(), None).unwrap());
        }
        "serialize" => {
            let descriptors = ensemble.descriptors(&options);
            measure(kind, iterations, || serde_json::to_string(&descriptors).unwrap());
        }
        "batch" => {
            let threads = args.get(5).map(|s| s.parse::<usize>()).transpose()?.unwrap_or(1);
            let pool = rayon::ThreadPoolBuilder::new().num_threads(threads).build()?;
            measure(kind, iterations, || pool.install(|| ensemble.conformers.par_iter()
                .map(|c| c.descriptors(&options)).collect::<Vec<_>>()));
        }
        _ => return Err("unknown measurement lane".into()),
    }
    Ok(())
}
