//! Source-bound observer of the explicit public certified-volume API.
use serde_json::{Value, json};
use steric_x::geometry::{AtomicSphere, CertifiedVolumeOptions, sphere_union_volume_certified};
#[cfg(feature = "profiling")]
#[global_allocator]
static ALLOCATOR: steric_x::profiling::TrackingAllocator = steric_x::profiling::TrackingAllocator;
fn main() -> Result<(), Box<dyn std::error::Error>> {
    let requests: Vec<Value> = serde_json::from_str(&std::fs::read_to_string(
        std::env::args().nth(1).ok_or("request path required")?,
    )?)?;
    for q in requests {
        let region: AtomicSphere = serde_json::from_value(q["region"].clone())?;
        let atoms: Vec<AtomicSphere> = serde_json::from_value(q["atoms"].clone())?;
        let options = CertifiedVolumeOptions {
            absolute_tolerance: q["tolerance"].as_f64().unwrap(),
            max_cells: q["budget"].as_u64().unwrap() as usize,
        };
        #[cfg(feature = "profiling")]
        let before = steric_x::profiling::allocation_snapshot();
        let result = sphere_union_volume_certified(region, &atoms, options);
        #[cfg(feature = "profiling")]
        let allocation = {
            let after = steric_x::profiling::allocation_snapshot();
            json!({"calls":after.allocation_calls-before.allocation_calls,
                "reallocations":after.reallocation_calls-before.reallocation_calls,
                "requested_bytes":after.requested_bytes-before.requested_bytes,
                "live_bytes_delta":after.live_requested_bytes as i128-before.live_requested_bytes as i128,
                "process_peak_live_bytes":after.peak_live_requested_bytes})
        };
        #[cfg(not(feature = "profiling"))]
        let allocation = Value::Null;
        match result {
            Err(e) => println!("{}", json!({"id":q["id"],"error":e.to_string()})),
            Ok(r) => println!("{}", json!({"id":q["id"],"volume":r.volume_angstrom3,
                "lower":r.lower_bound_angstrom3,"upper":r.upper_bound_angstrom3,
                "absolute_error_bound":r.absolute_error_bound_angstrom3,"relative_error_bound":r.relative_error_bound,
                "target_met":r.target_met,"status":r.status,"cells":r.cells_evaluated,
                "subdivisions":r.subdivisions,"geometric_evaluations":r.geometric_evaluations,
                "smooth_cells":r.smooth_cells,"range_cells":r.range_cells,"empty_cells":r.empty_cells,
                "peak_heap_cells":r.peak_heap_cells,"elapsed_ns":r.elapsed_seconds*1e9,
                "analytical":r.analytical,"allocations":allocation})),
        }
    }
    Ok(())
}
