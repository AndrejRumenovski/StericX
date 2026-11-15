//! Calls public continuous-volume APIs with externally specified inputs.
use serde::Deserialize;
use serde_json::json;
use std::{fs, hint::black_box, time::Instant};
use steric_x::geometry::{AtomicSphere, ContinuousVolumeMethod, sphere_union_volume};

#[cfg(feature = "profiling")]
#[global_allocator]
static ALLOCATOR: steric_x::profiling::TrackingAllocator = steric_x::profiling::TrackingAllocator;

#[derive(Deserialize)]
struct Case {
    id: String,
    region: AtomicSphere,
    atoms: Vec<AtomicSphere>,
    method: ContinuousVolumeMethod,
    samples: usize,
    #[serde(default = "one")]
    iterations: usize,
}
fn one() -> usize {
    1
}
fn main() -> Result<(), Box<dyn std::error::Error>> {
    let path = std::env::args().nth(1).ok_or("case JSON path required")?;
    let cases: Vec<Case> = serde_json::from_str(&fs::read_to_string(path)?)?;
    for case in cases {
        assert!(case.iterations > 0);
        black_box(sphere_union_volume(case.region, &case.atoms, case.method)?);
        let mut times = Vec::with_capacity(case.samples);
        let mut allocations = Vec::with_capacity(case.samples);
        for _ in 0..case.samples {
            #[cfg(feature = "profiling")]
            let before = steric_x::profiling::allocation_snapshot();
            let start = Instant::now();
            for _ in 0..case.iterations {
                black_box(sphere_union_volume(
                    case.region,
                    black_box(&case.atoms),
                    case.method,
                )?);
            }
            let elapsed = start.elapsed().as_nanos();
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
            let allocation = serde_json::Value::Null;
            times.push(elapsed);
            allocations.push(allocation);
        }
        let result = sphere_union_volume(case.region, &case.atoms, case.method)?;
        println!(
            "{}",
            json!({"id":case.id,"iterations":case.iterations,"elapsed_ns":times,"allocations":allocations,"result":result})
        );
    }
    Ok(())
}
