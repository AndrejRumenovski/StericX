//! Allocation diagnostics for one complete six-conformer API call.
#[cfg(feature = "profiling")]
#[global_allocator]
static ALLOCATOR: steric_x::profiling::TrackingAllocator = steric_x::profiling::TrackingAllocator;

#[cfg(feature = "profiling")]
fn main() -> Result<(), Box<dyn std::error::Error>> {
    use std::{fs, hint::black_box};
    use steric_x::ensemble::{DescriptorOptions, Ensemble, StatisticsOptions};
    use steric_x::profiling::allocation_snapshot;

    let args: Vec<String> = std::env::args().collect();
    if args.len() != 3 {
        return Err("usage: paton_ensemble_allocations INPUT OPTIONS_JSON".into());
    }
    let ensemble = Ensemble::from_file(&args[1])?;
    let options: DescriptorOptions = serde_json::from_str(&fs::read_to_string(&args[2])?)?;
    let run = || {
        let descriptors = ensemble.descriptors(&options);
        let statistics = descriptors
            .statistics(&StatisticsOptions::default(), None)
            .expect("benchmark conformers must have valid statistics");
        black_box((descriptors, statistics));
    };
    run();
    let before = allocation_snapshot();
    run();
    let after = allocation_snapshot();
    println!(
        "{}",
        serde_json::json!({
            "conformers": ensemble.conformers.len(),
            "measured_calls": 1,
            "allocation_calls": after.allocation_calls - before.allocation_calls,
            "reallocation_calls": after.reallocation_calls - before.reallocation_calls,
            "deallocation_calls": after.deallocation_calls - before.deallocation_calls,
            "requested_bytes": after.requested_bytes - before.requested_bytes,
            "live_requested_bytes_delta": after.live_requested_bytes as i128 - before.live_requested_bytes as i128,
            "failed_allocation_calls": after.failed_allocation_calls - before.failed_allocation_calls,
            "process_peak_live_requested_bytes": after.peak_live_requested_bytes,
            "scope": "One warmed descriptor/statistics call including result destruction; excludes input parsing and JSON encoding. Peak is process-wide through snapshot."
        })
    );
    Ok(())
}

#[cfg(not(feature = "profiling"))]
fn main() {
    eprintln!("Build this diagnostic with --features profiling");
    std::process::exit(2);
}
