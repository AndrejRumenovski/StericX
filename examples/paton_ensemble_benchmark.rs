//! In-memory ensemble timing adapter; no independent chemistry implementation.
use std::{error::Error, fs, hint::black_box, time::Instant};
use steric_x::ensemble::{DescriptorOptions, DescriptorOutcome, Ensemble, StatisticsOptions};

fn main() -> Result<(), Box<dyn Error>> {
    let args: Vec<String> = std::env::args().collect();
    if args.len() != 4 {
        return Err("usage: paton_ensemble_benchmark INPUT OPTIONS_JSON ITERATIONS".into());
    }
    let ensemble = Ensemble::from_file(&args[1])?;
    let options: DescriptorOptions = serde_json::from_str(&fs::read_to_string(&args[2])?)?;
    let iterations: usize = args[3].parse()?;
    if iterations == 0 {
        return Err("iterations must be positive".into());
    }
    let run = || {
        let descriptors = ensemble.descriptors(&options);
        let statistics: DescriptorOutcome<_> = descriptors
            .statistics(&StatisticsOptions::default(), None)
            .into();
        (descriptors, statistics)
    };
    black_box(run());
    let started = Instant::now();
    for _ in 0..iterations {
        black_box(run());
    }
    let elapsed_ns = started.elapsed().as_nanos();
    let (descriptors, statistics) = run();
    // Widen native f32 values exactly, matching the binding transport contract.
    let output = serde_json::json!({
        "elapsed_ns": elapsed_ns,
        "iterations": iterations,
        "result": {
            "descriptors": descriptors,
            "statistics": statistics,
            "sensitivity": null,
            "sensitivity_csv": null
        }
    });
    println!("{}", serde_json::to_string(&output)?);
    Ok(())
}
