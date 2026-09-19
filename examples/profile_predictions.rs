//! Export every prediction as native-endian f32 bytes for profiling equivalence checks.
//! This does not replace the user-facing `predict` workload or alter inference.

use std::error::Error;
use std::fs;
use std::path::Path;
use steric_x::{RegressXPredictor, SigPackReader};

#[cfg(feature = "profiling")]
#[global_allocator]
static PROFILE_ALLOCATOR: steric_x::profiling::TrackingAllocator =
    steric_x::profiling::TrackingAllocator;

fn main() -> Result<(), Box<dyn Error>> {
    #[cfg(feature = "profiling")]
    let profile = steric_x::profiling::Session::from_env();
    let result = export_predictions();
    #[cfg(feature = "profiling")]
    if let Some(profile) = profile {
        profile.finish()?;
    }
    result
}

fn export_predictions() -> Result<(), Box<dyn Error>> {
    steric_x::profile_scope!("orchestration", "profile_predictions::export");
    let args: Vec<_> = std::env::args_os().skip(1).collect();
    if args.len() != 3 {
        return Err("usage: profile_predictions DATA.sigpack WEIGHTS.json OUTPUT.f32".into());
    }
    let document: serde_json::Value = serde_json::from_slice(&fs::read(&args[1])?)?;
    let weights: [f32; 8] =
        serde_json::from_value(document.get("weights").unwrap_or(&document).clone())?;
    if weights.iter().any(|weight| !weight.is_finite()) {
        return Err("weights must be finite".into());
    }
    let records = SigPackReader::open(Path::new(&args[0]))?;
    let predictions = RegressXPredictor::new(weights).predict_batch(records.records());
    fs::write(&args[2], bytemuck::cast_slice(&predictions))?;
    Ok(())
}
