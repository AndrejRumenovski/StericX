use serde_json::{Value, json};
use steric_x::geometry::{
    AtomicSphere, CertifiedVolumeOptions, sphere_union_volume_certified_geometric,
    sphere_union_volume_certified_radius_profile, sphere_union_volume_certified_regions,
};
#[cfg(feature = "profiling")]
#[global_allocator]
static ALLOCATOR: steric_x::profiling::TrackingAllocator = steric_x::profiling::TrackingAllocator;
fn main() -> Result<(), Box<dyn std::error::Error>> {
    let qs: Vec<Value> = serde_json::from_str(&std::fs::read_to_string(
        std::env::args().nth(1).ok_or("request file")?,
    )?)?;
    for q in qs {
        let center: [f64; 3] = serde_json::from_value(q["center"].clone())?;
        let atoms: Vec<AtomicSphere> = serde_json::from_value(q["atoms"].clone())?;
        let radii: Vec<f64> = serde_json::from_value(q["radii"].clone())?;
        let opts = CertifiedVolumeOptions {
            absolute_tolerance: q["tolerance"].as_f64().unwrap(),
            max_cells: q["budget"].as_u64().unwrap() as usize,
        };
        #[cfg(feature = "profiling")]
        let before = steric_x::profiling::allocation_snapshot();
        let profile = sphere_union_volume_certified_radius_profile(center, &atoms, &radii, opts)?;
        #[cfg(feature = "profiling")]
        let allocations = {
            let after = steric_x::profiling::allocation_snapshot();
            json!({"calls":after.allocation_calls-before.allocation_calls,"reallocations":after.reallocation_calls-before.reallocation_calls,"requested_bytes":after.requested_bytes-before.requested_bytes,"live_bytes_delta":after.live_requested_bytes as i128-before.live_requested_bytes as i128,"process_peak_live_bytes":after.peak_live_requested_bytes})
        };
        #[cfg(not(feature = "profiling"))]
        let allocations = Value::Null;
        let mut controls = Vec::new();
        if q["controls"].as_bool().unwrap_or(false) {
            for &radius in &radii {
                let sphere = AtomicSphere { center, radius };
                controls.push(match sphere_union_volume_certified_regions(sphere,&atoms,opts) {
     Ok(r)=>json!({"regions":r,"total":sphere_union_volume_certified_geometric(sphere,&atoms,opts)?}),
     Err(e)=>json!({"error":e.to_string()}),
    });
            }
        }
        println!(
            "{}",
            json!({"id":q["id"],"profile":profile,"controls":controls,"allocations":allocations})
        );
    }
    Ok(())
}
