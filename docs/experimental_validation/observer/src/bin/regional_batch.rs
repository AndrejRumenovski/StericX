use rayon::prelude::*;
use serde_json::{Value, json};
use steric_x::geometry::{
    AtomicSphere, CertifiedVolumeOptions, sphere_union_volume_certified_regions,
};
fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<_> = std::env::args().collect();
    let qs: Vec<Value> = serde_json::from_str(&std::fs::read_to_string(&args[1])?)?;
    let threads: usize = args[2].parse()?;
    let cases: Vec<_> = qs
        .into_iter()
        .map(|q| {
            Ok::<_, Box<dyn std::error::Error>>((
                q["id"].clone(),
                serde_json::from_value::<AtomicSphere>(q["region"].clone())?,
                serde_json::from_value::<Vec<AtomicSphere>>(q["atoms"].clone())?,
                CertifiedVolumeOptions {
                    absolute_tolerance: q["tolerance"].as_f64().unwrap(),
                    max_cells: q["budget"].as_u64().unwrap() as usize,
                },
            ))
        })
        .collect::<Result<_, _>>()?;
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build()?;
    let start = std::time::Instant::now();
    let rows=pool.install(||cases.par_iter().map(|(id,region,atoms,opts)|json!({"id":id,"result":sphere_union_volume_certified_regions(*region,atoms,*opts).unwrap()})).collect::<Vec<_>>());
    println!(
        "{}",
        json!({"threads":threads,"elapsed_seconds":start.elapsed().as_secs_f64(),"rows":rows})
    );
    Ok(())
}
