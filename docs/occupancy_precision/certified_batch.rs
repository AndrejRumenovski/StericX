//! Parallel batches of the frozen serial certified API; ordered exact results.
use rayon::prelude::*;
use serde::Deserialize;
use serde_json::json;
use steric_x::geometry::{sphere_union_volume_certified, AtomicSphere, CertifiedVolumeOptions};

#[derive(Deserialize)]
struct Request {
    id: String,
    region: AtomicSphere,
    atoms: Vec<AtomicSphere>,
    tolerance: f64,
    budget: usize,
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<String> = std::env::args().collect();
    let requests: Vec<Request> = serde_json::from_str(&std::fs::read_to_string(&args[1])?)?;
    let threads: usize = args[2].parse()?;
    let pool = rayon::ThreadPoolBuilder::new().num_threads(threads).build()?;
    let calculate = || {
        requests.par_iter().map(|q| {
            sphere_union_volume_certified(q.region, &q.atoms, CertifiedVolumeOptions {
                absolute_tolerance: q.tolerance, max_cells: q.budget,
            }).expect("validated certified batch request")
        }).collect::<Vec<_>>()
    };
    std::hint::black_box(pool.install(calculate));
    let start = std::time::Instant::now();
    let result = pool.install(calculate);
    let elapsed_ns = start.elapsed().as_nanos();
    let rows: Vec<_> = requests.iter().zip(result).map(|(q, r)| {
        let mut value = serde_json::to_value(r).unwrap();
        value.as_object_mut().unwrap().remove("elapsed_seconds");
        json!({"id": q.id, "result": value})
    }).collect();
    println!("{}", json!({"threads":threads,"elapsed_ns":elapsed_ns,"rows":rows}));
    Ok(())
}
