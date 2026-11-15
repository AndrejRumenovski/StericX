//! Public-API volume workflow with separately reported warmed loop timing.
use rayon::prelude::*;
use serde_json::{Value, json};
use std::{hint::black_box, path::Path, time::Instant};
use steric_x::{BuriedVolumeCalculator, BuriedVolumeConfig, Molecule, bonded_neighbors};

fn calculate(path: &String) -> Value {
    let molecule = Molecule::from_xyz_file(Path::new(path)).unwrap();
    let donors: Vec<_> = molecule
        .atoms
        .iter()
        .enumerate()
        .filter(|(_, atom)| atom.element == "P")
        .map(|(i, _)| i)
        .collect();
    assert_eq!(donors.len(), 1);
    let donor = donors[0];
    let neighbors: Vec<_> = bonded_neighbors(&molecule, donor)
        .iter()
        .map(|(_, i)| *i)
        .collect();
    let config = BuriedVolumeConfig {
        center_distance: 2.28,
        ..Default::default()
    };
    let b = BuriedVolumeCalculator::compute_with_neighbors(
        &molecule,
        donor,
        neighbors.try_into().unwrap(),
        config,
    )
    .unwrap();
    json!({"file":path,"buried_volume":b.buried_volume,
        "percent_buried_volume":b.percent_buried_volume,"qvbur_min":b.qvbur_min,
        "qvbur_max":b.qvbur_max,"max_delta_qvbur":b.max_delta_qvbur,
        "ovbur_min":b.ovbur_min,"ovbur_max":b.ovbur_max,
        "near_vbur":b.near_vbur,"far_vbur":b.far_vbur})
}
fn main() {
    let mut args = std::env::args().skip(1);
    let iterations: usize = args.next().expect("warm iterations").parse().unwrap();
    let paths: Vec<_> = args.collect();
    let run = || paths.par_iter().map(calculate).collect::<Vec<_>>();
    let rows = run();
    let mut samples = Vec::with_capacity(iterations);
    for _ in 0..iterations {
        let start = Instant::now();
        let observed = black_box(run());
        samples.push(start.elapsed().as_nanos());
        assert_eq!(observed, rows);
    }
    serde_json::to_writer(
        std::io::stdout().lock(),
        &json!({"rows":rows,"warm_elapsed_ns":samples}),
    )
    .unwrap();
    println!();
}
