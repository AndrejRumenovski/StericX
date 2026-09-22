// End-to-end wrapper around the UNCHANGED public buried-volume API.
// No memoization, pre-parsing, numerical changes or private copied kernels.
use rayon::prelude::*;
use serde_json::{Value,json};
use steric_x::{Molecule, BuriedVolumeConfig, BuriedVolumeCalculator, bonded_neighbors};
fn calculate(path: &String) -> Value {
    let m = Molecule::from_xyz_file(std::path::Path::new(path)).unwrap();
    let donors: Vec<_> = m.atoms.iter().enumerate().filter(|(_,a)|a.element=="P").map(|(i,_)|i).collect();
    assert_eq!(donors.len(),1);
    let d=donors[0];
    let ns: Vec<_> = bonded_neighbors(&m,d).iter().map(|(_,i)|*i).collect();
    let cfg=BuriedVolumeConfig {center_distance:2.28,..Default::default()};
    let b=BuriedVolumeCalculator::compute_with_neighbors(&m,d,ns.try_into().unwrap(),cfg).unwrap();
    json!({"file":path,"buried_volume":b.buried_volume,"percent_buried_volume":b.percent_buried_volume,
           "qvbur_min":b.qvbur_min,"qvbur_max":b.qvbur_max,"max_delta_qvbur":b.max_delta_qvbur,
           "ovbur_min":b.ovbur_min,"ovbur_max":b.ovbur_max,"near_vbur":b.near_vbur,"far_vbur":b.far_vbur})
}
fn main() {
    let paths: Vec<_> = std::env::args().skip(1).collect();
    let rows: Vec<_> = paths.par_iter().map(calculate).collect();
    serde_json::to_writer_pretty(std::io::stdout().lock(), &rows).unwrap();
    println!();
}
