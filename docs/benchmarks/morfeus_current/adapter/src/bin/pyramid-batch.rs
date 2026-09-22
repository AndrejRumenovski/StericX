// Separate end-to-end pyramidalization API driver; no kernel changes.
use rayon::prelude::*;
use serde_json::json;
use steric_x::{Molecule,PyramidalizationCalculator,bonded_neighbors};
fn main() {
    let paths: Vec<_>=std::env::args().skip(1).collect();
    let rows: Vec<_>=paths.par_iter().map(|path| {
        let m=Molecule::from_xyz_file(std::path::Path::new(path)).unwrap();
        let donors: Vec<_>=m.atoms.iter().enumerate().filter(|(_,a)|a.element=="P").map(|(i,_)|i).collect();
        assert_eq!(donors.len(),1);
        let d=donors[0];
        let ns: Vec<_>=bonded_neighbors(&m,d).iter().map(|(_,i)|*i).collect();
        let p=PyramidalizationCalculator::compute(&m,d,ns.try_into().unwrap()).unwrap();
        json!({"file":path,"pyr_p":p.pyr_p,"pyr_alpha":p.pyr_alpha})
    }).collect();
    serde_json::to_writer_pretty(std::io::stdout().lock(),&rows).unwrap();
    println!();
}
