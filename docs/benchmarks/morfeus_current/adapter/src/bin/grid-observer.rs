#![allow(dead_code)]
use steric_x::Molecule;
use steric_x::geometry::covalent_radius;
#[macro_export]
macro_rules! profile_scope { ($($tokens:tt)*) => {}; }
mod copied { include!("../copied_buried_volume.rs"); }
fn main() -> Result<(), Box<dyn std::error::Error>> {
    for path in std::env::args().skip(1) {
        let m = Molecule::from_xyz_file(std::path::Path::new(&path))?;
        let d = m.atoms.iter().position(|a|a.element=="P").unwrap();
        let ns: Vec<_> = steric_x::bonded_neighbors(&m,d).iter().map(|(_,i)|*i).collect();
        let cfg = copied::BuriedVolumeConfig { center_distance:2.28, ..Default::default() };
        println!("{}",copied::grid_probe(&m,d,ns.try_into().unwrap(),cfg));
    }
    Ok(())
}
