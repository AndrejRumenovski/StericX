// Untimed public-API observation only. No descriptor implementation is copied.
use serde_json::json;
use steric_x::{Molecule, BuriedVolumeCalculator, BuriedVolumeConfig, SterimolCalculator,
    PyramidalizationCalculator, bonded_neighbors};
use steric_x::geometry::coordination_center_with_neighbors;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    for path in std::env::args().skip(1) {
        let m = Molecule::from_xyz_file(std::path::Path::new(&path))?;
        let donors: Vec<_> = m.atoms.iter().enumerate().filter(|(_,a)| a.element == "P").map(|(i,_)|i).collect();
        assert_eq!(donors.len(),1);
        let donor = donors[0];
        let neighbors: Vec<_> = bonded_neighbors(&m,donor).iter().map(|(_,i)|*i).collect();
        let neighbors: [usize;3] = neighbors.try_into().unwrap();
        let config = BuriedVolumeConfig { center_distance:2.28, ..Default::default() };
        let center = coordination_center_with_neighbors(&m,donor,neighbors,config)?;
        let b = BuriedVolumeCalculator::compute_with_neighbors(&m,donor,neighbors,config)?;
        let s = SterimolCalculator::compute_with_dummy(&m,donor,center)?;
        let p = PyramidalizationCalculator::compute(&m,donor,neighbors)?;
        println!("{}",json!({"file":path,"donor":donor,"neighbors":neighbors,
            "center":center.to_array(), "elements":m.atoms.iter().map(|a|a.element.clone()).collect::<Vec<_>>(),
            "coordinates":m.atoms.iter().map(|a|a.position.to_array()).collect::<Vec<_>>(),
            "radii":m.atoms.iter().map(|a|a.vdw_radius).collect::<Vec<_>>(),
            "sterimol_l":s.l+0.4,"sterimol_b1":s.b1,"sterimol_b5":s.b5,
            "pyr_p":p.pyr_p,"pyr_alpha":p.pyr_alpha,
            "buried_volume":b.buried_volume,"percent_buried_volume":b.percent_buried_volume,
            "qvbur_min":b.qvbur_min,"qvbur_max":b.qvbur_max,"max_delta_qvbur":b.max_delta_qvbur,
            "ovbur_min":b.ovbur_min,"ovbur_max":b.ovbur_max,"near_vbur":b.near_vbur,"far_vbur":b.far_vbur}));
    }
    Ok(())
}
