
pub fn policy_probe(m: &Molecule,d:usize,ns:[usize;3],cfg:BuriedVolumeConfig)->serde_json::Value {
    let center=coordination_center_with_neighbors(m,d,ns,cfg).unwrap();
    let grid=integration_grid(cfg);
    let points:Vec<_>=grid.iter().map(|p| {
        let index=p.to_array().map(|x| ((f64::from(x)+3.5)*31.0/7.0).round() as usize);
        serde_json::json!({"index":index,"coordinates":p.to_array(),"region":octant_index(*p)})
    }).collect();
    let frames:Vec<_>=ns.iter().map(|&plane| {
        let basis=coordinate_basis(m,d,plane,center).unwrap();
        let atoms=aligned_atoms(m,center,basis,cfg).unwrap();
        let occupied:Vec<_>=grid.iter().map(|p| atoms.iter().any(|a| p.distance_squared(a.position)<=a.radius_squared)).collect();
        let v=occupied_volumes(&grid,&atoms,cfg.sphere_radius);
        serde_json::json!({"plane":plane,"occupied":occupied,"buried_volume":v.buried_volume,"quadrants":v.quadrants,"octants":v.octants,"near_vbur":v.near_vbur,"far_vbur":v.far_vbur,
            "aligned_atoms":atoms.iter().map(|a|serde_json::json!({"position":a.position.to_array(),"radius_squared":a.radius_squared})).collect::<Vec<_>>()})
    }).collect();
    serde_json::json!({"donor":d,"neighbors":ns,"center":center.to_array(),"points":points,"frames":frames,"sphere_volume":sphere_volume(cfg.sphere_radius)})
}
