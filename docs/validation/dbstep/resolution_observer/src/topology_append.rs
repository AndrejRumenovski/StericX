// Additional topology lane only; original audit append stays unchanged.
pub fn audit_dump_with_neighbors(molecule: &Molecule, donor: usize, neighbors: [usize; 3], center: Option<Vec3>, config: BuriedVolumeConfig, include_points: bool) -> serde_json::Value {
    if let Err(e)=validate_config(config).and_then(|_|validate_neighbors(molecule,donor,neighbors)) {return serde_json::json!({"error":e.to_string()});}
    let center = match center { Some(c) => c, None => match coordination_center_with_neighbors(molecule, donor, neighbors, config) {Ok(c) => c, Err(e) => return serde_json::json!({"error":e.to_string()})}};
    let points = integration_grid(config);
    let mut orientations = Vec::new();
    for plane in neighbors {
        let basis = match coordinate_basis(molecule,donor,plane,center) {Ok(b) => b,Err(e)=>return serde_json::json!({"error":e.to_string(),"plane":plane})};
        let atoms=match aligned_atoms(molecule,center,basis,config) {Ok(a)=>a,Err(e)=>return serde_json::json!({"error":e.to_string()})};
        let v=occupied_volumes(&points,&atoms,config.sphere_radius);
        orientations.push(serde_json::json!({"plane":plane,"basis":[basis.x.to_array(),basis.y.to_array(),basis.z.to_array()],"aligned_atoms":atoms.iter().map(|a|serde_json::json!({"position":a.position.to_array(),"radius_squared":a.radius_squared})).collect::<Vec<_>>(),"buried_volume":v.buried_volume,"quadrants":v.quadrants,"octants":v.octants,"near_vbur":v.near_vbur,"far_vbur":v.far_vbur}));
    }
    serde_json::json!({"center":center.to_array(),"neighbors":neighbors,"grid_count":points.len(),"points":if include_points {Some(points.iter().map(|p|p.to_array()).collect::<Vec<_>>())} else {None},"orientations":orientations})
}
