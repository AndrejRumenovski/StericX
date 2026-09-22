// Untimed visibility probe appended to a byte-identical copy of current source.
pub fn grid_probe(m: &Molecule, d: usize, ns: [usize;3], cfg: BuriedVolumeConfig) -> serde_json::Value {
    let center = coordination_center_with_neighbors(m,d,ns,cfg).unwrap();
    let grid = integration_grid(cfg);
    let mut frames = Vec::new();
    for plane in ns {
        let basis = coordinate_basis(m,d,plane,center).unwrap();
        let atoms = aligned_atoms(m,center,basis,cfg).unwrap();
        let v = occupied_volumes(&grid,&atoms,cfg.sphere_radius);
        let mut populations = [0usize;8];
        let mut occupied = [0usize;8];
        for p in &grid {
            let index = octant_index(*p);
            populations[index] += 1;
            if atoms.iter().any(|a| p.distance_squared(a.position) <= a.radius_squared) {
                occupied[index] += 1;
            }
        }
        frames.push(serde_json::json!({"plane":plane,"grid_populations":populations,"occupied_counts":occupied,
            "buried_volume":v.buried_volume,"quadrants":v.quadrants,"octants":v.octants,"near_vbur":v.near_vbur,"far_vbur":v.far_vbur}));
    }
    serde_json::json!({"grid_count":grid.len(),"sphere_volume_f32":sphere_volume(cfg.sphere_radius),"frames":frames})
}
