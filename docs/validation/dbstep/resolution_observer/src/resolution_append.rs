// Observation only: all occupancy decisions call the unchanged production kernel.
// This adapter is the system under investigation, not an independent reference.

pub fn audit_grid_occupancy(
    atoms: &[(Vec3, f32)],
    config: BuriedVolumeConfig,
    include_points: bool,
) -> serde_json::Value {
    if let Err(error) = validate_config(config) {
        return serde_json::json!({"error": error.to_string()});
    }
    let points = integration_grid(config);
    audit_point_decisions(
        &points,
        atoms,
        config,
        include_points,
        "native_integration_grid",
    )
}

pub fn audit_raw_point_decisions(
    points: &[Vec3],
    atoms: &[(Vec3, f32)],
    config: BuriedVolumeConfig,
    include_points: bool,
) -> serde_json::Value {
    if let Err(error) = validate_config(config) {
        return serde_json::json!({"error": error.to_string()});
    }
    audit_point_decisions(points, atoms, config, include_points, "supplied_points")
}

fn audit_point_decisions(
    points: &[Vec3],
    input_atoms: &[(Vec3, f32)],
    config: BuriedVolumeConfig,
    include_points: bool,
    point_source: &str,
) -> serde_json::Value {
    let volume = sphere_volume(config.sphere_radius);
    if !volume.is_normal() || volume <= 0.0 {
        return serde_json::json!({"error": "point-decision observation requires a positive normal sphere volume"});
    }
    if points.is_empty() {
        return serde_json::json!({"error": "point-decision observation requires a nonempty point set"});
    }
    let atoms: Vec<_> = input_atoms
        .iter()
        .map(|(position, radius_squared)| AlignedAtom {
            position: *position,
            radius_squared: *radius_squared,
        })
        .collect();
    let aggregate = occupied_volumes(points, &atoms, config.sphere_radius);
    let mut mask = Vec::with_capacity(points.len());
    let mut bins = Vec::with_capacity(points.len());
    let mut quadrant_total = [0_usize; 4];
    let mut quadrant_occupied = [0_usize; 4];
    let mut octant_total = [0_usize; 8];
    let mut octant_occupied = [0_usize; 8];
    for point in points {
        let observed = occupied_volumes(std::slice::from_ref(point), &atoms, config.sphere_radius);
        // A one-point call returns either zero or the positive sphere volume.
        // No distance, radius comparison or occupancy predicate is copied here.
        if observed.buried_volume != 0.0 && observed.buried_volume.to_bits() != volume.to_bits() {
            return serde_json::json!({"error": "one-point occupancy did not produce zero or the sphere volume"});
        }
        let occupied = observed.buried_volume != 0.0;
        let quadrant = quadrant_index(*point);
        let octant = octant_index(*point);
        quadrant_total[quadrant] += 1;
        octant_total[octant] += 1;
        if occupied {
            quadrant_occupied[quadrant] += 1;
            octant_occupied[octant] += 1;
        }
        mask.push(occupied);
        bins.push([quadrant, octant]);
    }
    serde_json::json!({
        "point_source": point_source,
        "config": config,
        "atoms": atoms.iter().map(|atom| serde_json::json!({
            "position": atom.position.to_array(),
            "position_bits": atom.position.to_array().map(f32::to_bits),
            "radius_squared": atom.radius_squared,
            "radius_squared_bits": atom.radius_squared.to_bits(),
        })).collect::<Vec<_>>(),
        "grid_count": points.len(),
        "sphere_volume": volume,
        "sphere_volume_bits": volume.to_bits(),
        "aggregate": {
            "buried_volume": aggregate.buried_volume,
            "quadrants": aggregate.quadrants,
            "octants": aggregate.octants,
            "near_vbur": aggregate.near_vbur,
            "far_vbur": aggregate.far_vbur,
            "_bits": {
                "buried_volume": aggregate.buried_volume.to_bits(),
                "quadrants": aggregate.quadrants.map(f32::to_bits),
                "octants": aggregate.octants.map(f32::to_bits),
                "near_vbur": aggregate.near_vbur.to_bits(),
                "far_vbur": aggregate.far_vbur.to_bits(),
            },
        },
        "occupied_count": mask.iter().filter(|&&occupied| occupied).count(),
        "occupied_mask": mask,
        "point_bins": bins,
        "quadrant_total": quadrant_total,
        "quadrant_occupied": quadrant_occupied,
        "octant_total": octant_total,
        "octant_occupied": octant_occupied,
        "points": if include_points {
            Some(points.iter().map(|point| point.to_array()).collect::<Vec<_>>())
        } else {
            None
        },
        "point_bits": if include_points {
            Some(points.iter().map(|point| point.to_array().map(f32::to_bits)).collect::<Vec<_>>())
        } else {
            None
        },
    })
}
