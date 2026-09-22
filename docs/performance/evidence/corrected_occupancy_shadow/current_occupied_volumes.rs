fn occupied_volumes(sphere: &[Vec3], atoms: &[AlignedAtom], sphere_radius: f32) -> OccupiedVolumes {
    crate::profile_scope!("buried_volume", "buried_volume::occupied_volumes");
    struct RowCandidate {
        z: f32,
        radius_squared: f32,
        xy_squared: f32,
    }

    let mut occupied_total = 0_usize;
    let mut quadrant_total = [0_usize; 4];
    let mut quadrant_occupied = [0_usize; 4];
    let mut octant_total = [0_usize; 8];
    let mut octant_occupied = [0_usize; 8];
    let mut row_xy = None;
    let mut next_atom = 0;
    let mut candidates: Vec<RowCandidate> = Vec::with_capacity(atoms.len().min(64));
    for point in sphere {
        let xy = (point.x.to_bits(), point.y.to_bits());
        if row_xy != Some(xy) {
            row_xy = Some(xy);
            next_atom = 0;
            candidates.clear();
        }
        let quadrant = quadrant_index(*point);
        let octant = octant_index(*point);
        quadrant_total[quadrant] += 1;
        octant_total[octant] += 1;
        // glam::Vec3 uses (dx*dx + dy*dy) + dz*dz. A row reuses the
        // bit-identical first sum without changing that f32 operation order.
        let mut occupied = candidates.iter().any(|candidate| {
            let dz = point.z - candidate.z;
            candidate.xy_squared + dz * dz <= candidate.radius_squared
        });
        // The cache covers retained atoms from [0, next_atom), in original
        // order. Extend it only after all cached candidates miss this point;
        // a dense first-atom hit therefore never prepares the unused suffix.
        while !occupied && next_atom < atoms.len() {
            let atom = &atoms[next_atom];
            next_atom += 1;
            let dx = point.x - atom.position.x;
            let dy = point.y - atom.position.y;
            let xy_squared = dx * dx + dy * dy;
            // Adding a non-negative z square cannot turn this miss into a hit.
            // A NaN comparison is false and retains the candidate; inf <= inf
            // behavior is preserved by the final test.
            if xy_squared > atom.radius_squared {
                continue;
            }
            candidates.push(RowCandidate {
                z: atom.position.z,
                radius_squared: atom.radius_squared,
                xy_squared,
            });
            let dz = point.z - atom.position.z;
            occupied = xy_squared + dz * dz <= atom.radius_squared;
        }
        if occupied {
            occupied_total += 1;
            quadrant_occupied[quadrant] += 1;
            octant_occupied[octant] += 1;
        }
    }
    let volume = sphere_volume(sphere_radius);
    let quadrants = std::array::from_fn(|index| {
        occupied_fraction(quadrant_occupied[index], quadrant_total[index]) * volume / 4.0
    });
    let octants = std::array::from_fn(|index| {
        occupied_fraction(octant_occupied[index], octant_total[index]) * volume / 8.0
    });
    let near_vbur = octants[4..].iter().sum();
    let far_vbur = octants[..4].iter().sum();
    OccupiedVolumes {
        buried_volume: occupied_fraction(occupied_total, sphere.len()) * volume,
        quadrants,
        octants,
        near_vbur,
        far_vbur,
    }
}
