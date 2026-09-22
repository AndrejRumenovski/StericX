fn probe_instrumented_occupied_volumes(sphere: &[Vec3], atoms: &[AlignedAtom], sphere_radius: f32, context: &serde_json::Value) -> OccupiedVolumes {
    crate::profile_scope!("buried_volume", "buried_volume::occupied_volumes");
    let rows = probe_rows(sphere);
    let mut counts = ProbeCounts::default();
    let mut row_cursor = 0;
    let mut row_lengths = std::collections::BTreeMap::<usize,u64>::new();
    let mut eliminated_by_atom = vec![0_u64; atoms.len()];
    let mut first_hits_by_atom = vec![0_u64; atoms.len()+1];
    let mut witnesses = Vec::new();
    // Eager model is computed separately; these atoms are NOT added to the
    // actual lazy cache. It is an upper-work model, not an accepted algorithm.
    for row in &rows {
        *row_lengths.entry(row.end-row.start).or_default() += 1;
        for atom in atoms {
            counts.eager_xy_preparations += 1;
            let dx=sphere[row.start].x-atom.position.x;
            let dy=sphere[row.start].y-atom.position.y;
            let xy=dx*dx+dy*dy;
            if xy > atom.radius_squared { counts.eager_xy_rejections += 1; continue; }
            match probe_endpoint(row,atom,xy) {
                0 => counts.eager_nonfinite_fallback += 1,
                1 => counts.eager_inside_z_range += 1,
                2 => counts.eager_endpoint_checks += 1,
                3 => {counts.eager_endpoint_checks += 1;counts.eager_endpoint_rejections += 1;},
                _ => unreachable!(),
            }
        }
    }
    struct RowCandidate {
        atom_index: usize,
        shadow_skip: bool,
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
    for (point_index, point) in sphere.iter().enumerate() {
        counts.point_visits += 1;
        let mut hit_atom = atoms.len();
        let xy = (point.x.to_bits(), point.y.to_bits());
        if row_xy != Some(xy) {
            if row_xy.is_some() { row_cursor += 1; }
            counts.rows += 1;
            assert_eq!(rows[row_cursor].start, point_index);
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
            counts.cached_z_tests += 1;
            let dz = point.z - candidate.z;
            let accepted = candidate.xy_squared + dz * dz <= candidate.radius_squared;
            if candidate.shadow_skip {
                counts.eliminated_cached_z_tests += 1;
                eliminated_by_atom[candidate.atom_index] += 1;
                assert!(!accepted,"row-range rejected an occupied original cached test");
                counts.asserted_eliminated_tests_are_misses += 1;
                if witnesses.len()<8 { witnesses.push(serde_json::json!({"row":row_cursor,
                    "point_index":point_index,"atom":candidate.atom_index,"kind":"cached",
                    "point_bits":point.to_array().map(f32::to_bits),"xy_bits":candidate.xy_squared.to_bits(),
                    "radius_squared_bits":candidate.radius_squared.to_bits()})); }
            }
            if accepted {hit_atom=candidate.atom_index;counts.first_hits_cached += 1;}
            accepted
        });
        // The cache covers retained atoms from [0, next_atom), in original
        // order. Extend it only after all cached candidates miss this point;
        // a dense first-atom hit therefore never prepares the unused suffix.
        while !occupied && next_atom < atoms.len() {
            counts.xy_preparations += 1;
            let atom = &atoms[next_atom];
            next_atom += 1;
            let dx = point.x - atom.position.x;
            let dy = point.y - atom.position.y;
            let xy_squared = dx * dx + dy * dy;
            // Adding a non-negative z square cannot turn this miss into a hit.
            // A NaN comparison is false and retains the candidate; inf <= inf
            // behavior is preserved by the final test.
            if xy_squared > atom.radius_squared {
                counts.xy_rejections += 1;
                continue;
            }
            let endpoint = probe_endpoint(&rows[row_cursor],atom,xy_squared);
            match endpoint {
                0 => counts.nonfinite_fallback_lazy += 1,
                1 => counts.inside_z_range_lazy += 1,
                2 => counts.endpoint_checks_lazy += 1,
                3 => {counts.endpoint_checks_lazy += 1;counts.endpoint_rejections_lazy += 1;},
                _ => unreachable!(),
            }
            let shadow_skip = endpoint==3;
            candidates.push(RowCandidate {
                atom_index: next_atom-1,
                shadow_skip,
                z: atom.position.z,
                radius_squared: atom.radius_squared,
                xy_squared,
            });
            counts.new_candidate_z_tests += 1;
            let dz = point.z - atom.position.z;
            occupied = xy_squared + dz * dz <= atom.radius_squared;
            if shadow_skip {
                counts.eliminated_new_z_tests += 1;
                eliminated_by_atom[next_atom-1] += 1;
                assert!(!occupied,"row-range rejected an occupied original new test");
                counts.asserted_eliminated_tests_are_misses += 1;
                if witnesses.len()<8 { witnesses.push(serde_json::json!({"row":row_cursor,
                    "point_index":point_index,"atom":next_atom-1,"kind":"new",
                    "point_bits":point.to_array().map(f32::to_bits),"xy_bits":xy_squared.to_bits(),
                    "radius_squared_bits":atom.radius_squared.to_bits()})); }
            }
            if occupied {hit_atom=next_atom-1;counts.first_hits_new += 1;}
        }
        first_hits_by_atom[hit_atom] += 1;
        if !occupied {counts.point_misses += 1;}
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
    let result = OccupiedVolumes {
        buried_volume: occupied_fraction(occupied_total, sphere.len()) * volume,
        quadrants,
        octants,
        near_vbur,
        far_vbur,
    };
    let unchanged = occupied_volumes(sphere,atoms,sphere_radius);
    assert_eq!(probe_bits(&result),probe_bits(&unchanged),"diagnostic copy changed scientific bits");
    assert_eq!(counts.point_visits,counts.first_hits_cached+counts.first_hits_new+counts.point_misses);
    assert_eq!(counts.eliminated_cached_z_tests+counts.eliminated_new_z_tests,
        counts.asserted_eliminated_tests_are_misses);
    probe_record(serde_json::json!({"context":context,"counts":counts,
        "row_length_histogram":row_lengths,"eliminated_original_tests_by_atom":eliminated_by_atom,
        "original_first_hits_by_atom_or_miss":first_hits_by_atom,
        "elimination_witnesses":witnesses,"scientific_bits":probe_bits(&result),
        "matches_unmodified_current_function":true}));
    result
}
