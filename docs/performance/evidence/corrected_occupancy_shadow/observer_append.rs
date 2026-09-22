

// NONPRODUCTION DIAGNOSTIC ONLY. Everything above this marker is an EXACT
// byte-identical copy of the captured CURRENT production BV source.
// The original functions above still produce every CLI scientific output.
// The instrumented copy below always executes each original predicate.

#[derive(Default, serde::Serialize)]
struct ProbeCounts {
    rows: u64, point_visits: u64, xy_preparations: u64, xy_rejections: u64,
    cached_z_tests: u64, new_candidate_z_tests: u64,
    first_hits_cached: u64, first_hits_new: u64, point_misses: u64,
    endpoint_checks_lazy: u64, endpoint_rejections_lazy: u64,
    inside_z_range_lazy: u64, nonfinite_fallback_lazy: u64,
    eliminated_cached_z_tests: u64, eliminated_new_z_tests: u64,
    asserted_eliminated_tests_are_misses: u64,
    eager_xy_preparations: u64, eager_xy_rejections: u64,
    eager_endpoint_checks: u64, eager_endpoint_rejections: u64,
    eager_inside_z_range: u64, eager_nonfinite_fallback: u64,
}

struct ProbeRow { start: usize, end: usize, finite: bool, minimum: f32, maximum: f32 }

fn probe_rows(points: &[Vec3]) -> Vec<ProbeRow> {
    let mut rows = Vec::new(); let mut start = 0;
    while start < points.len() {
        let mut end = start + 1;
        while end < points.len() && points[end].x.to_bits() == points[start].x.to_bits()
            && points[end].y.to_bits() == points[start].y.to_bits() { end += 1; }
        let finite = points[start..end].iter().all(|p| p.is_finite());
        let minimum = points[start..end].iter().map(|p| p.z).fold(f32::INFINITY, f32::min);
        let maximum = points[start..end].iter().map(|p| p.z).fold(f32::NEG_INFINITY, f32::max);
        rows.push(ProbeRow { start, end, finite, minimum, maximum }); start = end;
    }
    rows
}

// 0: conservative fallback, 1: atom within the actual row z range,
// 2: closest actual endpoint hits, 3: closest actual endpoint misses.
fn probe_endpoint(row: &ProbeRow, atom: &AlignedAtom, xy: f32) -> u8 {
    if !row.finite || !atom.position.is_finite() || !atom.radius_squared.is_finite()
        || !xy.is_finite() { return 0; }
    let endpoint = if atom.position.z < row.minimum { row.minimum }
        else if atom.position.z > row.maximum { row.maximum } else { return 1; };
    let dz = endpoint - atom.position.z;
    if xy + dz * dz <= atom.radius_squared { 2 } else { 3 }
}

fn probe_bits(v: &OccupiedVolumes) -> Vec<u32> {
    let mut bits = vec![v.buried_volume.to_bits(), v.near_vbur.to_bits(), v.far_vbur.to_bits()];
    bits.extend(v.quadrants.map(f32::to_bits)); bits.extend(v.octants.map(f32::to_bits)); bits
}

type ProbeWriter = std::sync::Mutex<Option<std::io::BufWriter<std::fs::File>>>;
static PROBE_WRITER: std::sync::OnceLock<ProbeWriter> = std::sync::OnceLock::new();

fn probe_record(value: serde_json::Value) {
    use std::io::Write;
    let output = PROBE_WRITER.get_or_init(|| std::sync::Mutex::new(
        std::env::var_os("STERICX_ROW_RANGE_PROBE_PATH").map(|path| {
            std::io::BufWriter::new(std::fs::OpenOptions::new().write(true).create_new(true)
                .open(path).expect("fresh diagnostic sidecar"))
        })));
    if let Some(writer) = output.lock().unwrap().as_mut() {
        serde_json::to_writer(&mut *writer, &value).unwrap(); writer.write_all(b"\n").unwrap();
    }
}

pub fn row_range_probe_finish() {
    use std::io::Write;
    if let Some(writer) = PROBE_WRITER.get()
        && let Some(writer) = writer.lock().unwrap().as_mut() { writer.flush().unwrap(); }
}

pub fn row_range_probe_conformer(molecule: &Molecule, donor: usize, neighbors: [usize;3],
    config: BuriedVolumeConfig, source: &str, conformer: usize) {
    let center = coordination_center_with_neighbors(molecule, donor, neighbors, config).unwrap();
    let points = integration_grid(config);
    for (orientation, plane) in neighbors.into_iter().enumerate() {
        let basis = coordinate_basis(molecule, donor, plane, center).unwrap();
        let atoms = aligned_atoms(molecule, center, basis, config).unwrap();
        let context = serde_json::json!({"source":source,"conformer":conformer,
            "orientation":orientation,"plane":plane,"atoms":atoms.len(),"points":points.len()});
        let _ = probe_instrumented_occupied_volumes(&points, &atoms, config.sphere_radius, &context);
    }
}

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


#[cfg(test)]
mod row_range_shadow_tests {
    use super::*;
    #[test]
    fn row_range_shadow_preserves_boundaries_nonfinite_and_arbitrary_order() {
        let special = [0.0, -0.0, f32::from_bits(1), -f32::from_bits(1),
            f32::MIN_POSITIVE, -1., 1., 1e20, f32::MAX, f32::INFINITY,
            f32::NEG_INFINITY, f32::NAN];
        let points: Vec<_> = special.into_iter().map(|z| Vec3::new(0.,0.,z)).collect();
        for z in special { for r2 in special {
            let atoms = [AlignedAtom{position:Vec3::new(0.,0.,z),radius_squared:r2}];
            probe_instrumented_occupied_volumes(&points,&atoms,3.5,&serde_json::json!({"test":"nonfinite"}));
        }}
        let mut points = vec![Vec3::new(0.5,0.75,-0.25),Vec3::new(0.5,0.75,0.),
            Vec3::new(0.5,0.75,0.25),Vec3::new(-0.,0.,0.),Vec3::new(0.,0.,0.),
            Vec3::new(0.5,0.75,4.),Vec3::new(0.5,0.75,-4.)];
        for z in [-4.,-0.25,0.,0.25,4.] {
            for point in points.clone() {
                let p = Vec3::new(0.,0.,z); let r2=point.distance_squared(p);
                for radius_squared in [f32::from_bits(r2.to_bits().saturating_sub(1)),r2,
                    f32::from_bits(r2.to_bits()+1)] {
                    let atoms=[AlignedAtom{position:p,radius_squared}];
                    probe_instrumented_occupied_volumes(&points,&atoms,3.5,&serde_json::json!({"test":"boundary"}));
                    points.reverse();
                }
            }
        }
        let atoms=[AlignedAtom{position:Vec3::ZERO,radius_squared:1.0}];
        probe_instrumented_occupied_volumes(&[],&atoms,3.5,&serde_json::json!({"test":"empty"}));
        probe_instrumented_occupied_volumes(&points,&[],3.5,&serde_json::json!({"test":"no atoms"}));
        let finite_points=[Vec3::new(0.,0.,f32::MAX),Vec3::new(0.,0.,-f32::MAX),
            Vec3::new(0.,0.,1e20),Vec3::new(0.,0.,-1e20),Vec3::ZERO];
        let mut atoms=[AlignedAtom{position:Vec3::new(0.,0.,1.),radius_squared:1.},
            AlignedAtom{position:Vec3::new(0.,0.,-f32::MAX),radius_squared:f32::MAX},
            AlignedAtom{position:Vec3::ZERO,radius_squared:f32::NAN}];
        for _ in 0..atoms.len() {
            probe_instrumented_occupied_volumes(&finite_points,&atoms,3.5,&serde_json::json!({"test":"finite overflow and atom order"}));
            probe_instrumented_occupied_volumes(&points,&atoms,3.5,&serde_json::json!({"test":"atom order"}));
            atoms.rotate_left(1);
        }
    }
}

#[cfg(test)]
mod corrected_row_range_edge_tests {
    use super::*;
    #[test]
    fn finite_absorption_underflow_overflow_and_boundary_are_exercised() {
        let cases = [
            (vec![0.0, 0.00001], -0.0001, 1.0, 1.0, 2),
            (vec![-0.00001, 0.0], 0.0001, 1.0, 1.0, 2),
            (vec![0.0, f32::from_bits(1)], -f32::from_bits(1), 0.0, 0.0, 2),
            (vec![f32::MAX / 2.0, f32::MAX], -f32::MAX, 0.0, f32::MAX, 3),
            (vec![-f32::MAX, -f32::MAX / 2.0], f32::MAX, 0.0, f32::MAX, 3),
            (vec![1.0, 2.0], 0.0, 1.0, f32::from_bits(2.0_f32.to_bits()-1), 3),
            (vec![1.0, 2.0], 0.0, 1.0, 2.0, 2),
            (vec![1.0, 2.0], 0.0, 1.0, f32::from_bits(2.0_f32.to_bits()+1), 2),
            (vec![-2.0, 2.0, 0.0], 0.0, 0.0, 1.0, 1),
        ];
        for (index,(zs, az, x, r2, expected)) in cases.into_iter().enumerate() {
            let mut points: Vec<_> = zs.into_iter().map(|z| Vec3::new(x,0.0,z)).collect();
            let atom=AlignedAtom{position:Vec3::new(0.0,0.0,az),radius_squared:r2};
            for reverse in [false,true] {
                if reverse {points.reverse();}
                let rows=probe_rows(&points);
                let dx=points[0].x-atom.position.x;
                let dy=points[0].y-atom.position.y;
                assert_eq!(probe_endpoint(&rows[0],&atom,dx*dx+dy*dy),expected);
                probe_instrumented_occupied_volumes(&points,std::slice::from_ref(&atom),3.5,
                    &serde_json::json!({"test":"finite focused endpoint","case":index,"reverse":reverse}));
            }
        }
        let points=[Vec3::new(0.0,0.0,2.0),Vec3::new(0.0,0.0,f32::NAN)];
        let atom=AlignedAtom{position:Vec3::ZERO,radius_squared:1.0};
        assert_eq!(probe_endpoint(&probe_rows(&points)[0],&atom,0.0),0);
        let row=probe_rows(&[Vec3::new(0.0,0.0,2.0)]);
        for (position,radius_squared,xy) in [
            (Vec3::ZERO,f32::INFINITY,0.0),
            (Vec3::new(0.0,0.0,f32::INFINITY),1.0,0.0),
            (Vec3::ZERO,1.0,f32::INFINITY),
            (Vec3::ZERO,f32::NAN,0.0),
        ] { assert_eq!(probe_endpoint(&row[0],&AlignedAtom{position,radius_squared},xy),0); }
    }
}
