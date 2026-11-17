// Observation adapter appended to a byte-identical production-source copy.
// Counter replay is an algorithmic model, never a hardware-counter measurement.
pub fn occupancy_probe(input: &serde_json::Value) -> serde_json::Value {
    use std::{hint::black_box, time::Instant};
    let config = BuriedVolumeConfig {
        sphere_radius: input["radius"].as_f64().unwrap_or(3.5) as f32,
        density: input["density"].as_f64().unwrap_or(0.01) as f32,
        ..BuriedVolumeConfig::default()
    };
    let atoms: Vec<_> = input["atoms"].as_array().unwrap().iter().map(|a| AlignedAtom {
        position: Vec3::from_array(std::array::from_fn(|i| a["position"][i].as_f64().unwrap() as f32)),
        radius_squared: a["radius_squared"].as_f64().unwrap() as f32,
    }).collect();
    let start = Instant::now();
    let grid = integration_grid(config);
    let grid_ns = start.elapsed().as_nanos();
    let mut rows = Vec::new();
    let mut mismatches = 0;
    for row in grid.chunk_by(|a,b| a.x.to_bits()==b.x.to_bits() && a.y.to_bits()==b.y.to_bits()) {
        let mut xy_tests = 0;
        let mut candidates = 0;
        let mut z_tests = 0;
        let mut simd_vectors = 0;
        let mut scalar_tail = 0;
        let mut redundant_hits = 0;
        let mut hits = 0;
        for block in row.chunks(64) {
            let z: Vec<_> = block.iter().map(|p| p.z).collect();
            let full = u64::MAX >> (64-block.len());
            let mut mask = 0u64;
            for atom in &atoms {
                xy_tests += 1;
                let dx = row[0].x-atom.position.x;
                let dy = row[0].y-atom.position.y;
                let xy = dx*dx+dy*dy;
                if xy > atom.radius_squared { continue; }
                candidates += 1;
                z_tests += block.len();
                simd_vectors += block.len()/8;
                scalar_tail += block.len()%8;
                let expected = z.iter().enumerate().fold(0u64, |mask,(i,&p)| {
                    let dz = p-atom.position.z;
                    mask | (u64::from(xy+dz*dz<=atom.radius_squared)<<i)
                });
                #[cfg(any(target_arch="x86", target_arch="x86_64"))]
                if std::is_x86_feature_detected!("avx2") {
                    let actual = unsafe { row_mask_avx2(&z,atom.position.z,xy,atom.radius_squared) };
                    mismatches += usize::from(actual != expected);
                }
                redundant_hits += (mask & expected).count_ones() as usize;
                mask |= expected;
                if mask==full {break;}
            }
            for (i,point) in block.iter().enumerate() {
                let brute = atoms.iter().any(|a| point.distance_squared(a.position)<=a.radius_squared);
                mismatches += usize::from(brute != (mask & (1u64<<i)!=0));
            }
            hits += mask.count_ones() as usize;
        }
        rows.push(serde_json::json!({"points":row.len(),"xy_tests":xy_tests,"candidates":candidates,
            "z_tests":z_tests,"simd_vectors":simd_vectors,"scalar_tail":scalar_tail,
            "redundant_hits":redundant_hits,"hits":hits}));
    }
    let iterations = input["iterations"].as_u64().unwrap_or(100) as usize;
    let mut timing = Vec::new();
    for rep in 0..5 {
        for scalar in if rep%2==0 {[false,true]} else {[true,false]} {
            let start = Instant::now();
            for _ in 0..iterations {
                black_box(if scalar {occupied_volumes_scalar(black_box(&grid),black_box(&atoms),config.sphere_radius)}
                    else {occupied_volumes(black_box(&grid),black_box(&atoms),config.sphere_radius)});
            }
            timing.push(serde_json::json!({"rep":rep,"warmup":rep==0,"scalar":scalar,
                "elapsed_ns":start.elapsed().as_nanos(),"iterations":iterations}));
        }
    }
    serde_json::json!({"id":input["id"],"atoms":atoms.len(),"points":grid.len(),
        "rows":rows,"mask_mismatches":mismatches,"grid_ns":grid_ns,
        "grid_capacity_bytes":grid.capacity()*std::mem::size_of::<Vec3>(),"timing":timing})
}
