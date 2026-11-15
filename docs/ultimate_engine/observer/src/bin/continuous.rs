//! Measures explicit continuous/finite-scan APIs; no reference equations.
use glam::Vec3;
use serde::{Deserialize, Serialize};
use serde_json::json;
use std::{fs, hint::black_box, time::Instant};
use steric_x::geometry::{Atom, Molecule, SterimolCalculator};

#[derive(Deserialize)]
struct InputAtom { position: [f32;3], radius: f32 }
#[derive(Deserialize)]
struct Case { id: String, atoms: Vec<InputAtom>, base: usize, dummy: [f32;3] }

fn measurement<T:Serialize>(method:&str, iterations:usize, run:impl Fn()->T) -> serde_json::Value {
    black_box(run());
    let mut samples = Vec::new();
    for _ in 0..5 {
        let start = Instant::now();
        for _ in 0..iterations { black_box(run()); }
        samples.push(start.elapsed().as_nanos());
    }
    json!({"method":method,"iterations":iterations,"elapsed_ns":samples,"result":run()})
}

fn main() -> Result<(),Box<dyn std::error::Error>> {
    let path = std::env::args().nth(1).ok_or("case JSON path required")?;
    let cases:Vec<Case> = serde_json::from_str(&fs::read_to_string(path)?)?;
    for case in cases {
        let molecule = Molecule { atoms:case.atoms.iter().map(|a| Atom {
            element:"C".into(),position:Vec3::from_array(a.position),vdw_radius:a.radius,
        }).collect() };
        let dummy = Vec3::from_array(case.dummy);
        let mut rows = vec![
            measurement("continuous_events_f64",20, || SterimolCalculator::compute_continuous_with_dummy(&molecule,case.base,dummy).unwrap()),
            measurement("default_360_f32",20, || SterimolCalculator::compute_with_dummy(&molecule,case.base,dummy).unwrap()),
        ];
        for n in [361,1441,5761,23041] {
            rows.push(measurement(&format!("endpoint_scan_{n}_f32"),3, ||
                SterimolCalculator::compute_with_dummy_with_sampling(&molecule,case.base,dummy,n).unwrap()));
        }
        println!("{}",json!({"id":case.id,"measurements":rows}));
    }
    Ok(())
}
