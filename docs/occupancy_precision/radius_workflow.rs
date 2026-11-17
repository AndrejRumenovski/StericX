//! Public radius profile, with an untimed copied-predicate observation mode.
#![allow(dead_code)]
use rayon::prelude::*;
use serde_json::{Value,json};
use steric_x::{Molecule,BuriedVolumeConfig,BuriedVolumeCalculator,bonded_neighbors};
use steric_x::geometry::covalent_radius;
#[macro_export] macro_rules! profile_scope { ($($tokens:tt)*) => {}; }
mod copied {include!("../../policy_copied_volume.rs");}
fn calculate(path:&String)->Value {
    let m=Molecule::from_xyz_file(std::path::Path::new(path)).unwrap();
    let d=m.atoms.iter().position(|a|a.element=="P").unwrap();
    let ns:Vec<_>=bonded_neighbors(&m,d).iter().map(|(_,i)|*i).collect();let ns:[usize;3]=ns.try_into().unwrap();
    let rows:Vec<_>=[2.5_f64,3.5,4.5].into_iter().map(|r|{
        let cfg=BuriedVolumeConfig{sphere_radius:r as f32,density:(0.01*(r/3.5).powi(3)) as f32,center_distance:2.28,..Default::default()};
        json!({"radius":r,"values":BuriedVolumeCalculator::compute_with_neighbors(&m,d,ns,cfg).unwrap()})
    }).collect();json!({"file":path,"profile":rows})
}
fn main(){
    let mut args=std::env::args().skip(1);let mode=args.next().unwrap();
    if mode=="probe" {
        for path in args {
            let m=Molecule::from_xyz_file(std::path::Path::new(&path)).unwrap();let d=m.atoms.iter().position(|a|a.element=="P").unwrap();
            let ns:Vec<_>=bonded_neighbors(&m,d).iter().map(|(_,i)|*i).collect();let ns:[usize;3]=ns.try_into().unwrap();
            for r in [2.5_f64,3.5,4.5] {let cfg=copied::BuriedVolumeConfig{sphere_radius:r as f32,density:(0.01*(r/3.5).powi(3)) as f32,center_distance:2.28,..Default::default()};println!("{}",json!({"file":path,"radius":r,"probe":copied::policy_probe(&m,d,ns,cfg)}));}
        }return;
    }
    let iterations:usize=mode.parse().unwrap();let paths:Vec<_>=args.collect();let run=||paths.par_iter().map(calculate).collect::<Vec<_>>();let rows=run();let mut times=Vec::new();
    for _ in 0..iterations {let start=std::time::Instant::now();let r=std::hint::black_box(run());times.push(start.elapsed().as_nanos());assert_eq!(r,rows);}
    println!("{}",json!({"rows":rows,"warm_elapsed_ns":times}));
}
