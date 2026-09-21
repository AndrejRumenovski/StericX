//! Measurement adapter for frozen StericX. NOT an independent reference.
use glam::Vec3;
use serde_json::{json, Value, Map};
use std::io::{self, BufRead, Write};
pub use steric_x::{Molecule, covalent_radius};
use steric_x::{Atom, BuriedVolumeCalculator as BV, BuriedVolumeConfig as Config,
    BuriedVolumeParams as Params, SterimolCalculator, SterimolParams,
    PyramidalizationCalculator, EyringKineticLink as Eyring};

#[macro_export]
macro_rules! profile_scope { ($($arg:tt)*) => {}; }
#[path = "../../geometry/frozen_buried_volume.rs"]
#[allow(dead_code)]
mod observed_bv;
mod model_observation;

fn num(v: &Value) -> f32 {
    if let Some(x)=v.as_f64() { x as f32 } else {
        match v.as_str().unwrap_or("") { "NaN"=>f32::NAN,"Infinity"|"inf"=>f32::INFINITY,"-Infinity"|"-inf"=>f32::NEG_INFINITY, s=>s.parse().unwrap_or(f32::NAN) }
    }
}
fn field(v:&Value,key:&str,default:f32)->f32 { v.get(key).map(num).unwrap_or(default) }
fn idx(v:&Value,key:&str,default:usize)->usize { v[key].as_u64().map(|x|x as usize).unwrap_or(default) }
fn vec3(v:&Value)->Vec3 { Vec3::new(num(&v[0]),num(&v[1]),num(&v[2])) }
fn floats(pairs:&[(&str,f32)])->Value {
    let mut out=Map::new(); let mut bits=Map::new(); let mut nonfinite=Map::new();
    for &(key,x) in pairs {out.insert(key.into(),json!(x));bits.insert(key.into(),json!(format!("{:08x}",x.to_bits())));if !x.is_finite(){nonfinite.insert(key.into(),json!(x.to_string()));}}
    out.insert("_bits".into(),json!(bits));out.insert("_nonfinite".into(),json!(nonfinite));json!(out)
}
fn sterimol(p:SterimolParams)->Value {floats(&[("l",p.l),("b1",p.b1),("b5",p.b5)])}
fn bv(p:Params)->Value {floats(&[("buried_volume",p.buried_volume),("percent_buried_volume",p.percent_buried_volume),("qvbur_min",p.qvbur_min),("qvbur_max",p.qvbur_max),("max_delta_qvbur",p.max_delta_qvbur),("ovbur_min",p.ovbur_min),("ovbur_max",p.ovbur_max),("near_vbur",p.near_vbur),("far_vbur",p.far_vbur)])}
fn config(v:&Value)->Config {let d=Config::default();Config{sphere_radius:field(v,"sphere_radius",d.sphere_radius),density:field(v,"density",d.density),center_distance:field(v,"center_distance",d.center_distance),radii_scale:field(v,"radii_scale",d.radii_scale),include_hydrogens:v["include_hydrogens"].as_bool().unwrap_or(d.include_hydrogens)}}
fn copied_config(c:Config)->observed_bv::BuriedVolumeConfig {observed_bv::BuriedVolumeConfig{sphere_radius:c.sphere_radius,density:c.density,center_distance:c.center_distance,radii_scale:c.radii_scale,include_hydrogens:c.include_hydrogens}}
fn geometry(v:&Value)->Value {
    let atoms=v["atoms"].as_array().expect("atoms array").iter().map(|a|{let mut atom=Atom::new(a["element"].as_str().expect("element string"),vec3(&a["position"]));if let Some(r)=a.get("radius"){atom.vdw_radius=num(r)} atom}).collect();
    let m=Molecule{atoms}; let d=idx(v,"donor",0); let r=idx(v,"reference",1);let c=config(&v["config"]);let mut out=json!({});
    out["atoms"]=json!(m.atoms.iter().map(|a|json!({"element":a.element,"position":a.position.to_array(),"vdw_radius":a.vdw_radius,"covalent_radius":covalent_radius(&a.element)})).collect::<Vec<_>>());
    out["sterimol_bond"]=sterimol(SterimolCalculator::compute(&m,idx(v,"attach",d),idx(v,"sterimol_neighbor",r)));
    let neighbors=if d<m.atoms.len(){steric_x::bonded_neighbors(&m,d)}else{vec![]};out["bonded_neighbors"]=json!(neighbors.iter().map(|(distance_squared,index)|json!({"index":index,"distance_squared":distance_squared})).collect::<Vec<_>>());
    let n=if let Some(a)=v["neighbors"].as_array(){if a.len()==3{Some([a[0].as_u64().unwrap() as usize,a[1].as_u64().unwrap() as usize,a[2].as_u64().unwrap() as usize])}else{None}}else if neighbors.len()==3 {Some([neighbors[0].1,neighbors[1].1,neighbors[2].1])}else{None};
    out["pyramidalization"]=match n {Some(n)=>{let p=PyramidalizationCalculator::compute(&m,d,n);floats(&[("pyr_p",p.pyr_p),("pyr_alpha",p.pyr_alpha)])},None=>json!({"error":"observer requires exactly three explicit or inferred neighbors"})};
    let center=if let Some(center)=v.get("center"){Ok(vec3(center))}else{steric_x::coordination_center(&m,d,r,c).map_err(|e|e.to_string())};
    match center {Ok(center)=>{out["center"]=json!(center.to_array());let p=SterimolCalculator::compute_with_dummy(&m,d,center);out["sterimol_dummy_raw"]=sterimol(p);out["sterimol_coordination"]=sterimol(SterimolParams{l:p.l+0.40,..p});},Err(ref e)=>{out["center"]=json!({"error":e});out["sterimol_dummy_raw"]=json!({"error":e});out["sterimol_coordination"]=json!({"error":e});}}
    if !v["sterimol_only"].as_bool().unwrap_or(false) {
        let result=if let Some(x)=v.get("center"){BV::compute_with_center(&m,d,r,vec3(x),c)}else{BV::compute(&m,d,r,c)};
        out["buried_volume"]=match result{Ok(p)=>bv(p),Err(e)=>json!({"error":e.to_string()})};
        if v["dump"].as_bool().unwrap_or(false){out["dump"]=observed_bv::audit_dump(&m,d,r,v.get("center").map(vec3),copied_config(c),v["include_points"].as_bool().unwrap_or(false));}
    }
    out
}
fn kinetics(v:&Value)->Value {let d=num(&v["ddg"]);let t=num(&v["temperature"]);let (major,minor)=Eyring::calculate_enantiomeric_ratio(d,t);let p=Eyring::product_ratio(d,t);floats(&[("input_ddg",d),("input_temperature",t),("rate",Eyring::calculate_rate_constant(d,t)),("major_percent",major),("minor_percent",minor),("r_percent",p.percent_r),("s_percent",p.percent_s),("ee_percent",p.ee_percent)])}
fn aggregate(v:&Value)->Value {
    let cs=v["conformers"].as_array().expect("conformers array").iter().map(|p|Params{buried_volume:field(p,"buried_volume",0.0),percent_buried_volume:field(p,"percent_buried_volume",0.0),qvbur_min:field(p,"qvbur_min",0.0),qvbur_max:field(p,"qvbur_max",0.0),max_delta_qvbur:field(p,"max_delta_qvbur",0.0),ovbur_min:field(p,"ovbur_min",0.0),ovbur_max:field(p,"ovbur_max",0.0),near_vbur:field(p,"near_vbur",0.0),far_vbur:field(p,"far_vbur",0.0)}).collect::<Vec<_>>();
    let weights=v["weights"].as_array().expect("weights array").iter().map(num).collect::<Vec<_>>();
    match BV::aggregate(&cs,&weights){Err(e)=>json!({"error":e.to_string()}),Ok(e)=>{let mut o=floats(&[("vbur_boltz",e.vbur_boltz),("vbur_min",e.vbur_min),("vbur_max",e.vbur_max),("vbur_delta",e.vbur_delta),("qvbur_min_boltz",e.qvbur_min_boltz),("qvbur_max_boltz",e.qvbur_max_boltz),("max_delta_qvbur_boltz",e.max_delta_qvbur_boltz),("max_delta_qvbur_min",e.max_delta_qvbur_min),("max_delta_qvbur_max",e.max_delta_qvbur_max),("max_delta_qvbur_delta",e.max_delta_qvbur_delta),("max_delta_qvbur_vburminconf",e.max_delta_qvbur_vburminconf),("near_vbur_boltz",e.near_vbur_boltz),("far_vbur_boltz",e.far_vbur_boltz)]);o["conformer_count"]=json!(e.conformer_count);o}}
}
fn observe(v:&Value)->Value {match v["op"].as_str().unwrap_or("geometry") {"geometry"=>geometry(v),"kinetics"=>kinetics(v),"aggregate"=>aggregate(v),"model"=>model_observation::observe(v),"raw_occupancy"=>{let points=v["points"].as_array().unwrap().iter().map(vec3).collect::<Vec<_>>();let atoms=v["atoms"].as_array().unwrap().iter().map(|a|(vec3(&a["position"]),num(&a["radius_squared"]))).collect::<Vec<_>>();observed_bv::audit_raw(&points,&atoms,field(v,"sphere_radius",3.5))},_=>json!({"error":"unknown operation"})}}
fn main() {
    let stdin=io::stdin();let mut stdout=io::BufWriter::new(io::stdout().lock());
    for line in stdin.lock().lines(){let line=line.expect("read stdin");if line.trim().is_empty(){continue}let result=std::panic::catch_unwind(||{let v:Value=serde_json::from_str(&line).expect("JSON request");let mut out=observe(&v);out["id"]=v["id"].clone();out});let out=result.unwrap_or_else(|_|json!({"error":"observer caught panic; see stderr"}));writeln!(stdout,"{}",out).expect("write stdout");stdout.flush().expect("flush");}
}
