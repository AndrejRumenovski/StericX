//! Identical observer timing boundary for public precision methods.
use serde_json::{Value,json};
use steric_x::geometry::{AtomicSphere,CertifiedVolumeOptions,CertifiedVolumeResult,ContinuousVolumeMethod,ContinuousVolumeResult,sphere_union_volume,sphere_union_volume_certified};
#[cfg(feature="profiling")]
#[global_allocator]
static ALLOCATOR:steric_x::profiling::TrackingAllocator=steric_x::profiling::TrackingAllocator;
enum Output { Certified(CertifiedVolumeResult), Legacy(ContinuousVolumeResult) }
fn main()->Result<(),Box<dyn std::error::Error>> {
    let qs:Vec<Value>=serde_json::from_str(&std::fs::read_to_string(std::env::args().nth(1).ok_or("request path")?)?)?;
    for q in qs {
        let region:AtomicSphere=serde_json::from_value(q["region"].clone())?;let atoms:Vec<AtomicSphere>=serde_json::from_value(q["atoms"].clone())?;
        let method=q["kind"].as_str().unwrap();
        let run=||->Result<Output,steric_x::geometry::DescriptorError>{
            if method=="certified" {Ok(Output::Certified(sphere_union_volume_certified(region,&atoms,CertifiedVolumeOptions{absolute_tolerance:q["tolerance"].as_f64().unwrap(),max_cells:q["budget"].as_u64().unwrap() as usize})?))}
            else {let m=if method=="reference" {ContinuousVolumeMethod::Reference{absolute_tolerance:q["tolerance"].as_f64().unwrap(),max_boxes:q["budget"].as_u64().unwrap() as usize}}else{ContinuousVolumeMethod::MidpointSlices{subdivisions:q["subdivisions"].as_u64().unwrap() as usize}};Ok(Output::Legacy(sphere_union_volume(region,&atoms,m)?))}
        };
        let warmups=q["warmups"].as_u64().unwrap_or(1);let samples=q["samples"].as_u64().unwrap_or(1);let iterations=q["iterations"].as_u64().unwrap_or(1);assert!(iterations>0);
        for _ in 0..warmups {std::hint::black_box(run()?);}
        let mut times=Vec::new();let mut allocations=Vec::new();let mut result=None;
        for _ in 0..samples {
            #[cfg(feature="profiling")]let before=steric_x::profiling::allocation_snapshot();
            let start=std::time::Instant::now();for _ in 0..iterations {result=Some(std::hint::black_box(run()?));}times.push(start.elapsed().as_nanos());
            #[cfg(feature="profiling")]let allocation={let a=steric_x::profiling::allocation_snapshot();json!({"calls":a.allocation_calls-before.allocation_calls,"reallocations":a.reallocation_calls-before.reallocation_calls,"requested_bytes":a.requested_bytes-before.requested_bytes,"process_peak_live_bytes":a.peak_live_requested_bytes})};
            #[cfg(not(feature="profiling"))]let allocation=Value::Null;
            allocations.push(allocation);
        }
        let result=match match result {Some(r)=>r,None=>run()?} {Output::Certified(r)=>json!({"kind":method,"result":r}),Output::Legacy(r)=>json!({"kind":method,"result":r})};
        println!("{}",json!({"id":q["id"],"iterations":iterations,"elapsed_ns":times,"allocations":allocations,"output":result}));
    }Ok(())
}
