//! `stericx simulate`: Eyring rate and enantiomeric product distribution.

use crate::output::{print_memory_metrics, resident_memory_bytes};
use std::error::Error;
use std::mem::size_of;
use std::time::Instant;
use steric_x::EyringKineticLink;

pub(crate) fn simulate_command(ddg_kcal: f32, temp_k: f32) -> Result<(), Box<dyn Error>> {
    if !ddg_kcal.is_finite() {
        return Err("--ddg must be finite".into());
    }
    if !temp_k.is_finite() || temp_k <= 0.0 {
        return Err("--temp must be a positive finite temperature".into());
    }

    let total_started = Instant::now();
    let rss_start = resident_memory_bytes();
    let rate = EyringKineticLink::calculate_rate_constant(ddg_kcal, temp_k);
    let (major_percent, minor_percent) =
        EyringKineticLink::calculate_enantiomeric_ratio(ddg_kcal, temp_k);
    let ee_percent = EyringKineticLink::calculate_enantiomeric_excess(ddg_kcal, temp_k);
    let distribution = EyringKineticLink::product_ratio(ddg_kcal, temp_k);
    let total_time = total_started.elapsed();

    println!("command=simulate");
    println!("ddg_kcal_mol={ddg_kcal:.7}");
    println!("temperature_k={temp_k:.2}");
    println!("rate_constant_s^-1={rate:.7e}");
    println!("major_enantiomer_percent={major_percent:.4}");
    println!("minor_enantiomer_percent={minor_percent:.4}");
    println!("percent_r={:.4}", distribution.percent_r);
    println!("percent_s={:.4}", distribution.percent_s);
    println!("ee_percent={ee_percent:.4}");
    println!("calculations_performed=3");
    println!("input_bytes={}", 2 * size_of::<f32>());
    println!("output_bytes={}", 6 * size_of::<f32>());
    println!("total_microseconds={:.3}", total_time.as_secs_f64() * 1e6);
    print_memory_metrics(rss_start, resident_memory_bytes());
    Ok(())
}
