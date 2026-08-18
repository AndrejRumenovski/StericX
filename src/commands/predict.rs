//! `stericx predict`: parallel SIMD regression over a mapped `.sigpack` matrix.

use crate::output::{millis, print_memory_metrics, resident_memory_bytes};
use serde::Deserialize;
use std::error::Error;
use std::fs;
use std::mem::size_of;
use std::path::Path;
use std::time::Instant;
use steric_x::model::MODEL_FEATURE_COUNT;
use steric_x::{RegressXPredictor, SigPackReader};

#[derive(Debug, Deserialize)]
#[serde(untagged)]
enum WeightDocument {
    Array([f32; MODEL_FEATURE_COUNT]),
    Object { weights: [f32; MODEL_FEATURE_COUNT] },
}

pub(crate) fn load_weights_json(path: &Path) -> Result<[f32; MODEL_FEATURE_COUNT], Box<dyn Error>> {
    let contents = fs::read_to_string(path)?;
    let document: WeightDocument = serde_json::from_str(&contents).map_err(|error| {
        format!(
            "weight file {} must be an eight-number JSON array or an object with an eight-number `weights` array: {error}",
            path.display()
        )
    })?;
    let weights = match document {
        WeightDocument::Array(weights) | WeightDocument::Object { weights } => weights,
    };
    if weights.iter().any(|weight| !weight.is_finite()) {
        return Err(format!("weight file {} contains NaN or infinity", path.display()).into());
    }
    Ok(weights)
}

pub(crate) fn predict_command(data: &Path, weights_path: &Path) -> Result<(), Box<dyn Error>> {
    let total_started = Instant::now();
    let rss_start = resident_memory_bytes();
    if !data.is_file() {
        return Err(format!("sigpack file does not exist: {}", data.display()).into());
    }
    if !weights_path.is_file() {
        return Err(format!("weight file does not exist: {}", weights_path.display()).into());
    }

    let weights_started = Instant::now();
    let weights = load_weights_json(weights_path)?;
    let weights_time = weights_started.elapsed();
    let predictor = RegressXPredictor::new(weights);

    let mapping_started = Instant::now();
    let reader = SigPackReader::open(data)?;
    let mapping_time = mapping_started.elapsed();
    let records = reader.records();
    if records.is_empty() {
        return Err("sigpack matrix contains no records".into());
    }

    let inference_started = Instant::now();
    let predictions = predictor.predict_batch(records);
    let inference_time = inference_started.elapsed();
    let mse = predictions
        .iter()
        .zip(records)
        .map(|(prediction, record)| {
            let residual = f64::from(*prediction) - f64::from(record.exp_ddg);
            residual * residual
        })
        .sum::<f64>()
        / records.len() as f64;
    let total_time = total_started.elapsed();
    let inference_throughput =
        records.len() as f64 / inference_time.as_secs_f64().max(f64::EPSILON);

    println!("command=predict");
    println!("records_predicted={}", records.len());
    println!("rayon_threads={}", rayon::current_num_threads());
    println!("sigpack_input={}", data.display());
    println!("weights_input={}", weights_path.display());
    println!("mapped_bytes={}", fs::metadata(data)?.len());
    println!(
        "prediction_buffer_bytes={}",
        predictions.len() * size_of::<f32>()
    );
    println!("weights_load_ms={:.3}", millis(weights_time));
    println!("memory_map_ms={:.3}", millis(mapping_time));
    println!("prediction_latency_ms={:.3}", millis(inference_time));
    println!("total_ms={:.3}", millis(total_time));
    println!("throughput_records_per_second={inference_throughput:.1}");
    println!("mse_kcal2_per_mol2={mse:.8}");
    for (index, prediction) in predictions.iter().take(5).enumerate() {
        println!(
            "prediction_preview[{index}]={prediction:.7},experimental={:.7}",
            records[index].exp_ddg
        );
    }
    print_memory_metrics(rss_start, resident_memory_bytes());
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::test_support::temporary_directory;

    #[test]
    fn weight_loader_accepts_array_and_object_json() {
        let directory = temporary_directory("weights");
        fs::create_dir_all(&directory).unwrap();
        let array_path = directory.join("array.json");
        let object_path = directory.join("object.json");
        fs::write(&array_path, "[0,1,2,3,4,5,6,7]").unwrap();
        fs::write(&object_path, r#"{"weights":[7,6,5,4,3,2,1,0]}"#).unwrap();

        assert_eq!(load_weights_json(&array_path).unwrap()[7], 7.0);
        assert_eq!(load_weights_json(&object_path).unwrap()[0], 7.0);
        fs::remove_dir_all(directory).unwrap();
    }
}
