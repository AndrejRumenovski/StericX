//! Thin JSONL observer for independent population verification. No reference math.
use serde::Deserialize;
use serde_json::json;
use std::io::{self, BufRead};
use steric_x::ensemble::{EnergyMetadata, Population};

#[derive(Deserialize)]
struct Request {
    id: String,
    energies: Vec<EnergyMetadata>,
    temperature_k: f64,
    assumptions: String,
}

fn main() {
    for line in io::stdin().lock().lines() {
        let response = match line {
            Ok(line) => match serde_json::from_str::<Request>(&line) {
                Ok(request) => {
                    let bits = request
                        .energies
                        .iter()
                        .map(|e| e.value.to_bits())
                        .collect::<Vec<_>>();
                    match Population::boltzmann(
                        &request.energies,
                        request.temperature_k,
                        &request.assumptions,
                    ) {
                        Ok(population) => {
                            json!({"id":request.id,"status":"success","population":population,
                            "energy_bits":bits,"temperature_bits":request.temperature_k.to_bits()})
                        }
                        Err(error) => json!({"id":request.id,"status":"error","message":error}),
                    }
                }
                Err(error) => json!({"status":"error","message":error.to_string()}),
            },
            Err(error) => json!({"status":"error","message":error.to_string()}),
        };
        println!("{response}");
    }
}
