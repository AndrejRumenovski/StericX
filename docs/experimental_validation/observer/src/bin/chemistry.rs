//! Research adapter: authoritative Rust descriptors, ensemble statistics and weights.
use serde_json::{Value, json};
use steric_x::ensemble::{
    DescriptorOptions, EnergyMetadata, Ensemble, Population, StatisticsOptions, SterimolAxis,
};
use steric_x::geometry::{
    BuriedVolumeCalculator, CertifiedVolumeOptions, SterimolCalculator,
    coordination_center_with_neighbors, detect_trivalent_donor,
};
fn main() -> Result<(), Box<dyn std::error::Error>> {
    let requests: Vec<Value> = serde_json::from_str(&std::fs::read_to_string(
        std::env::args().nth(1).ok_or("request file")?,
    )?)?;
    for q in requests {
        let run = || -> Result<Value, Box<dyn std::error::Error>> {
            let mut ensemble = Ensemble {
                conformers: Vec::new(),
            };
            for f in q["files"].as_array().ok_or("files")? {
                let one = Ensemble::from_file(f["path"].as_str().ok_or("path")?)?;
                if one.conformers.len() != 1 {
                    return Err("one original conformer per file required".into());
                }
                let mut c = one.conformers.into_iter().next().unwrap();
                c.index = ensemble.conformers.len();
                c.identifier = f["identifier"].as_str().ok_or("identifier")?.into();
                if !f["energy"].is_null() {
                    c.energy = Some(serde_json::from_value::<EnergyMetadata>(
                        f["energy"].clone(),
                    )?);
                }
                ensemble.conformers.push(c);
            }
            let opts = DescriptorOptions {
                axis: SterimolAxis::Coordination,
                ..Default::default()
            };
            let descriptors = ensemble.descriptors(&opts);
            let statistics = descriptors.statistics(&StatisticsOptions::default(), None)?;
            let mut populations = Vec::new();
            if ensemble.conformers.iter().all(|c| c.energy.is_some()) {
                let energies: Vec<_> = ensemble
                    .conformers
                    .iter()
                    .filter_map(|c| c.energy.clone())
                    .collect();
                for t in [298.15, 353.15, 400.] {
                    let population = Population::boltzmann(
                        &energies,
                        t,
                        "Approximate force-field/electronic-energy weighting, not measured populations. Enumerated retained conformers treated as distinct states with degeneracy one; incomplete search and missing entropy/solvent are limitations.",
                    )?;
                    populations.push(json!({"temperature_k":t,"statistics":descriptors.statistics(&StatisticsOptions::default(),Some(&population))?}));
                }
            }
            let mut precision = Vec::new();
            for c in &ensemble.conformers {
                let result = || -> Result<Value, Box<dyn std::error::Error>> {
                    let frame = c.frame.as_ref().map_err(|s| s.as_str())?;
                    let topology = detect_trivalent_donor(
                        &frame.molecule,
                        frame.bonds.as_deref(),
                        "P",
                        None,
                        None,
                        true,
                    )?;
                    let center = coordination_center_with_neighbors(
                        &frame.molecule,
                        topology.donor_idx,
                        topology.substituents,
                        opts.volume,
                    )?;
                    let continuous = SterimolCalculator::compute_continuous_with_dummy(
                        &frame.molecule,
                        topology.donor_idx,
                        center,
                    )?;
                    let mut regional = Vec::new();
                    if q["certified"].as_bool().unwrap_or(false) {
                        for plane in topology.substituents {
                            let prepared = BuriedVolumeCalculator::prepare_certified_frame(
                                &frame.molecule,
                                topology.donor_idx,
                                plane,
                                center,
                                f64::from(opts.volume.radii_scale),
                                opts.volume.include_hydrogens,
                            )?;
                            let volumes = prepared.compute_regions(
                                f64::from(opts.volume.sphere_radius),
                                CertifiedVolumeOptions {
                                    absolute_tolerance: 0.05,
                                    max_cells: 1048576,
                                },
                            )?;
                            regional.push(json!({"plane_atom_index":plane,"prepared_frame":prepared,"volumes":volumes}));
                        }
                    }
                    Ok(
                        json!({"continuous_sterimol_raw_l":continuous,"center":center.to_array(),"donor":topology.donor_idx,"neighbors":topology.substituents,"certified_regions":regional}),
                    )
                }();
                precision.push(match result {
                    Ok(v) => json!({"identifier":c.identifier,"result":v}),
                    Err(e) => json!({"identifier":c.identifier,"error":e.to_string()}),
                });
            }
            Ok(
                json!({"descriptors":descriptors,"sensitivity":statistics.sensitivity(q["id"].as_str().unwrap_or("unknown")),"population_analyses":populations,"precision":precision}),
            )
        };
        println!(
            "{}",
            match run() {
                Ok(v) =>
                    json!({"id":q["id"],"ligand_id":q["ligand_id"],"lane":q["lane"],"result":v}),
                Err(e) => json!({"id":q["id"],"error":e.to_string()}),
            }
        );
    }
    Ok(())
}
