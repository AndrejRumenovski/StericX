//! Binding/serialization only: every descriptor, statistic and population calls Rust core.
use pyo3::{create_exception, exceptions::PyValueError, prelude::*};
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use steric_x::{
    density::{DensityDeclaration, DensityGrid},
    ensemble::{
        DescriptorOptions, DescriptorOutcome, EnergyMetadata, Ensemble, Population,
        StatisticsOptions,
    },
};

create_exception!(_native, InputError, PyValueError);
create_exception!(_native, PopulationError, PyValueError);
create_exception!(_native, CalculationError, PyValueError);
create_exception!(_native, DensityError, PyValueError);

fn decode<T: serde::de::DeserializeOwned>(text: &str) -> PyResult<T> {
    serde_json::from_str(text).map_err(|e| InputError::new_err(e.to_string()))
}
fn encode<T: Serialize>(value: &T) -> PyResult<String> {
    // Value widens native f32 to exact f64 before JSON formatting. Direct f32
    // string formatting would lose that identity on Python's f64 JSON decode.
    let value =
        serde_json::to_value(value).map_err(|e| CalculationError::new_err(e.to_string()))?;
    serde_json::to_string(&value).map_err(|e| CalculationError::new_err(e.to_string()))
}
#[derive(Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case", deny_unknown_fields)]
enum PopulationRequest {
    Supplied {
        weights: Vec<f64>,
        provenance: String,
    },
    Boltzmann {
        temperature_k: f64,
        assumptions: String,
    },
}
impl PopulationRequest {
    fn calculate(self, ensemble: &Ensemble) -> PyResult<Population> {
        match self {
            Self::Supplied {
                weights,
                provenance,
            } => Population::supplied(&weights, &provenance),
            Self::Boltzmann {
                temperature_k,
                assumptions,
            } => ensemble.boltzmann_population(temperature_k, &assumptions),
        }
        .map_err(PopulationError::new_err)
    }
}
#[pyclass]
struct EnsembleEngine {
    inner: Ensemble,
}

#[derive(Deserialize)]
struct CoordinateRecord {
    xyz: String,
    identifier: String,
    conformer_id: i64,
}
#[pymethods]
impl EnsembleEngine {
    #[staticmethod]
    fn from_coordinate_records(records_json: &str, bonds_json: &str) -> PyResult<Self> {
        let records: Vec<CoordinateRecord> = decode(records_json)?;
        let bonds: Vec<[usize; 2]> = decode(bonds_json)?;
        if records.is_empty() {
            return Err(InputError::new_err("coordinate records must be nonempty"));
        }
        let mut conformers = Vec::new();
        for (index, record) in records.into_iter().enumerate() {
            let parsed =
                Ensemble::from_xyz_str(&record.xyz, "rdkit").map_err(InputError::new_err)?;
            if parsed.conformers.len() != 1 {
                return Err(InputError::new_err(
                    "each coordinate record must contain one geometry",
                ));
            }
            let mut conformer = parsed.conformers.into_iter().next().unwrap();
            if let Ok(frame) = &mut conformer.frame {
                let mut seen = std::collections::BTreeSet::new();
                for &[a, b] in &bonds {
                    if a == b
                        || a >= frame.molecule.atoms.len()
                        || b >= frame.molecule.atoms.len()
                        || !seen.insert((a.min(b), a.max(b)))
                    {
                        return Err(InputError::new_err(
                            "invalid or duplicate coordinate-record bond",
                        ));
                    }
                }
                frame.bonds = Some(bonds.clone());
            }
            conformer.index = index;
            conformer.identifier = record.identifier;
            conformer.source_format = "rdkit_coordinates".into();
            conformer
                .properties
                .push(steric_x::ensemble::SourceProperty {
                    name: "STERICX_RDKIT_CONFORMER_ID".into(),
                    value: record.conformer_id.to_string(),
                });
            conformers.push(conformer);
        }
        Ok(Self {
            inner: Ensemble { conformers },
        })
    }
    #[staticmethod]
    fn from_file(py: Python<'_>, path: String) -> PyResult<Self> {
        py.detach(|| Ensemble::from_file(path))
            .map(|inner| Self { inner })
            .map_err(InputError::new_err)
    }
    #[staticmethod]
    fn from_text(py: Python<'_>, text: String, format: &str, source: String) -> PyResult<Self> {
        let result = match format {
            "sdf" | "mol" => py.detach(|| Ensemble::from_sdf_str(&text, &source)),
            "xyz" => py.detach(|| Ensemble::from_xyz_str(&text, &source)),
            _ => Err("supported coordinate formats are sdf, mol and xyz".into()),
        };
        result
            .map(|inner| Self { inner })
            .map_err(InputError::new_err)
    }
    fn __len__(&self) -> usize {
        self.inner.conformers.len()
    }
    fn set_energies(&mut self, energies_json: &str) -> PyResult<()> {
        let energies: Vec<Option<EnergyMetadata>> = decode(energies_json)?;
        if energies.len() != self.inner.conformers.len() {
            return Err(PopulationError::new_err(
                "energies must match every ordered conformer",
            ));
        }
        for energy in energies.iter().flatten() {
            energy.validate().map_err(PopulationError::new_err)?;
        }
        for (row, energy) in self.inner.conformers.iter_mut().zip(energies) {
            row.energy = energy;
        }
        Ok(())
    }
    fn records(&self) -> PyResult<String> {
        encode(&self.inner.conformers.iter().map(|c|serde_json::json!({"index":c.index,"identifier":c.identifier,
            "source":c.source,"source_format":c.source_format,"source_geometry":c.source_geometry,"properties":c.properties,
            "energy":c.energy,"geometry_error":c.frame.as_ref().err()})).collect::<Vec<_>>())
    }
    fn descriptors(&self, py: Python<'_>, options_json: &str) -> PyResult<String> {
        let options: DescriptorOptions = decode(options_json)?;
        py.detach(|| encode(&self.inner.descriptors(&options)))
    }
    fn analyze(
        &self,
        py: Python<'_>,
        options_json: &str,
        statistics_json: &str,
        population_json: &str,
        ligand_id: Option<String>,
    ) -> PyResult<String> {
        let options: DescriptorOptions = decode(options_json)?;
        let statistics: StatisticsOptions = decode(statistics_json)?;
        let request: Option<PopulationRequest> = decode(population_json)?;
        let population = request.map(|r| r.calculate(&self.inner)).transpose()?;
        py.detach(|| {
            let descriptors=self.inner.descriptors(&options);
            let stats=descriptors.statistics(&statistics,population.as_ref());
            let sensitivity=ligand_id.as_ref().and_then(|id|stats.as_ref().ok().map(|s|s.sensitivity(id)));
            let csv=sensitivity.as_ref().map(|s|s.to_csv()).transpose().map_err(CalculationError::new_err)?;
            let outcome:DescriptorOutcome<_>=stats.into();
            encode(&serde_json::json!({"descriptors":descriptors,"statistics":outcome,"sensitivity":sensitivity,"sensitivity_csv":csv}))
        })
    }
    fn volume_scan(
        &self,
        py: Python<'_>,
        index: usize,
        radii: Vec<f64>,
        options_json: &str,
    ) -> PyResult<String> {
        let options: DescriptorOptions = decode(options_json)?;
        let conformer = self
            .inner
            .conformers
            .get(index)
            .ok_or_else(|| InputError::new_err("conformer index out of bounds"))?;
        py.detach(|| {
            encode(
                &conformer
                    .buried_volume_scan(&radii, &options)
                    .map_err(CalculationError::new_err)?,
            )
        })
    }
}
#[pyclass]
struct DensityEngine {
    inner: DensityGrid,
}
#[pymethods]
impl DensityEngine {
    #[staticmethod]
    fn from_file(py: Python<'_>, path: String, declaration_json: &str) -> PyResult<Self> {
        let declaration: DensityDeclaration = decode(declaration_json)?;
        py.detach(|| DensityGrid::from_cube_file(path, declaration))
            .map(|inner| Self { inner })
            .map_err(|e| DensityError::new_err(e.to_string()))
    }
    fn metadata(&self) -> PyResult<String> {
        encode(self.inner.metadata())
    }
    fn volume(
        &self,
        py: Python<'_>,
        center: [f64; 3],
        radius: f64,
        threshold: f64,
    ) -> PyResult<String> {
        py.detach(|| {
            encode(
                &self
                    .inner
                    .buried_volume(center, radius, threshold)
                    .map_err(|e| DensityError::new_err(e.to_string()))?,
            )
        })
    }
    fn scan(
        &self,
        py: Python<'_>,
        center: [f64; 3],
        radii: Vec<f64>,
        threshold: f64,
    ) -> PyResult<String> {
        py.detach(|| {
            encode(
                &self
                    .inner
                    .buried_volume_scan(center, &radii, threshold)
                    .map_err(|e| DensityError::new_err(e.to_string()))?,
            )
        })
    }
}
#[pyfunction]
fn batch(
    py: Python<'_>,
    paths: Vec<String>,
    options_json: &str,
    threads: Option<usize>,
) -> PyResult<String> {
    let options: DescriptorOptions = decode(options_json)?;
    if threads == Some(0) {
        return Err(InputError::new_err("threads must be positive"));
    }
    py.detach(|| {
        let work = || {
            paths
                .par_iter()
                .map(|path| {
                    let result: DescriptorOutcome<_> = Ensemble::from_file(path)
                        .map(|e| e.descriptors(&options))
                        .into();
                    serde_json::json!({"path":path,"result":result})
                })
                .collect::<Vec<_>>()
        };
        let values = match threads {
            Some(n) => rayon::ThreadPoolBuilder::new()
                .num_threads(n)
                .build()
                .map_err(|e| InputError::new_err(e.to_string()))?
                .install(work),
            None => work(),
        };
        encode(&values)
    })
}
#[pymodule]
fn _native(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<EnsembleEngine>()?;
    m.add_class::<DensityEngine>()?;
    m.add_function(wrap_pyfunction!(batch, m)?)?;
    m.add("InputError", m.py().get_type::<InputError>())?;
    m.add("PopulationError", m.py().get_type::<PopulationError>())?;
    m.add("CalculationError", m.py().get_type::<CalculationError>())?;
    m.add("DensityError", m.py().get_type::<DensityError>())?;
    m.add("__version__", env!("CARGO_PKG_VERSION"))?;
    Ok(())
}
